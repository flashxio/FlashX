#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <pthread.h>
#include <assert.h>

#include "cache.h"
#include "associative_cache.h"
#include "global_cached_private.h"
#include "part_global_cached_private.h"

const int MIN_NUM_PROCESS_REQ = 100;
const int BUF_SIZE = 100;

/**
 * This callback is used in global cache IO.
 * When a request is complete in the global cache IO, this callback
 * will be invoked.
 */
class for_global_callback: public callback
{
	part_global_cached_io *io;
public:
	for_global_callback(part_global_cached_io *io) {
		this->io = io;
	}

	int invoke(io_request *request);
};

int for_global_callback::invoke(io_request *req)
{
	printf("a request is complete\n");
	io_reply rep(req, true, 0);
	io->reply(req, &rep, 1);
	return 0;
}

bool part_global_cached_io::set_callback(callback *cb)
{
	this->cb = cb;
	return true;
}

int part_global_cached_io::init() {
#if NUM_NODES > 1
	/* let's bind the thread to a specific node first. */
	struct bitmask *nodemask = numa_allocate_cpumask();
	numa_bitmask_clearall(nodemask);
	printf("thread %d is associated to node %d\n", thread_id, get_group_id());
	numa_bitmask_setbit(nodemask, get_group_id());
	numa_bind(nodemask);
	numa_set_strict(1);
	numa_set_bind_policy(1);
#endif

	request_queue = new thread_safe_FIFO_queue<io_request>(NUMA_REQ_QUEUE_SIZE);
	reply_queue = new thread_safe_FIFO_queue<io_reply>(NUMA_REPLY_QUEUE_SIZE);

	/* 
	 * there is a global lock for all threads.
	 * so this lock makes sure cache initialization is serialized
	 */
	pthread_mutex_lock(&init_mutex);
	thread_group *group = &groups[group_idx];
	if (group->cache == NULL) {
		/* Each cache has their own memory managers */
		group->cache = global_cached_io::create_cache(cache_type, cache_size,
				// TODO I need to set the node id right.
				-1, num_groups);
	}
	num_finish_init++;
	pthread_mutex_unlock(&init_mutex);
	global_cached_io::init();

	pthread_mutex_lock(&wait_mutex);
	while (num_finish_init < nthreads) {
		pthread_cond_wait(&cond, &wait_mutex);
	}
	pthread_mutex_unlock(&wait_mutex);
	pthread_mutex_lock(&wait_mutex);
	pthread_cond_broadcast(&cond);
	pthread_mutex_unlock(&wait_mutex);
	printf("thread %d finishes initialization\n", thread_id);

	/*
	 * we have to initialize senders here
	 * because we need to make sure other threads have 
	 * initialized all queues.
	 */
	/* there is a request sender for each node. */
	req_senders = (msg_sender<io_request> **) numa_alloc_local(
			sizeof(msg_sender<io_request> *) * num_groups);
	for (int i = 0; i < num_groups; i++) {
		thread_safe_FIFO_queue<io_request> *queues[groups[i].nthreads];
		for (int j = 0; j < groups[i].nthreads; j++)
			queues[j] = groups[i].ios[j]->request_queue;
		req_senders[i] = new msg_sender<io_request>(NUMA_MSG_CACHE_SIZE,
				queues, groups[i].nthreads);
	}
	/* 
	 * there is a reply sender for each thread.
	 * therefore, there is only one queue for a sender.
	 */
	reply_senders = (msg_sender<io_reply> **) numa_alloc_local(
			sizeof(msg_sender<io_reply> *) * nthreads);
	int idx = 0;
	for (int i = 0; i < num_groups; i++) {
		for (int j = 0; j < groups[i].nthreads; j++) {
			thread_safe_FIFO_queue<io_reply> *queues[1];
			assert(idx == groups[i].ios[j]->thread_id);
			queues[0] = groups[i].ios[j]->reply_queue;
			reply_senders[idx++] = new msg_sender<io_reply>(NUMA_MSG_CACHE_SIZE,
					queues, 1);
		}
	}

	return 0;
}

part_global_cached_io::part_global_cached_io(int num_groups,
		io_interface *underlying, int idx, long cache_size,
		int cache_type): global_cached_io(underlying) {
	this->thread_id = idx;
	this->cb = NULL;
	remote_reads = 0;
	//		assert(nthreads % num_groups == 0);
	this->num_groups = num_groups;
	this->group_idx = group_id(idx, num_groups);
	this->cache_size = cache_size / num_groups;
	this->cache_type = cache_type;
	processed_requests = 0;
	finished_threads = 0;
	req_senders = NULL;
	reply_senders = NULL;

	printf("cache is partitioned\n");
	printf("thread id: %d, group id: %d, num groups: %d\n", idx, group_idx, num_groups);

	if (groups == NULL) {
		pthread_mutex_init(&init_mutex, NULL);
		pthread_mutex_init(&wait_mutex, NULL);
		pthread_cond_init(&cond, NULL);
		num_finish_init = 0;
		groups = new thread_group[num_groups];
		for (int i = 0; i < num_groups; i++) {
			groups[i].id = i;
			groups[i].nthreads = nthreads / num_groups;
			if (nthreads % num_groups)
				groups[i].nthreads++;
			groups[i].ios = new part_global_cached_io*[groups[i].nthreads];
			groups[i].cache = NULL;
			for (int j = 0; j < groups[i].nthreads; j++)
				groups[i].ios[j] = NULL;
		}
	}

	/* assign a thread to a group. */
	thread_group *group = NULL;
	for (int i = 0; i < num_groups; i++) {
		if (groups[i].id == group_idx) {
			group = &groups[i];
			break;
		}
	}
	assert (group);
	int i = thread_idx(idx, num_groups);
	assert (group->ios[i] == NULL);
	group->ios[i] = this;

	global_cached_io::set_callback(new for_global_callback(this));
}

/**
 * send replies to the thread that sent the requests.
 */
int part_global_cached_io::reply(io_request *requests,
		io_reply *replies, int num) {
	for (int i = 0; i < num; i++) {
		part_global_cached_io *io = (part_global_cached_io *) requests[i].get_io();
		int thread_id = io->thread_id;
		int num_sent = reply_senders[thread_id]->send_cached(&replies[i]);
		if (num_sent == 0) {
			// TODO the buffer is already full.
			// discard the request for now.
			printf("the reply buffer for thread %d is already full\n", thread_id);
			continue;
		}
	}
	// TODO We shouldn't flush here, but we need to flush at some point.
#if 0
	for (int i = 0; i < nthreads; i++)
		reply_senders[i]->flush();
#endif
	return 0;
}

/* distribute requests to nodes. */
int part_global_cached_io::distribute_reqs(io_request *requests, int num) {
	int num_sent = 0;
	for (int i = 0; i < num; i++) {
		int idx = hash_req(&requests[i]);
		assert (idx < num_groups);
		if (idx != get_group_id())
			remote_reads++;
		int ret = req_senders[idx]->send_cached(&requests[i]);
		if (ret == 0)
			break;
		num_sent++;
	}
	for (int i = 0; i < num_groups; i++)
		req_senders[i]->flush();
	return num_sent;
}

/* process the requests sent to this thread */
int part_global_cached_io::process_requests(int max_nreqs) {
	int num_processed = 0;
	io_request local_reqs[BUF_SIZE];
	while (!request_queue->is_empty() && num_processed < max_nreqs) {
		int num = request_queue->fetch(local_reqs, BUF_SIZE);
		for (int i = 0; i < num; i++) {
			io_request *req = &local_reqs[i];
			assert(req->get_offset() >= 0);
			// The thread will be blocked by the global cache IO
			// if too many requests flush in.
			int ret = global_cached_io::access(req, 1);
			if (ret < 0) {
				fprintf(stderr, "part global cache can't issue a request\n");
				io_reply rep(req, false, errno);
				reply(req, &rep, 1);
			}
		}
		num_processed += num;
	}
	processed_requests += num_processed;
	return num_processed;
}

/* 
 * process the replies and return the number
 * of bytes that have been accessed.
 */
int part_global_cached_io::process_replies(int max_nreplies) {
	int num_processed = 0;
	io_reply local_replies[BUF_SIZE];
	int size = 0;
	while(!reply_queue->is_empty() && num_processed < max_nreplies) {
		int num = reply_queue->fetch(local_replies, BUF_SIZE);
		for (int i = 0; i < num; i++) {
			io_reply *reply = &local_replies[i];
			int ret = process_reply(reply);
			if (ret >= 0)
				size += ret;
		}
		num_processed += num;
	}
	return num_processed;
}

int part_global_cached_io::process_reply(io_reply *reply) {
	int ret = -1;
	if (reply->is_success()) {
		ret = reply->get_size();
	}
	else {
		fprintf(stderr, "access error: %s\n",
				strerror(reply->get_status()));
	}
	io_request req(reply->get_buf(), reply->get_offset(), reply->get_size(),
			// It doesn't really matter what node id is specified
			// for the request. The request is just used for notifying
			// the user code of the completion of the request.
			reply->get_access_method(), get_thread(thread_id)->get_io(), -1);
	cb->invoke(&req);
	return ret;
}

ssize_t part_global_cached_io::access(io_request *requests, int num) {
	int num_recv = 0;
	int num_sent = 0;
	// If we can't sent the requests to the destination node,
	// we have to busy wait.
	while (num - num_sent > 0) {
		int ret = distribute_reqs(&requests[num_sent], num - num_sent);
		if (ret > 0)
			num_sent += ret;
		/*
		 * let's process up to twice as many requests as demanded by the upper layer,
		 * I hope this can help load balancing problem.
		 */
		if (num == 0)
			num = MIN_NUM_PROCESS_REQ;
		/* we need to process them at least once. */
		process_requests(num * 2);
		ret = process_replies(num * 4);
		num_recv += ret;
	}
	return num_recv;
}

void part_global_cached_io::cleanup() {
	int num = 0;
	printf("thread %d: start to clean up\n", thread_id);
	for (int i = 0; i < num_groups; i++) {
		for (int j = 0; j < groups[i].nthreads; j++)
			if (groups[i].ios[j])
				__sync_fetch_and_add(&groups[i].ios[j]->finished_threads, 1);
	}
	while (!request_queue->is_empty()
			|| !reply_queue->is_empty()
			/*
			 * if finished_threads == nthreads,
			 * then all threads have reached the point.
			 */
			|| finished_threads < nthreads) {
		process_requests(200);
		process_replies(200);
		num++;
	}
	printf("thread %d processed %ld requests\n", thread_id, processed_requests);
}

thread_group *part_global_cached_io::groups;
pthread_mutex_t part_global_cached_io::init_mutex;
int part_global_cached_io::num_finish_init;
pthread_mutex_t part_global_cached_io::wait_mutex;
pthread_cond_t part_global_cached_io::cond;
