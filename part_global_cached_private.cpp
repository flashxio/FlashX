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
#include "parameters.h"

const int MIN_NUM_PROCESS_REQ = 100;
const int BUF_SIZE = 100;

#if 0
thread_group::~thread_group()
{
	for (size_t i = 0; i < process_request_threads.size(); i++) {
		thread *t = process_request_threads[i];
		t->stop();
		// TODO I need to free memory of the thread here.
	}
	for (size_t i = 0; i < process_reply_threads.size(); i++) {
		thread *t = process_reply_threads[i];
		t->stop();
		// TODO I need to free memory of the thread here.
	}
	delete request_queue;
	delete reply_queue;
}
#endif

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
	io_reply rep(req, true, 0);
	io->reply(req, &rep, 1);
	return 0;
}

/**
 * This thread runs on the node with cache.
 * It is dedicated to processing requests.
 * The application threads send requests to the node with cache
 * by placing the requests in the queue. This thread fetches requests
 * from the queue and process them. If the queue is empty, the thread
 * will sleep.
 */
class process_request_thread: public thread
{
	part_global_cached_io *io;
public:
	process_request_thread(part_global_cached_io *io): thread(
			std::string("process_request_thread-") + itoa(io->get_node_id()),
			// We don't use blocking mode of the thread because
			// we are using blocking queue.
			io->get_node_id(), false) {
		this->io = io;
	}
	void run();
};

/**
 * This thread runs on the node with the application threads.
 * It is dedicated to processing replies. If the queue is empty,
 * the thread will sleep.
 */
class process_reply_thread: public thread
{
	part_global_cached_io *io;
public:
	process_reply_thread(part_global_cached_io *io): thread(
			std::string("process_reply_thread-") + itoa(io->get_node_id()),
			io->get_node_id(), false) {
		this->io = io;
	}
	void run();
};

void process_request_thread::run()
{
	io->process_requests(NUMA_NUM_PROCESS_MSGS);
}

void process_reply_thread::run()
{
	io->process_replies(NUMA_NUM_PROCESS_MSGS);
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

	/* 
	 * there is a global lock for all threads.
	 * so this lock makes sure cache initialization is serialized
	 */
	pthread_mutex_lock(&init_mutex);
	thread_group *group = local_group;
	if (group->cache == NULL) {
		/* Each cache has their own memory managers */
		group->cache = global_cached_io::create_cache(cache_type, cache_size,
				// TODO I need to set the node id right.
				-1, 1);
		group->request_queue = new blocking_FIFO_queue<io_request>("request_queue", 
				NUMA_REQ_QUEUE_SIZE);
		group->reply_queue = new blocking_FIFO_queue<io_reply>("reply_queue", 
				NUMA_REPLY_QUEUE_SIZE);
		// Create processing threads.
		for (int i = 0; i < NUMA_NUM_PROCESS_THREADS; i++) {
			thread *t = new process_request_thread(this);
			t->start();
			printf("create a request processing thread, blocking: %d\n", t->is_blocking());
			group->process_request_threads.push_back(t);
			t = new process_reply_thread(this);
			t->start();
			printf("create a reply processing thread, blocking: %d\n", t->is_blocking());
			group->process_reply_threads.push_back(t);
		}
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
	for (std::tr1::unordered_map<int, struct thread_group>::const_iterator it
			= groups.begin(); it != groups.end(); it++) {
		thread_safe_FIFO_queue<io_request> *q[1];
		q[0] = it->second.request_queue;
		req_senders.insert(std::pair<int, msg_sender<io_request> *>(it->first,
					new msg_sender<io_request>(NUMA_MSG_CACHE_SIZE, q, 1)));
		thread_safe_FIFO_queue<io_reply> *q1[1];
		q1[0] = it->second.reply_queue;
		reply_senders.insert(std::pair<int, msg_sender<io_reply> *>(it->first,
					new msg_sender<io_reply>(NUMA_MSG_CACHE_SIZE, q1, 1)));
	}

	return 0;
}

part_global_cached_io::part_global_cached_io(int num_groups,
		io_interface *underlying, int idx, long cache_size,
		int cache_type, access_mapper *mapper): global_cached_io(underlying) {
	this->mapper = mapper->clone();
	this->thread_id = idx;
	this->final_cb = NULL;
	this->my_cb = NULL;
	remote_reads = 0;
	//		assert(nthreads % num_groups == 0);
	this->num_groups = num_groups;
	this->group_idx = underlying->get_node_id();
	this->cache_size = (long) (cache_size * ((double) underlying->get_local_size()
		/ underlying->get_size()));
	this->cache_type = cache_type;
	processed_requests = 0;
	finished_threads = 0;

	printf("cache is partitioned\n");
	printf("thread id: %d, group id: %d, num groups: %d, cache size: %ld\n",
			idx, group_idx, num_groups, this->cache_size);

	if (groups.size() == 0) {
		pthread_mutex_init(&init_mutex, NULL);
		pthread_mutex_init(&wait_mutex, NULL);
		pthread_cond_init(&cond, NULL);
		num_finish_init = 0;
	}
	// If the group hasn't been added to the table yet.
	if (groups.find(group_idx) == groups.end()) {
		struct thread_group group;
		group.id = group_idx;
		group.nthreads = nthreads / num_groups;
		if (nthreads % num_groups)
			group.nthreads++;
		group.ios = new part_global_cached_io*[group.nthreads];
		group.cache = NULL;
		for (int j = 0; j < group.nthreads; j++)
			group.ios[j] = NULL;
		groups.insert(std::pair<int, struct thread_group>(group_idx,
					group));
	}

	/* assign a thread to a group. */
	local_group = &groups[group_idx];
	int i = thread_idx(idx, num_groups);
	assert (local_group->ios[i] == NULL);
	local_group->ios[i] = this;

	my_cb = new for_global_callback(this);
	global_cached_io::set_callback(my_cb);
}

/**
 * send replies to the thread that sent the requests.
 */
int part_global_cached_io::reply(io_request *requests,
		io_reply *replies, int num) {
	for (int i = 0; i < num; i++) {
		part_global_cached_io *io = (part_global_cached_io *) requests[i].get_io();
		int group_id = io->group_idx;
		int num_sent = reply_senders[group_id]->send_cached(&replies[i]);
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
		if (idx != get_group_id())
			remote_reads++;
		int ret = req_senders[idx]->send_cached(&requests[i]);
		if (ret == 0)
			break;
		num_sent++;
	}
	for (std::tr1::unordered_map<int, msg_sender<io_request> *>::const_iterator it
			= req_senders.begin(); it != req_senders.end(); it++)
		it->second->flush();
	return num_sent;
}

/* process the requests sent to this thread */
int part_global_cached_io::process_requests(int max_nreqs) {
	int num_processed = 0;
	thread_safe_FIFO_queue<io_request> *request_queue
		= local_group->request_queue;
	io_request local_reqs[BUF_SIZE];
	while (num_processed < max_nreqs) {
		int num = request_queue->fetch(local_reqs, BUF_SIZE);
		for (int i = 0; i < num; i++) {
			io_request *req = &local_reqs[i];
			assert(req->get_offset() >= 0);
			assert(req->get_size() > 0);
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
	thread_safe_FIFO_queue<io_reply> *reply_queue = local_group->reply_queue;
	io_reply local_replies[BUF_SIZE];
	while(num_processed < max_nreplies) {
		int num = reply_queue->fetch(local_replies, BUF_SIZE);
		for (int i = 0; i < num; i++) {
			io_reply *reply = &local_replies[i];
			process_reply(reply);
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
	final_cb->invoke(&req);
	return ret;
}

ssize_t part_global_cached_io::access(io_request *requests, int num) {
	int num_sent = 0;
	// If we can't sent the requests to the destination node,
	// we have to busy wait.
	while (num - num_sent > 0) {
		int ret = distribute_reqs(&requests[num_sent], num - num_sent);
		if (ret > 0)
			num_sent += ret;
	}
	return num_sent;
}

void part_global_cached_io::cleanup() {
	printf("thread %d: start to clean up\n", thread_id);
	for (std::tr1::unordered_map<int, struct thread_group>::const_iterator it
			= groups.begin(); it != groups.end(); it++) {
		for (int j = 0; j < it->second.nthreads; j++)
			if (it->second.ios[j])
				__sync_fetch_and_add(&it->second.ios[j]->finished_threads, 1);
	}
	while (!local_group->request_queue->is_empty()
			|| !local_group->reply_queue->is_empty()
			/*
			 * if finished_threads == nthreads,
			 * then all threads have reached the point.
			 */
			|| finished_threads < nthreads) {
	}
	printf("thread %d processed %ld requests\n", thread_id, processed_requests);
}

std::tr1::unordered_map<int, thread_group> part_global_cached_io::groups;
pthread_mutex_t part_global_cached_io::init_mutex;
int part_global_cached_io::num_finish_init;
pthread_mutex_t part_global_cached_io::wait_mutex;
pthread_cond_t part_global_cached_io::cond;
