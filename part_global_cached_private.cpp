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

int part_global_cached_private::thread_init() {
	/* let's bind the thread to a specific node first. */
	struct bitmask *nodemask = numa_allocate_cpumask();
	numa_bitmask_clearall(nodemask);
	printf("thread %d is associated to node %d\n", idx, get_group_id());
	numa_bitmask_setbit(nodemask, get_group_id());
	numa_bind(nodemask);
	numa_set_strict(1);
	numa_set_bind_policy(1);

	read_private::thread_init();
	request_queue = new bulk_queue<io_request>(REQ_QUEUE_SIZE);
	reply_queue = new bulk_queue<io_reply>(REPLY_QUEUE_SIZE);

	/* 
	 * there is a global lock for all threads.
	 * so this lock makes sure cache initialization is serialized
	 */
	pthread_mutex_lock(&init_mutex);
	thread_group *group = &groups[group_idx];
	if (group->cache == NULL) {
		/* this allocates all pages for the cache. */
		page::allocate_cache(cache_size);
		group->cache = global_cached_private::create_cache(cache_type, cache_size, manager);
	}
	num_finish_init++;
	pthread_mutex_unlock(&init_mutex);

	pthread_mutex_lock(&wait_mutex);
	while (num_finish_init < nthreads) {
		pthread_cond_wait(&cond, &wait_mutex);
	}
	pthread_mutex_unlock(&wait_mutex);
	pthread_mutex_lock(&wait_mutex);
	pthread_cond_broadcast(&cond);
	pthread_mutex_unlock(&wait_mutex);
	printf("thread %d finishes initialization\n", idx);

	/*
	 * we have to initialize senders here
	 * because we need to make sure other threads have 
	 * initialized all queues.
	 */
	/* there is a request sender for each node. */
	req_senders = (msg_sender<io_request> **) numa_alloc_local(
			sizeof(msg_sender<io_request> *) * num_groups);
	for (int i = 0; i < num_groups; i++) {
		bulk_queue<io_request> *queues[groups[i].nthreads];
		for (int j = 0; j < groups[i].nthreads; j++)
			queues[j] = groups[i].threads[j]->request_queue;
		req_senders[i] = new msg_sender<io_request>(BUF_SIZE, queues, groups[i].nthreads);
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
			bulk_queue<io_reply> *queues[1];
			assert(idx == groups[i].threads[j]->idx);
			queues[0] = groups[i].threads[j]->reply_queue;
			reply_senders[idx++] = new msg_sender<io_reply>(BUF_SIZE, queues, 1);
		}
	}

	return 0;
}

part_global_cached_private::part_global_cached_private(int num_groups,
		const char *names[], int num, long size, int idx,
		long cache_size, int entry_size, int cache_type,
		memory_manager *manager): global_cached_private(names,
			num, size, idx, entry_size) {
	this->manager = manager;
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
			groups[i].threads = new part_global_cached_private*[groups[i].nthreads];
			groups[i].cache = NULL;
			for (int j = 0; j < groups[i].nthreads; j++)
				groups[i].threads[j] = NULL;
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
	assert (group->threads[i] == NULL);
	group->threads[i] = this;
}

/**
 * send replies to the thread that sent the requests.
 */
int part_global_cached_private::reply(io_request *requests,
		io_reply *replies, int num) {
	for (int i = 0; i < num; i++) {
		part_global_cached_private *thread
			= (part_global_cached_private *) requests[i].get_thread();
		int thread_id = thread->idx;
		int num_sent = reply_senders[thread_id]->send_cached(&replies[i]);
		if (num_sent == 0) {
			// TODO the buffer is already full.
			// discard the request for now.
			printf("the reply buffer for thread %d is already full\n", thread_id);
			continue;
		}
	}
	for (int i = 0; i < nthreads; i++)
		reply_senders[i]->flush();
	return 0;
}

/* distribute requests to nodes. */
void part_global_cached_private::distribute_reqs(io_request *requests, int num) {
	for (int i = 0; i < num; i++) {
		int idx = hash_req(&requests[i]);
		assert (idx < num_groups);
		if (idx != get_group_id())
			remote_reads++;
		int num_sent = req_senders[idx]->send_cached(&requests[i]);
		// TODO if we fail to send the requests to the specific node,
		// we should rehash it and give it to another node.
		if (num_sent == 0) {
			// TODO the buffer is already full.
			// discard the request for now.
			printf("the request buffer for group %d is already full\n", idx);
			continue;
		}
	}
	for (int i = 0; i < num_groups; i++)
		req_senders[i]->flush();
}

/* process the requests sent to this thread */
int part_global_cached_private::process_requests(int max_nreqs) {
	int num_processed = 0;
	io_request local_reqs[BUF_SIZE];
	io_reply local_replies[BUF_SIZE];
	while (!request_queue->is_empty() && num_processed < max_nreqs) {
		int num = request_queue->fetch(local_reqs, BUF_SIZE);
		for (int i = 0; i < num; i++) {
			io_request *req = &local_reqs[i];
			// TODO will it be better if I collect all data
			// and send them back to the initiator in blocks?
			int access_method = req->get_access_method();
			assert(req->get_offset() >= 0);
			int ret = global_cached_private::access(req->get_buf(),
					req->get_offset(), req->get_size(), access_method);
			local_replies[i] = io_reply(req, ret >= 0, errno);
		}
		num_processed += num;
		reply(local_reqs, local_replies, num);
	}
	processed_requests += num_processed;
	return num_processed;
}

/* 
 * process the replies and return the number
 * of bytes that have been accessed.
 */
int part_global_cached_private::process_replies(int max_nreplies) {
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

int part_global_cached_private::process_reply(io_reply *reply) {
	int access_method = reply->get_access_method();
	int ret = -1;
	if (reply->is_success()) {
		read_bytes += reply->get_size();
		if (access_method == READ) {
			assert(*(unsigned long *) reply->get_buf()
					== reply->get_offset() / sizeof(long));
		}
		ret = reply->get_size();
		buf->free_entry(reply->get_buf());
	}
	else {
		fprintf(stderr, "access error: %s\n",
				strerror(reply->get_status()));
	}
	return ret;
}

ssize_t part_global_cached_private::access(io_request *requests,
		int num, int access_method) {
	distribute_reqs(requests, num);
	/*
	 * let's process up to twice as many requests as demanded by the upper layer,
	 * I hope this can help load balancing problem.
	 */
	int num_recv = 0;
	/* we need to process them at least once. */
	process_requests(num * 2);
	num_recv += process_replies(num * 4);
	while (buf->is_full()) {
		process_requests(num * 2);
		num_recv += process_replies(num * 4);
	}
	return num_recv;
}

void part_global_cached_private::cleanup() {
	int num = 0;
	printf("thread %d: start to clean up\n", idx);
	for (int i = 0; i < num_groups; i++) {
		for (int j = 0; j < groups[i].nthreads; j++)
			if (groups[i].threads[j])
				__sync_fetch_and_add(&groups[i].threads[j]->finished_threads, 1);
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
	printf("thread %d processed %ld requests\n", idx, processed_requests);
}

thread_group *part_global_cached_private::groups;
pthread_mutex_t part_global_cached_private::init_mutex;
int part_global_cached_private::num_finish_init;
pthread_mutex_t part_global_cached_private::wait_mutex;
pthread_cond_t part_global_cached_private::cond;
