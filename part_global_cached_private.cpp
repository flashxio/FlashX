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
 * This wraps global cached IO to handle requests in the cache
 * of the local node.
 */
class node_cached_io: public global_cached_io
{
	// thread id <-> msg sender
	pthread_key_t replier_key;
	const std::tr1::unordered_map<int, struct thread_group> *groups;
	const struct thread_group *local_group;
	long processed_requests;
	long num_requests;

	pthread_t processing_thread_id;

	const struct thread_group *get_group(int group_id) {
		std::tr1::unordered_map<int, struct thread_group>::const_iterator it
			= groups->find(get_node_id());
		assert(it != groups->end());
		return &it->second;
	}
public:
	node_cached_io(io_interface *underlying, const std::tr1::unordered_map<int,
			struct thread_group> *groups);

	virtual page_cache *get_global_cache() {
		return local_group->cache;
	}

	long get_num_requests() const {
		return num_requests;
	}

	std::tr1::unordered_map<int, msg_sender<io_reply> *> * init_repliers();

	int process_requests(int max_nreqs);

	int reply(io_request *requests, io_reply *replies, int num);

	void flush_replies();

	virtual void notify_completion(io_request *req);
};

node_cached_io::node_cached_io(io_interface *underlying,
		const std::tr1::unordered_map<int, struct thread_group> *groups):
	global_cached_io(underlying->clone())
{
	this->groups = groups;
	local_group = get_group(underlying->get_node_id());
	assert(local_group);
	processed_requests = 0;
	// This IO instance is created inside the right thread.
	global_cached_io::init();

	processing_thread_id = 0;
	pthread_key_create(&replier_key, NULL);
}

std::tr1::unordered_map<int, msg_sender<io_reply> *> *node_cached_io::init_repliers()
{
	std::tr1::unordered_map<int, msg_sender<io_reply> *> *reply_senders
		= new std::tr1::unordered_map<int, msg_sender<io_reply> *>();
	for (std::tr1::unordered_map<int, struct thread_group>::const_iterator it
			= groups->begin(); it != groups->end(); it++) {
		const struct thread_group *group = &it->second;
		for (int j = 0; j < group->nthreads; j++) {
			// We need a sender for each parted global IO because we need to
			// make sure replies are sent to the right the IO instance that
			// issued the request.
			fifo_queue<io_reply> *q1[1];
			q1[0] = group->ios[j]->reply_queue;
			reply_senders->insert(std::pair<int, msg_sender<io_reply> *>(
						group->ios[j]->thread_id,
						new msg_sender<io_reply>(NUMA_MSG_CACHE_SIZE, q1, 1)));
		}
	}
	pthread_setspecific(replier_key, reply_senders);
	return reply_senders;
}

void node_cached_io::flush_replies()
{
	// TODO we need a way to flush replies.
#if 0
	for (std::tr1::unordered_map<int, msg_sender<io_reply> *>::const_iterator it
			= reply_senders.begin(); it != reply_senders.end(); it++) {
		msg_sender<io_reply> *sender = it->second;
		sender->flush_all();
	}
#endif
}

/* process the requests sent to this thread */
int node_cached_io::process_requests(int max_nreqs)
{
	if (processing_thread_id == 0)
		processing_thread_id = pthread_self();
	assert(processing_thread_id == pthread_self());
	int num_processed = 0;
	fifo_queue<io_request> *request_queue = local_group->request_queue;
	io_request local_reqs[NUMA_REQ_BUF_SIZE];
	while (num_processed < max_nreqs) {
		int num = request_queue->fetch(local_reqs, NUMA_REQ_BUF_SIZE);
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
	num_requests += num_processed;
	return num_processed;
}

/**
 * send replies to the thread that sent the requests.
 */
int node_cached_io::reply(io_request *requests, io_reply *replies, int num)
{
	std::tr1::unordered_map<int, msg_sender<io_reply> *> *reply_senders
		= (std::tr1::unordered_map<int, msg_sender<io_reply> *> *) pthread_getspecific(replier_key);
	if (reply_senders == NULL)
		reply_senders = init_repliers();
	for (int i = 0; i < num; i++) {
		part_global_cached_io *io = (part_global_cached_io *) requests[i].get_io();
		// If the reply is sent to the thread on a different node.
		if (io->get_node_id() != this->get_node_id()) {
			int num_sent = (*reply_senders)[io->thread_id]->send_cached(&replies[i]);
			// We use blocking queues here, so the send must succeed.
			assert(num_sent > 0);
		}
		else {
			io->process_reply(&replies[i]);
			io->processed_replies.inc(1);
		}
	}
	return 0;
}

void node_cached_io::notify_completion(io_request *req)
{
	io_reply rep(req, true, 0);
	reply(req, &rep, 1);
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
	node_cached_io *io;
public:
	process_request_thread(node_cached_io *io): thread(
			std::string("process_request_thread-") + itoa(io->get_node_id()),
			// We don't use blocking mode of the thread because
			// we are using blocking queue.
			io->get_node_id(), false) {
		this->io = io;
	}
	node_cached_io *get_io() const {
		return io;
	}
	void run();

	int get_cache_hits() const {
		return io->get_cache_hits();
	}

	long get_num_requests() const {
		return io->get_num_requests();
	}
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
	std::vector<node_cached_io *> node_ios;
	pthread_mutex_lock(&init_mutex);
	thread_group *group = &groups[group_idx];
	if (group->cache == NULL) {
		/* Each cache has their own memory managers */
		group->cache = global_cached_io::create_cache(cache_type, cache_size,
				group_idx, 1);
		group->request_queue = new blocking_FIFO_queue<io_request>("request_queue", 
				NUMA_REQ_QUEUE_SIZE);
		// Create processing threads.
		for (int i = 0; i < NUMA_NUM_PROCESS_THREADS; i++) {
			node_cached_io *io = new node_cached_io(underlying->clone(), &groups);
			node_ios.push_back(io);
			thread *t = new process_request_thread(io);
			t->start();
			printf("create a request processing thread, blocking: %d\n", t->is_blocking());
			group->process_request_threads.push_back(t);
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
		const struct thread_group *group = &it->second;
		fifo_queue<io_request> *q[1];
		q[0] = group->request_queue;
		req_senders.insert(std::pair<int, msg_sender<io_request> *>(it->first,
					new msg_sender<io_request>(NUMA_MSG_CACHE_SIZE, q, 1)));
	}

	return 0;
}

part_global_cached_io::part_global_cached_io(int num_groups,
		io_interface *underlying, int idx, long cache_size,
		int cache_type, file_mapper *mapper): global_cached_io(
			underlying->clone()) {
	this->mapper = mapper->clone();
	this->thread_id = idx;
	this->final_cb = NULL;
	remote_reads = 0;
	this->group_idx = underlying->get_node_id();
	this->cache_size = (long) (cache_size * ((double) underlying->get_local_size()
		/ underlying->get_size()));
	this->cache_type = cache_type;
	this->underlying = underlying;

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
	struct thread_group *group = &groups[group_idx];
	int i = thread_idx(idx, num_groups);
	assert (group->ios[i] == NULL);
	group->ios[i] = this;

	reply_queue = new blocking_FIFO_queue<io_reply>("reply_queue", 
			NUMA_REPLY_QUEUE_SIZE);
	// Create a thread for processing replies.
	reply_processor = new process_reply_thread(this);
	reply_processor->start();
}

/* distribute requests to nodes. */
int part_global_cached_io::distribute_reqs(io_request *requests, int num) {
	int num_sent = 0;
	for (int i = 0; i < num; i++) {
		int idx = hash_req(&requests[i]);
		if (idx != get_group_id()) {
			remote_reads++;
			int ret = req_senders[idx]->send_cached(&requests[i]);
			if (ret == 0)
				break;
		}
		else {
			int ret = global_cached_io::access(&requests[i], 1);
			if (ret < 0) {
				fprintf(stderr, "part global cache can't issue a request\n");
				io_reply rep(&requests[i], false, errno);
				process_reply(&rep);
				processed_replies.inc(1);
			}
		}
		num_sent++;
	}
	for (std::tr1::unordered_map<int, msg_sender<io_request> *>::const_iterator it
			= req_senders.begin(); it != req_senders.end(); it++)
		it->second->flush();
	return num_sent;
}

/* 
 * process the replies and return the number
 * of bytes that have been accessed.
 */
int part_global_cached_io::process_replies(int max_nreplies) {
	int num_processed = 0;
	io_reply local_replies[NUMA_REPLY_BUF_SIZE];
	while(num_processed < max_nreplies) {
		int num = reply_queue->fetch(local_replies, NUMA_REPLY_BUF_SIZE);
		for (int i = 0; i < num; i++) {
			io_reply *reply = &local_replies[i];
			process_reply(reply);
		}
		num_processed += num;
	}
	processed_replies.inc(num_processed);
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
		if (ret > 0) {
			// This variable is only accessed in one thread, so we don't
			// need to protect it.
			processed_requests.inc(ret);
			num_sent += ret;
		}
	}
	return num_sent;
}

void part_global_cached_io::cleanup()
{
	int num_threads = 0;
	for (std::tr1::unordered_map<int, struct thread_group>::const_iterator it
			= groups.begin(); it != groups.end(); it++) {
		num_threads += it->second.nthreads;
	}
	printf("thread %d of %d: start to clean up\n", thread_id, num_threads);

	// First make sure all requests have been flushed for processing.
	for (std::tr1::unordered_map<int, msg_sender<io_request> *>::const_iterator it
			= req_senders.begin(); it != req_senders.end(); it++) {
		msg_sender<io_request> *sender = it->second;
		sender->flush_all();
	}

	// Make sure all threads have finished issuing requests.
	num_finish_issuing_threads.inc(1);
	while (num_finish_issuing_threads.get() < num_threads) {
		usleep(1000 * 10);
	}

	// Now we know no more requests will be put in the request queues.
	// Make sure all request queues empty.
	bool empty;
	do {
		empty = true;
		for (std::tr1::unordered_map<int, struct thread_group>::const_iterator it
				= groups.begin(); it != groups.end(); it++) {
			empty &= it->second.request_queue->is_empty();
		}
		usleep(1000 * 10);
	} while (!empty);

	// Now all requests have been issued to the underlying IOs.
	// Make sure to clean the queues in its own underlying IO.
	underlying->cleanup();

	// Now we need to wait for all requests issued to the disks
	// to be completed.
//	while (processed_requests.get() > processed_replies.get()) {
//		usleep(1000 * 10);
	struct thread_group *group = &groups[group_idx];
	for (size_t i = 0; i < group->process_request_threads.size(); i++) {
		process_request_thread *t
			= (process_request_thread *) group->process_request_threads[i];
		t->get_io()->flush_replies();
	}
//	}
	
	// Let's just exit together.
	num_finished_threads.inc(1);
	while (num_finished_threads.get() < num_threads) {
		usleep(1000 * 10);
	}
	printf("thread %d processed %d requests and %d replies\n", thread_id,
			processed_requests.get(), processed_replies.get());
}

#ifdef STATISTICS
void part_global_cached_io::print_stat()
{
	static long tot_remote_reads = 0;
	static int seen_threads = 0;
	tot_remote_reads += remote_reads;
	seen_threads++;
	if (seen_threads == nthreads) {
		printf("there are %ld requests sent to the remote nodes\n",
				tot_remote_reads);
		int tot_hits = 0;
		for (std::tr1::unordered_map<int, struct thread_group>::const_iterator
				it = groups.begin(); it != groups.end(); it++) {
			const struct thread_group *group = &it->second;
			printf("cache on node %d\n", group->id);
			group->cache->print_stat();
		}
		for (std::tr1::unordered_map<int, struct thread_group>::const_iterator
				it = groups.begin(); it != groups.end(); it++) {
			const struct thread_group *group = &it->second;
			int tot_group_hits = 0;
			for (size_t i = 0; i < group->process_request_threads.size(); i++) {
				process_request_thread *thread = (process_request_thread *) group
					->process_request_threads[i];
				printf("group %d thread %ld gets %ld requests\n", group->id, i,
						thread->get_num_requests());
				tot_group_hits += thread->get_cache_hits();
			}
			printf("group %d gets %d hits\n", group->id, tot_group_hits);
			tot_hits += tot_group_hits;
		}
		printf("There are %d cache hits\n", tot_hits);
	}
}
#endif

std::tr1::unordered_map<int, thread_group> part_global_cached_io::groups;
pthread_mutex_t part_global_cached_io::init_mutex;
int part_global_cached_io::num_finish_init;
pthread_mutex_t part_global_cached_io::wait_mutex;
pthread_cond_t part_global_cached_io::cond;
atomic_integer part_global_cached_io::num_finish_issuing_threads;
atomic_integer part_global_cached_io::num_finished_threads;
