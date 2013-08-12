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
#include <limits.h>

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
 * This wraps global cached IO to handle requests sent from remote threads.
 */
class node_cached_io: public global_cached_io
{
	// thread id <-> msg sender
	pthread_key_t replier_key;
	const struct thread_group *local_group;
	long processed_requests;
	long num_requests;
	io_request local_msg_reqs[NUMA_REQ_BUF_SIZE];
	io_request local_reqs[NUMA_REQ_BUF_SIZE];
	io_reply local_reply_buf[REPLY_BUF_SIZE];

	pthread_t processing_thread_id;
public:
	node_cached_io(io_interface *underlying, struct thread_group *local_group);

	virtual page_cache *get_global_cache() {
		return local_group->cache;
	}

	long get_num_requests() const {
		return num_requests;
	}

	int process_requests(int max_nreqs);

	// Requests are from remote threads, so we need to send replies
	// to them for notification.
	virtual void notify_completion(io_request *req);
	virtual void notify_completion(io_request *reqs[], int num);
};

#if 0
/**
 * This thread runs on the node with the application threads.
 * It is dedicated to processing replies. If the queue is empty,
 * the thread will sleep.
 */
class process_reply_thread: public thread
{
	part_global_cached_io *io;
	blocking_FIFO_queue<io_reply> *reply_queue;
public:
	process_reply_thread(part_global_cached_io *io,
			blocking_FIFO_queue<io_reply> *reply_queue): thread(
			std::string("process_reply_thread-") + itoa(io->get_node_id()),
			io->get_node_id(), false) {
		this->io = io;
		this->reply_queue = reply_queue;
	}
	void run();

	blocking_FIFO_queue<io_reply> *get_reply_queue() {
		return reply_queue;
	}
};
#endif

node_cached_io::node_cached_io(io_interface *underlying,
		struct thread_group *local_group):
	global_cached_io(underlying->clone())
{
	assert(local_group);
	this->local_group = local_group;

	processed_requests = 0;
	// This IO instance is created inside the right thread.
	global_cached_io::init();

	processing_thread_id = 0;
	pthread_key_create(&replier_key, NULL);
}

/* process the requests sent to this thread */
int node_cached_io::process_requests(int max_nreqs)
{
	if (processing_thread_id == 0)
		processing_thread_id = pthread_self();
	assert(processing_thread_id == pthread_self());
	int num_processed = 0;
	fifo_queue<io_request> *request_queue = local_group->request_queue;
	while (num_processed < max_nreqs) {
		int num = request_queue->fetch(local_msg_reqs, NUMA_REQ_BUF_SIZE);
		for (int i = 0; i < num; i++)
			local_reqs[i] = local_msg_reqs[i];
		global_cached_io::access(local_reqs, num);
		num_processed += num;
	}
	num_requests += num_processed;
	return num_processed;
}

/*
 * We sort requests according to the io instance where they will be delivered to.
 */
void sort_requests(io_request *reqs[], int num,
		std::tr1::unordered_map<io_interface *, std::vector<io_request *> > &req_map)
{
	for (int i = 0; i < num; i++) {
		io_request *req = reqs[i];
		io_interface *io = req->get_io();
		std::vector<io_request *> *vec;
		std::tr1::unordered_map<io_interface *,
			std::vector<io_request *> >::iterator it = req_map.find(io);
		if (it == req_map.end()) {
			req_map.insert(std::pair<io_interface *, std::vector<io_request *> >(
						io, std::vector<io_request *>()));
			vec = &req_map[io];
		}
		else
			vec = &it->second;
		vec->push_back(req);
	}
}

void node_cached_io::notify_completion(io_request *reqs[], int num)
{
	std::tr1::unordered_map<io_interface *, std::vector<io_request *> > req_map;
	sort_requests(reqs, num, req_map);
	int node_id = local_group->id;
	for (std::tr1::unordered_map<io_interface *,
			std::vector<io_request *> >::iterator it = req_map.begin();
			it != req_map.end(); it++) {
		part_global_cached_io *io = (part_global_cached_io *) it->first;
		std::vector<io_request *> *vec = &it->second;

		for (size_t i = 0; i < vec->size(); i++)
			local_reply_buf[i] = io_reply(vec->at(i), true, 0);
		// The reply must be sent to the thread on a different node.
		assert(io->get_node_id() != this->get_node_id());
		int num_sent;
		if ((int) vec->size() > NUMA_REPLY_CACHE_SIZE)
			num_sent = io->get_reply_sender(node_id)->send(local_reply_buf,
					vec->size());
		else
			num_sent = io->get_reply_sender(node_id)->send_cached(
					local_reply_buf, vec->size());
		// We use blocking queues here, so the send must succeed.
		assert(num_sent == (int) vec->size());
	}
}

void node_cached_io::notify_completion(io_request *req)
{
	io_reply rep(req, true, 0);
	part_global_cached_io *io = (part_global_cached_io *) req->get_io();
	// The reply must be sent to the thread on a different node.
	assert(io->get_node_id() != this->get_node_id());
	int node_id = local_group->id;
	int num_sent = io->get_reply_sender(node_id)->send_cached(&rep);
	// We use blocking queues here, so the send must succeed.
	assert(num_sent > 0);
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

	int get_num_fast_process() const {
		return io->get_num_fast_process();
	}

	long get_num_requests() const {
		return io->get_num_requests();
	}
};

void process_request_thread::run()
{
	io->process_requests(NUMA_REQ_BUF_SIZE);
}

#if 0
void process_reply_thread::run()
{
	int max_nreplies = NUMA_REPLY_BUF_SIZE;
	int num_processed = 0;
	io_reply local_replies[NUMA_REPLY_BUF_SIZE];
	while(num_processed < max_nreplies) {
		int num = reply_queue->fetch(local_replies, NUMA_REPLY_BUF_SIZE);
		for (int i = 0; i < num; i++) {
			io_reply *reply = &local_replies[i];
			((part_global_cached_io *) reply->get_io())->process_reply(reply);
		}
		num_processed += num;
	}
}
#endif

int part_global_cached_io::init() {
	global_cached_io::init();

	/*
	 * we have to initialize senders here
	 * because we need to make sure other threads have 
	 * initialized all queues.
	 */
	for (std::tr1::unordered_map<int, struct thread_group>::const_iterator it
			= groups.begin(); it != groups.end(); it++) {
		const struct thread_group *group = &it->second;
		req_senders.insert(std::pair<int, request_sender *>(it->first,
					new request_sender(group->request_queue, NUMA_REQ_CACHE_SIZE)));
	}

	return 0;
}

// This function initializes all processing threads.
int part_global_cached_io::init_io_system(const std::vector<int> &node_id_array,
		std::map<int, io_interface *> &underlyings,
		const cache_config *cache_conf)
{
	pthread_mutex_lock(&init_mutex);
	if (groups.size() == 0) {
		for (unsigned i = 0 ; i < node_id_array.size(); i++) {
			int node_id = node_id_array[i];
			groups.insert(std::pair<int, struct thread_group>(node_id, thread_group()));
			struct thread_group *group = &groups[node_id];

			group->id = node_id;
			group->cache = cache_conf->create_cache_on_node(node_id);
			group->request_queue = new blocking_FIFO_queue<io_request>("request_queue", 
					NUMA_REQ_QUEUE_SIZE, NUMA_REQ_QUEUE_SIZE);
			io_interface *underlying = underlyings[node_id];
			assert(underlying);
			// Create processing threads.
			for (int i = 0; i < NUMA_NUM_PROCESS_THREADS; i++) {
				node_cached_io *io = new node_cached_io(underlying->clone(), group);
				thread *t = new process_request_thread(io);
				t->start();
				group->process_request_threads.push_back(t);
			}
		}
	}
	pthread_mutex_unlock(&init_mutex);
	return 0;
}

part_global_cached_io::part_global_cached_io(int num_groups,
		io_interface *underlying, const cache_config *config): global_cached_io(
			underlying->clone())
{
	processed_requests = 0;;
	sent_requests = 0;
	processed_replies = 0;

	this->final_cb = NULL;
	remote_reads = 0;
	this->group_idx = underlying->get_node_id();
	this->underlying = underlying;
	this->cache_conf = config;

	long cache_size = cache_conf->get_part_size(underlying->get_node_id());
	printf("thread id: %d, group id: %d, num groups: %d, cache size: %ld\n",
			get_io_idx(), group_idx, num_groups, cache_size);

	/* assign a thread to a group. */
	struct thread_group *group = &groups[group_idx];
	this->local_group = group;
	reply_queue = new blocking_FIFO_queue<io_reply>("reply_queue", 
			// We don't want that adding replies is blocked, so we allow
			// the reply queue to be expanded to an arbitrary size.
			NUMA_REPLY_QUEUE_SIZE, INT_MAX / sizeof(io_reply));

	fifo_queue<io_reply> *q1[1];
	q1[0] = reply_queue;
	for (unsigned i = 0; i < groups.size(); i++) {
		struct thread_group *group = &groups[i];
		if ((int) reply_senders.size() <= group->id)
			reply_senders.resize(group->id + 1);
		reply_senders[group->id] = new thread_safe_msg_sender<io_reply>(
				NUMA_REPLY_CACHE_SIZE, q1, 1);
	}
}

/* distribute requests to nodes. */
int part_global_cached_io::distribute_reqs(io_request *requests, int num) {
	int num_sent = 0;
	int num_local_reqs = 0;
	for (int i = 0; i < num; i++) {
		assert(requests[i].within_1page());
		int idx = cache_conf->page2cache(requests[i].get_offset());
		if (idx != get_group_id()) {
			remote_reads++;
			io_request req(requests[i]);
			int ret = req_senders[idx]->send_cached(&req);
			sent_requests++;
			assert(ret > 0);
		}
		else {
			assert(num_local_reqs < REQ_BUF_SIZE);
			local_req_buf[num_local_reqs++] = requests[i];
		}
		num_sent++;
	}

	global_cached_io::access(local_req_buf, num_local_reqs);
	// TODO how to deal with error of access.
	int num_remaining = 0;
	for (std::tr1::unordered_map<int, request_sender *>::const_iterator it
			= req_senders.begin(); it != req_senders.end(); it++) {
		it->second->flush(false);
		num_remaining += it->second->get_num_remaining();
	}
	return num_remaining;
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
			reply->get_access_method(), this, -1);
	io_request *reqs[1] = {&req};
	final_cb->invoke(reqs, 1);
	return ret;
}

int part_global_cached_io::process_replies()
{
	int num_processed = 0;
	while (!reply_queue->is_empty()) {
		int num = reply_queue->non_blocking_fetch(local_reply_buf, REPLY_BUF_SIZE);
		for (int i = 0; i < num; i++) {
			io_reply *reply = &local_reply_buf[i];
			this->process_reply(reply);
		}
		num_processed += num;
	}
	processed_replies += num_processed;
	return num_processed;
}

void part_global_cached_io::access(io_request *requests, int num, io_status status[])
{
	// TODO I'll write status to the status array later.
	int num_sent = 0;
	int num_remaining = distribute_reqs(&requests[num_sent], num - num_sent);
	// This variable is only accessed in one thread, so we don't
	// need to protect it.
	processed_requests += num;
	num_sent += num;

	// Let's process all replies before proceeding.
	process_replies();
	while (num_remaining > 0) {
		num_remaining = 0;
		for (std::tr1::unordered_map<int, request_sender *>::const_iterator it
				= req_senders.begin(); it != req_senders.end(); it++) {
			it->second->flush(true);
			num_remaining += it->second->get_num_remaining();
		}
	}
}

void part_global_cached_io::cleanup()
{
	// First make sure all requests have been flushed for processing.
	for (std::tr1::unordered_map<int, request_sender *>::const_iterator it
			= req_senders.begin(); it != req_senders.end(); it++) {
		request_sender *sender = it->second;
		// Send a flush request to force request process threads to flush all requests.
		io_request flush_req;
		sender->send_cached(&flush_req);
		sender->flush_all();
	}

	// Now all requests have been issued to the underlying IOs.
	// Make sure to clean the queues in its own underlying IO.
	underlying->cleanup();

	// Now we need to wait for all requests issued to the disks
	// to be completed.
	while (sent_requests > processed_replies) {
		// flush all replies in the reply senders, so this IO can process
		// the remaining replies.
		for (unsigned i = 0; i < reply_senders.size(); i++) {
			thread_safe_msg_sender<io_reply> *sender = reply_senders[i];
			if (sender)
				sender->flush_all();
		}
		process_replies();

		usleep(1000 * 10);
	}

	printf("thread %d processed %ld requests (%ld remote requests) and %ld replies\n",
			get_io_idx(), processed_requests, sent_requests, processed_replies);
}

#ifdef STATISTICS
void part_global_cached_io::print_stat(int nthreads)
{
	static long tot_remote_reads = 0;
	static int seen_threads = 0;
	static int tot_hits = 0;
	static int tot_fast_process = 0;
	tot_remote_reads += remote_reads;
	seen_threads++;
	tot_hits += get_cache_hits();
	tot_fast_process += this->get_num_fast_process();
	if (seen_threads == nthreads) {
		printf("there are %ld requests sent to the remote nodes\n",
				tot_remote_reads);
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
				tot_fast_process += thread->get_num_fast_process();
			}
			printf("group %d gets %d hits\n", group->id, tot_group_hits);
			tot_hits += tot_group_hits;
		}
		printf("There are %d cache hits\n", tot_hits);
		printf("There are %d requests processed in the fast path\n", tot_fast_process);
	}
}
#endif

int part_global_cached_io::preload(off_t start, long size)
{
	if (size > cache_conf->get_size()) {
		fprintf(stderr, "we can't preload data larger than the cache size\n");
		exit(1);
	}

	assert(ROUND_PAGE(start) == start);
	for (long offset = start; offset < start + size; offset += PAGE_SIZE) {
		off_t old_off = -1;
		// We only preload data to the local cache.
		if (cache_conf->page2cache(offset) != get_group_id())
			continue;
		thread_safe_page *p = (thread_safe_page *) local_group->cache->search(
					ROUND_PAGE(offset), old_off);
		// This is mainly for testing. I don't need to really read data from disks.
		if (!p->data_ready()) {
			p->set_io_pending(false);
			p->set_data_ready(true);
		}
		p->dec_ref();
	}
	return 0;
}

std::tr1::unordered_map<int, thread_group> part_global_cached_io::groups;
pthread_mutex_t part_global_cached_io::init_mutex = PTHREAD_MUTEX_INITIALIZER;
