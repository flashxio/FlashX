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

#include <set>

#include "cache.h"
#include "associative_cache.h"
#include "global_cached_private.h"
#include "part_global_cached_private.h"
#include "parameters.h"

class process_request_thread;
struct thread_group
{
	int id;
	page_cache *cache;
	std::vector<process_request_thread *> process_request_threads;
	const io_interface *underlying;
	blocking_FIFO_queue<io_request> *request_queue;

	thread_group() {
		id = -1;
		cache = NULL;
		request_queue = NULL;
		underlying = NULL;
	}
};

/**
 * The IO processing table is created for each virtual file.
 * Therefore, each virtual file has a group of dedicated threads created
 * for processing I/O requests. This is a simpler solution compared with
 * the alternative that all virtual files share the same set of I/O
 * processing threads. The alternative requires many table lookups and
 * sorting IO requests for each file, etc. This solution relies on
 * the thread scheduling in the OS. It should works well when there are
 * only a few files open. These are threads, so the context switch should
 * be rather cheap.
 * Another benefit of this solution is that no locks are required.
 */
class part_io_process_table
{
	std::map<int, struct thread_group> groups;
	const cache_config *cache_conf;
	int num_groups;
public:
	part_io_process_table(std::map<int, io_interface *> &underlyings,
			const cache_config *config);

	~part_io_process_table();

	std::tr1::unordered_map<int, request_sender *> create_req_senders(
			int node_id) const;
	void destroy_req_senders(const std::tr1::unordered_map<int,
			request_sender *> &) const;

	const cache_config *get_cache_config() const {
		return cache_conf;
	}

	const struct thread_group *get_thread_group(int node_id) const {
		std::map<int, struct thread_group>::const_iterator it
			= groups.find(node_id);
		return &it->second;
	}

	page_cache *get_cache(int node_id) {
		std::map<int, struct thread_group>::const_iterator it
			= groups.find(node_id);
		return it->second.cache;
	}

	const io_interface *get_underlying_io(int node_id) const {
		std::map<int, struct thread_group>::const_iterator it
			= groups.find(node_id);
		return it->second.underlying;
	}

	int get_num_groups() const {
		return num_groups;
	}

	std::set<int> get_node_ids() const {
		std::set<int> node_ids;
		for (std::map<int, struct thread_group>::const_iterator it
				= groups.begin(); it != groups.end(); it++)
			node_ids.insert(it->first);
		return node_ids;
	}
};

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
	const struct thread_group *local_group;

	// thread id <-> msg sender
	pthread_key_t replier_key;
	long processed_requests;
	long num_requests;
	io_request local_msg_reqs[NUMA_REQ_BUF_SIZE];
	io_request local_reqs[NUMA_REQ_BUF_SIZE];
	io_reply local_reply_buf[REPLY_BUF_SIZE];
	pthread_t processing_thread_id;

	node_cached_io(io_interface *underlying, struct thread_group *local_group);
	~node_cached_io() {
	}
public:
	static node_cached_io *create(io_interface *underlying,
			struct thread_group *local_group) {
		assert(underlying->get_node_id() >= 0);
		void *addr = numa_alloc_onnode(sizeof(node_cached_io),
				underlying->get_node_id());
		return new(addr) node_cached_io(underlying, local_group);
	}

	static void destroy(node_cached_io *io) {
		io->~node_cached_io();
		numa_free(io, sizeof(*io));
	}

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
	global_cached_io(underlying, local_group->cache)
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
	blocking_FIFO_queue<io_request> *request_queue = local_group->request_queue;
	while (num_processed < max_nreqs) {
		int num = request_queue->fetch(local_msg_reqs, NUMA_REQ_BUF_SIZE,
				true, true);
		// We have been interrupted from waiting for IO requests.
		// Maybe it's a signal for stopping the thread.
		if (num == 0)
			break;

		for (int i = 0; i < num; i++)
			local_reqs[i] = local_msg_reqs[i];
		global_cached_io::access(local_reqs, num);
		flush_requests();
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

	process_request_thread(node_cached_io *io): thread(
			std::string("process_request_thread-") + itoa(io->get_node_id()),
			// We don't use blocking mode of the thread because
			// we are using blocking queue.
			io->get_node_id(), false) {
		this->io = io;
	}
public:
	static process_request_thread *create(node_cached_io *io) {
		assert(io->get_node_id() >= 0);
		void *addr = numa_alloc_onnode(sizeof(process_request_thread),
				io->get_node_id());
		return new(addr) process_request_thread(io);
	}

	static void destroy(process_request_thread *t,
			const struct thread_group *group) {
		while (!t->has_exit()) {
			t->stop();
			group->request_queue->wakeup();
			usleep(10000);
		}
		t->join();
		t->~process_request_thread();
		numa_free(t, sizeof(*t));
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

part_io_process_table::part_io_process_table(
		std::map<int, io_interface *> &underlyings,
		const cache_config *config)
{
	cache_conf = config;
	num_groups = underlyings.size();
	for (std::map<int, io_interface *>::const_iterator it = underlyings.begin();
			it != underlyings.end(); it++) {
		int node_id = it->first;
		const io_interface *underlying = it->second;
		assert(node_id == underlying->get_node_id());

		struct thread_group group;
		group.id = node_id;
		group.underlying = underlying;
		group.cache = cache_conf->create_cache_on_node(node_id);
		group.request_queue = blocking_FIFO_queue<io_request>::create(node_id, "request_queue", 
				NUMA_REQ_QUEUE_SIZE, NUMA_REQ_QUEUE_SIZE);
		assert(underlying);
		groups.insert(std::pair<int, struct thread_group>(node_id, group));

		struct thread_group *groupp = &groups[node_id];
		// Create processing threads.
		for (int i = 0; i < NUMA_NUM_PROCESS_THREADS; i++) {
			node_cached_io *io = node_cached_io::create(underlying->clone(), groupp);
			process_request_thread *t = process_request_thread::create(io);
			t->start();
			groupp->process_request_threads.push_back(t);
		}
	}
}

part_io_process_table::~part_io_process_table()
{
	for (std::map<int, struct thread_group>::const_iterator it
			= groups.begin(); it != groups.end(); it++) {
		const struct thread_group *group = &it->second;
		for (int i = 0; i < NUMA_NUM_PROCESS_THREADS; i++) {
			process_request_thread *t = group->process_request_threads[i];
			node_cached_io *io = t->get_io();
			process_request_thread::destroy(t, group);
			node_cached_io::destroy(io);
		}
		cache_conf->destroy_cache_on_node(group->cache);
		blocking_FIFO_queue<io_request>::destroy(group->request_queue);
	}
}

std::tr1::unordered_map<int, request_sender *>
part_io_process_table::create_req_senders(int node_id) const
{
	std::tr1::unordered_map<int, request_sender *> req_senders;
	// Initialize the request senders.
	for (std::map<int, struct thread_group>::const_iterator it
			= groups.begin(); it != groups.end(); it++) {
		const struct thread_group *group = &it->second;
		assert(group->id == it->first);
		req_senders.insert(std::pair<int, request_sender *>(group->id,
					request_sender::create(node_id, group->request_queue,
						NUMA_REQ_CACHE_SIZE)));
	}
	return req_senders;
}

void part_io_process_table::destroy_req_senders(
		const std::tr1::unordered_map<int, request_sender *> &req_senders) const
{
	for (std::tr1::unordered_map<int, request_sender *>::const_iterator it
			= req_senders.begin(); it != req_senders.end(); it++) {
		request_sender::destroy(it->second);
	}
}

int part_global_cached_io::init()
{
	global_cached_io::init();

	return 0;
}

int part_global_cached_io::close_file(part_io_process_table *table)
{
	delete table;
	return 0;
}

part_io_process_table *part_global_cached_io::open_file(
		std::map<int, io_interface *> &underlyings,
		const cache_config *config)
{
	return new part_io_process_table(underlyings, config);
}

part_global_cached_io::part_global_cached_io(int node_id,
		part_io_process_table *table): global_cached_io(
			table->get_underlying_io(node_id)->clone(),
			table->get_cache(node_id))
{
	processed_requests = 0;;
	sent_requests = 0;
	processed_replies = 0;

	this->final_cb = NULL;
	remote_reads = 0;
	this->cache_conf = table->get_cache_config();
	this->global_table = table;

#ifdef DEBUG
	long cache_size = cache_conf->get_part_size(node_id);
	printf("thread id: %d, group id: %d, num groups: %d, cache size: %ld\n",
			get_io_idx(), node_id, table->get_num_groups(), cache_size);
#endif

	/* assign a thread to a group. */
	this->local_group = table->get_thread_group(node_id);
	reply_queue = blocking_FIFO_queue<io_reply>::create(node_id, "reply_queue", 
			// We don't want that adding replies is blocked, so we allow
			// the reply queue to be expanded to an arbitrary size.
			NUMA_REPLY_QUEUE_SIZE, INT_MAX / sizeof(io_reply));

	std::set<int> node_ids = table->get_node_ids();
	for (std::set<int>::const_iterator it = node_ids.begin();
			it != node_ids.end(); it++) {
		int node_id = *it;
		if ((int) reply_senders.size() <= node_id)
			reply_senders.resize(node_id + 1);
		reply_senders[node_id] = thread_safe_msg_sender<io_reply>::create(node_id,
				NUMA_REPLY_CACHE_SIZE, reply_queue);
	}

	req_senders = table->create_req_senders(node_id);
}

part_global_cached_io::~part_global_cached_io()
{
	blocking_FIFO_queue<io_reply>::destroy(reply_queue);
	for (unsigned i = 0; i < reply_senders.size(); i++) {
		if (reply_senders[i])
			thread_safe_msg_sender<io_reply>::destroy(reply_senders[i]);
	}
	global_table->destroy_req_senders(req_senders);
}

/* distribute requests to nodes. */
int part_global_cached_io::distribute_reqs(io_request *requests, int num) {
	int num_sent = 0;
	int num_local_reqs = 0;
	for (int i = 0; i < num; i++) {
		assert(requests[i].within_1page());
		int idx = cache_conf->page2cache(requests[i].get_offset());
		if (idx != get_node_id()) {
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
	global_cached_io::cleanup();

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
		global_cached_io::cleanup();

		usleep(1000 * 10);
	}

#ifdef STATISTICS
	printf("thread %d processed %ld requests (%ld remote requests) and %ld replies\n",
			get_io_idx(), processed_requests, sent_requests, processed_replies);
#endif
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
		std::set<int> node_ids = global_table->get_node_ids();
		for (std::set<int>::const_iterator it = node_ids.begin();
				it != node_ids.end(); it++) {
			int node_id = *it;
			const struct thread_group *group = global_table->get_thread_group(node_id);
			printf("cache on node %d\n", group->id);
			group->cache->print_stat();
		}
		for (std::set<int>::const_iterator it = node_ids.begin();
				it != node_ids.end(); it++) {
			int node_id = *it;
			const struct thread_group *group = global_table->get_thread_group(node_id);
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
		if (cache_conf->page2cache(offset) != get_node_id())
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
