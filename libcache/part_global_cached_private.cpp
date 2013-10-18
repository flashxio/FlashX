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
#include "cache_config.h"
#include "slab_allocator.h"
#include "thread.h"
#include "remote_access.h"

// The size of a request >= sizeof(io_request).
const int NUMA_REQ_BUF_SIZE = NUMA_MSG_SIZE / sizeof(io_request);
// The size of a reply >= sizeof(io_reply).
const int NUMA_REPLY_BUF_SIZE = NUMA_MSG_SIZE / sizeof(io_reply);

const int NUMA_REQ_QUEUE_SIZE = 2000;
const int NUMA_REQ_MSG_QUEUE_SIZE = NUMA_REQ_QUEUE_SIZE / (
		NUMA_MSG_SIZE / sizeof(io_request));
const int NUMA_REPLY_QUEUE_SIZE = 1024;
const int NUMA_REPLY_MSG_QUEUE_SIZE = NUMA_REPLY_QUEUE_SIZE / (
		NUMA_MSG_SIZE / sizeof(io_reply));
const int MAX_REQS_TO_SINGLE_SENDER = IO_MSG_SIZE;

class process_request_thread;
struct thread_group
{
	int id;
	page_cache *cache;
	std::vector<process_request_thread *> process_request_threads;
	slab_allocator *msg_allocator;

	thread_group() {
		id = -1;
		cache = NULL;
		msg_allocator = NULL;
	}
};

class group_request_sender
{
	// It points to the current sender for sending requests.
	unsigned int curr_idx;

	std::vector<process_request_thread *> threads;
	std::vector<request_sender *> senders;
	int num_reqs;

	group_request_sender(slab_allocator *alloc,
			const std::vector<process_request_thread *> &threads);
public:
	static group_request_sender *create(slab_allocator *alloc,
			const std::vector<process_request_thread *> &threads);

	static void destroy(group_request_sender *s) {
		s->~group_request_sender();
		numa_free(s, sizeof(*s));
	}

	int flush();

	int send_cached(io_request *msgs, int num = 1) {
		num_reqs += num;
		if (senders[curr_idx]->get_num_remaining() < MAX_REQS_TO_SINGLE_SENDER)
			return senders[curr_idx]->send_cached(msgs, num);
		else {
			curr_idx++;
			curr_idx = curr_idx % threads.size();
			return senders[curr_idx]->send_cached(msgs, num);
		}
	}

	int get_num_remaining() const {
		int tmp = 0;
		for (unsigned i = 0; i < senders.size(); i++)
			tmp += senders[i]->get_num_remaining();
		assert(tmp == num_reqs);
		return num_reqs;
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
	part_io_process_table(const std::vector<disk_read_thread *> &io_threads,
			file_mapper *mapper, const cache_config *config);

	~part_io_process_table();

	std::tr1::unordered_map<int, group_request_sender *> create_req_senders(
			int node_id) const;
	void destroy_req_senders(const std::tr1::unordered_map<int,
			group_request_sender *> &) const;

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
	// This message buffer is used for copying remote messages to
	// the local memory.
	message<io_request> local_msgs[MSG_BUF_SIZE];
	pthread_t processing_thread_id;
	msg_queue<io_request> *request_queue;

	node_cached_io(io_interface *underlying,
			struct thread_group *local_group);
	~node_cached_io() {
		msg_queue<io_request>::destroy(request_queue);
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

	msg_queue<io_request> *get_queue() const {
		return request_queue;
	}
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
		struct thread_group *local_group): global_cached_io(
			underlying->get_thread(), underlying, local_group->cache)
{
	assert(local_group);
	this->local_group = local_group;
	assert(local_group->id == underlying->get_node_id());

	request_queue = msg_queue<io_request>::create(underlying->get_node_id(),
			"request_queue", NUMA_REQ_MSG_QUEUE_SIZE, INT_MAX, false);

	processed_requests = 0;
	// This IO instance is created inside the right thread.
	global_cached_io::init();

	processing_thread_id = 0;
	pthread_key_create(&replier_key, NULL);

	for (int i = 0; i < MSG_BUF_SIZE; i++) {
		message<io_request> tmp(local_group->msg_allocator, false);
		local_msgs[i] = tmp;
	}
}

/* process the requests sent to this thread */
int node_cached_io::process_requests(int max_nreqs)
{
	if (processing_thread_id == 0)
		processing_thread_id = pthread_self();
	assert(processing_thread_id == pthread_self());
	int num_processed = 0;
	message<io_request> tmp_msgs[MSG_BUF_SIZE];
	// This request buffer has requests more than a message can contain.
	io_request local_reqs[NUMA_REQ_BUF_SIZE];
	while (num_processed < max_nreqs) {
		int num = request_queue->fetch(tmp_msgs, MSG_BUF_SIZE);
		if (num == 0)
			break;

		// We need to copy the message buffers to the local memory.
		for (int i = 0; i < num; i++) {
			tmp_msgs[i].copy_to(local_msgs[i]);
		}
		for (int i = 0; i < num; i++) {
			printf("msg has %d reqs\n", local_msgs[i].get_num_objs());
			while (!local_msgs[i].is_empty()) {
				int num_reqs = local_msgs[i].get_next_objs(local_reqs,
						min(NUMA_REQ_BUF_SIZE,
							global_cached_io::get_remaining_io_slots()));
				global_cached_io::access(local_reqs, num_reqs);
				num_processed += num_reqs;
			}
		}
	}
	// We don't want the thread to be blocked here. The thread is activated
	// in two cases: new requests are sent to its request queue; or some
	// requests have been completed.
	if (global_cached_io::num_pending_ios() > 0)
		global_cached_io::wait4complete(0);

	num_requests += num_processed;
	return num_processed;
}

void part_global_cached_io::notify_upper(io_request *reqs[], int num)
{
	for (int i = 0; i < num; i++) {
		io_request *req = reqs[i];
		io_interface *io = req->get_io();
		assert(io == this);
		if (io->get_callback())
			io->get_callback()->invoke(&req, 1);
	}
}

/**
 * This method is invoked in two places: the underlying IO instance and
 * node_cached_io instance.
 * If it's invoked by the underlying IO, the request was issued by the current
 * thread, so it can just notify the application code.
 * Otherwise, it needs to send a reply message to the thread that issues
 * the request.
 */
void part_global_cached_io::notify_completion(io_request *reqs[], int num)
{
	io_request *local_reqs[num];
	// This buffer is used for sending replies.
#ifdef MEMCHECK
	io_reply *local_reply_buf = new io_reply[num];
#else
	io_reply local_reply_buf[num];
#endif
	int node_id = local_group->id;
	int num_remote = 0;
	int num_local = 0;

	for (int i = 0; i < num; i++) {
		int idx = cache_conf->page2cache(reqs[i]->get_offset());
		if (idx != get_node_id())
			local_reply_buf[num_remote++] = io_reply(reqs[i], true, 0);
		else {
			assert(reqs[i]->get_io() == this);
			local_reqs[num_local++] = reqs[i];
		}
	}

	// The reply must be sent to the thread on a different node.
	int num_sent;
//	if (num_remote > NUMA_REPLY_CACHE_SIZE)
		num_sent = get_reply_sender(node_id)->send(local_reply_buf, num_remote);
//	else
//		num_sent = get_reply_sender(node_id)->send_cached(
//				local_reply_buf, num_remote);
	// We use blocking queues here, so the send must succeed.
	assert(num_sent == num_remote);
#ifdef MEMCHECK
	delete local_reply_buf;
#endif
	get_thread()->activate();
	if (num_local > 0)
		notify_upper(local_reqs, num_local);
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

	process_request_thread(int node_id): thread(
			std::string("process_request_thread-") + itoa(node_id), node_id) {
		this->io = NULL;
	}
public:
	static process_request_thread *create(int node_id) {
		void *addr = numa_alloc_onnode(sizeof(process_request_thread),
				node_id);
		return new(addr) process_request_thread(node_id);
	}

	static void destroy(process_request_thread *t,
			const struct thread_group *group) {
		while (!t->has_exit()) {
			t->stop();
			usleep(10000);
		}
		t->join();
		t->~process_request_thread();
		numa_free(t, sizeof(*t));
	}

	void set_io(node_cached_io *io) {
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

group_request_sender *group_request_sender::create(slab_allocator *alloc,
		const std::vector<process_request_thread *> &threads)
{
	assert(threads.size() > 0);
	int node_id = threads[0]->get_node_id();
	for (unsigned i = 1; i < threads.size(); i++)
		assert(node_id == threads[i]->get_node_id());
	void *addr = numa_alloc_onnode(sizeof(group_request_sender), node_id);
	return new(addr) group_request_sender(alloc, threads);
}

group_request_sender::group_request_sender(slab_allocator *alloc,
		const std::vector<process_request_thread *> &threads)
{
	senders.resize(threads.size());
	for (unsigned i = 0; i < threads.size(); i++) {
		request_sender *sender = request_sender::create(threads[i]->get_node_id(),
				alloc, threads[i]->get_io()->get_queue());
		senders[i] = sender;
	}
	this->threads = threads;
	curr_idx = random() % threads.size();
	num_reqs = 0;
}

int group_request_sender::flush()
{
	int ret = 0;
	for (unsigned i = 0; i < senders.size(); i++) {
		if (senders[i]->get_num_remaining() > 0) {
			ret += senders[i]->flush();
			assert(senders[i]->get_num_remaining() == 0);
			// After we send requests to the remote thread, we need to
			// activate the thread to process them.
			threads[i]->activate();
		}
	}
	num_reqs = 0;
	return ret;
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
		const std::vector<disk_read_thread *> &io_threads,
		file_mapper *mapper, const cache_config *config)
{
	int num_ssds = io_threads.size();
	std::vector<int> node_ids;
	config->get_node_ids(node_ids);
	cache_conf = config;
	num_groups = node_ids.size();
	for (unsigned i = 0; i < node_ids.size(); i++) {
		int node_id = node_ids[i];
		struct thread_group group;

		group.id = node_id;
		group.cache = cache_conf->create_cache_on_node(node_id,
				MAX_NUM_FLUSHES_PER_FILE * num_ssds / config->get_num_caches());
		group.msg_allocator = new slab_allocator(NUMA_MSG_SIZE,
				NUMA_MSG_SIZE * 128, INT_MAX, node_id);
		groups.insert(std::pair<int, struct thread_group>(node_id, group));

		struct thread_group *groupp = &groups[node_id];
		// Create processing threads.
		for (int i = 0; i < params.get_numa_num_process_threads(); i++) {
			process_request_thread *t = process_request_thread::create(node_id);
			node_cached_io *io = node_cached_io::create(
					new remote_disk_access(io_threads, mapper, t), groupp);
			t->set_io(io);
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
		for (unsigned i = 0; i < group->process_request_threads.size(); i++) {
			process_request_thread *t = group->process_request_threads[i];
			node_cached_io *io = t->get_io();
			process_request_thread::destroy(t, group);
			node_cached_io::destroy(io);
		}
		cache_conf->destroy_cache_on_node(group->cache);
	}
}

std::tr1::unordered_map<int, group_request_sender *>
part_io_process_table::create_req_senders(int node_id) const
{
	std::tr1::unordered_map<int, group_request_sender *> req_senders;
	// Initialize the request senders.
	for (std::map<int, struct thread_group>::const_iterator it
			= groups.begin(); it != groups.end(); it++) {
		const struct thread_group *group = &it->second;
		assert(group->id == it->first);
		req_senders.insert(std::pair<int, group_request_sender *>(group->id,
					group_request_sender::create(group->msg_allocator,
						group->process_request_threads)));
	}
	return req_senders;
}

void part_io_process_table::destroy_req_senders(
		const std::tr1::unordered_map<int, group_request_sender *> &req_senders) const
{
	for (std::tr1::unordered_map<int, group_request_sender *>::const_iterator it
			= req_senders.begin(); it != req_senders.end(); it++) {
		group_request_sender::destroy(it->second);
	}
}

int part_global_cached_io::init()
{
	underlying->init();

	return 0;
}

int part_global_cached_io::close_file(part_io_process_table *table)
{
	delete table;
	return 0;
}

part_io_process_table *part_global_cached_io::open_file(
		const std::vector<disk_read_thread *> &io_threads,
		file_mapper *mapper, const cache_config *config)
{
	return new part_io_process_table(io_threads, mapper, config);
}

part_global_cached_io::part_global_cached_io(io_interface *underlying,
		part_io_process_table *table): io_interface(underlying->get_thread())
{
	int node_id = underlying->get_node_id();
	this->underlying = new global_cached_io(underlying->get_thread(), underlying,
			table->get_cache(node_id));
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
	reply_queue = msg_queue<io_reply>::create(node_id, "reply_queue", 
			// We don't want that adding replies is blocked, so we allow
			// the reply queue to be expanded to an arbitrary size.
			NUMA_REPLY_MSG_QUEUE_SIZE, INT_MAX / sizeof(io_reply), false);

	std::set<int> node_ids = table->get_node_ids();
	for (std::set<int>::const_iterator it = node_ids.begin();
			it != node_ids.end(); it++) {
		int node_id = *it;
		if ((int) reply_senders.size() <= node_id)
			reply_senders.resize(node_id + 1);
		reply_senders[node_id] = thread_safe_msg_sender<io_reply>::create(
				node_id, local_group->msg_allocator, reply_queue);
	}

	req_senders = table->create_req_senders(node_id);

	for (int i = 0; i < MSG_BUF_SIZE; i++) {
		message<io_reply> tmp(local_group->msg_allocator, false);
		local_reply_msgs[i] = tmp;
	}
}

part_global_cached_io::~part_global_cached_io()
{
	msg_queue<io_reply>::destroy(reply_queue);
	for (unsigned i = 0; i < reply_senders.size(); i++) {
		if (reply_senders[i])
			thread_safe_msg_sender<io_reply>::destroy(reply_senders[i]);
	}
	global_table->destroy_req_senders(req_senders);
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
	message<io_reply> tmp_msgs[MSG_BUF_SIZE];
	// The reply buffer size should be able to contain more replies than
	// a reply message.
	io_reply local_reply_buf[NUMA_REPLY_BUF_SIZE];
	while (!reply_queue->is_empty()) {
		int num = reply_queue->fetch(tmp_msgs, MSG_BUF_SIZE);
		assert(num <= MSG_BUF_SIZE);
		for (int i = 0; i < num; i++) {
			tmp_msgs[i].copy_to(local_reply_msgs[i]);
		}
		for (int i = 0; i < num; i++) {
			int num_replies = local_reply_msgs[i].get_next_objs(
					local_reply_buf, NUMA_REPLY_BUF_SIZE);
			assert(local_reply_msgs[i].is_empty());
			for (int j = 0; j < num_replies; j++) {
				io_reply *reply = &local_reply_buf[j];
				this->process_reply(reply);
			}
			num_processed += num_replies;
		}
	}
	processed_replies += num_processed;
	return num_processed;
}

void part_global_cached_io::access(io_request *requests, int num, io_status status[])
{
	ASSERT_EQ(get_thread(), thread::get_curr_thread());
	// TODO I'll write status to the status array later.
	int num_sent = 0;
	int num_local_reqs = 0;
	// Distribute requests to the corresponding nodes.
	for (int i = 0; i < num; i++) {
		assert(requests[i].within_1page());
		int idx = cache_conf->page2cache(requests[i].get_offset());
		if (idx != get_node_id()) {
			remote_reads++;
			// TODO why do this?
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

	underlying->access(local_req_buf, num_local_reqs);
	// TODO how to deal with error of access.
	
	// This variable is only accessed in one thread, so we don't
	// need to protect it.
	processed_requests += num;

	// Let's process all replies before proceeding.
	process_replies();
}

void part_global_cached_io::cleanup()
{
	// First make sure all requests have been flushed for processing.
	for (std::tr1::unordered_map<int, group_request_sender *>::const_iterator it
			= req_senders.begin(); it != req_senders.end(); it++) {
		group_request_sender *sender = it->second;
		// Send a flush request to force request process threads to flush all requests.
		io_request flush_req(true);
		sender->send_cached(&flush_req);
		sender->flush();
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
		underlying->cleanup();

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
	tot_hits += underlying->get_cache_hits();
	tot_fast_process += underlying->get_num_fast_process();
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

int part_global_cached_io::wait4complete(int num_to_complete)
{
	// We first send all requests in the local buffer to their destinations.
	flush_requests();

	int pending = num_pending_ios();
	num_to_complete = min(pending, num_to_complete);
	// We process the completed requests sent to remote nodes.
	process_replies();
	// We process the completed requests in the local IO instance.
	// But we don't want to be blocked.
	underlying->wait4complete(0);
	while (pending - num_pending_ios() < num_to_complete) {
		get_thread()->wait();
		process_replies();
		underlying->wait4complete(0);
	}
	return pending - num_pending_ios();
}

void part_global_cached_io::flush_requests()
{
	for (std::tr1::unordered_map<int, group_request_sender *>::const_iterator it
			= req_senders.begin(); it != req_senders.end(); it++) {
		it->second->flush();
		int ret = it->second->get_num_remaining();
		assert(ret == 0);
	}
	underlying->flush_requests();
}
