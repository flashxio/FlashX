/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

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
#include "debugger.h"
#include "timer.h"

const int NUMA_REQ_QUEUE_SIZE = 2000;
const int NUMA_REQ_MSG_QUEUE_SIZE = NUMA_REQ_QUEUE_SIZE / (
		NUMA_MSG_SIZE / sizeof(io_request));
const int NUMA_REPLY_QUEUE_SIZE = 1024;
const int NUMA_REPLY_MSG_QUEUE_SIZE = NUMA_REPLY_QUEUE_SIZE / (
		NUMA_MSG_SIZE / sizeof(io_reply));

class process_request_thread;
class underlying_io_thread;
struct thread_group
{
	int id;
	page_cache *cache;
	std::vector<process_request_thread *> process_request_threads;
	// This processes the requests that are issued by global_cached_io and
	// are supposed to be sent to the I/O threads.
	underlying_io_thread *underlying_thread;
	slab_allocator *msg_allocator;

	thread_group() {
		id = -1;
		cache = NULL;
		msg_allocator = NULL;
		underlying_thread = NULL;
	}

	void print_state();
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
		return senders[curr_idx]->send_cached(msgs, num);
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
	part_io_process_table(const std::vector<disk_io_thread *> &io_threads,
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

	int get_num_threads() const {
		int nthreads = 0;
		for (std::map<int, struct thread_group>::const_iterator it
				= groups.begin(); it != groups.end(); it++)
			nthreads += it->second.process_request_threads.size();
		return nthreads;
	}

	std::set<int> get_node_ids() const {
		std::set<int> node_ids;
		for (std::map<int, struct thread_group>::const_iterator it
				= groups.begin(); it != groups.end(); it++)
			node_ids.insert(it->first);
		return node_ids;
	}

	/**
	 * Print the statistics info of request processing threads.
	 * `nthreads': the number of user threads.
	 */
	void print_stat(int num_user_threads);

	void print_state() {
		for (std::map<int, struct thread_group>::iterator it
				= groups.begin(); it != groups.end(); it++) {
			it->second.print_state();
		}
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
	// The size of a request >= sizeof(io_request).
	static const int NUMA_REQ_BUF_SIZE = 1024;

	const struct thread_group *local_group;

	// thread id <-> msg sender
	pthread_key_t replier_key;
	long processed_requests;
	long num_requests;

	// This message buffer is used for copying remote messages to
	// the local memory.
	message<io_request> local_msgs[MSG_BUF_SIZE];
	// This request buffer has requests more than a message can contain.
	io_request local_reqs[NUMA_REQ_BUF_SIZE];
	// A temporary message buffer
	message<io_request> tmp_msg_buf[MSG_BUF_SIZE];

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

	int process_requests();

	msg_queue<io_request> *get_queue() const {
		return request_queue;
	}

	void print_state() {
		printf("node cached io %d has %d pending requests and %d reqs in the queue\n",
				get_io_idx(), num_pending_ios(), request_queue->get_num_objs());
		global_cached_io::print_state();
	}

	// This method is used for IO requests that need to be sent to disks.
	// Thanks to the IO proxy, it runs in the underlying IO thread,
	// which is on the same NUMA node as the request process thread where
	// the IO instance belongs to.
	virtual void notify_completion(io_request *reqs[], int num) {
		assert(thread::get_curr_thread()->get_node_id() == get_thread()->get_node_id());
		stack_array<io_request> req_copies(num);
		for (int i = 0; i < num; i++) {
			req_copies[i] = *reqs[i];
			assert(req_copies[i].get_io());
		}
		process_disk_completed_requests(req_copies.data(), num);
		process_completed_requests();
		if (has_pending_requests())
			get_thread()->activate();
	}

	int get_file_id() const {
		assert(0);
		return -1;
	}
};

static void notify_gcached_io(io_interface *io, io_request *reqs[], int num)
{
	// TODO we may to more than this since we run in the same NUMA node
	// as the IO instance.
	io->notify_completion(reqs, num);
}

/**
 * This IO proxy is responsible for sending requests to disks and gather
 * completion notifications.
 */
class underlying_io_proxy: public remote_io
{
public:
	underlying_io_proxy(const std::vector<disk_io_thread *> &remotes,
			file_mapper *mapper, thread *curr_thread): remote_io(
				remotes, mapper, curr_thread) {
	}

	// This runs in underlying_io_thread.
	virtual int process_completed_requests(io_request reqs[], int num) {
		stack_array<io_request *, 256> req_ptrs(num);
		for (int i = 0; i < num; i++) {
			reqs[i].set_io((io_interface *) reqs[i].get_user_data());
			reqs[i].set_user_data(NULL);
			req_ptrs[i] = &reqs[i];
		}
		process_reqs_on_io(req_ptrs.data(), num, notify_gcached_io);
		return num;
	}
};

class underlying_io_thread: public thread
{
	static const int MSG_BUF_SIZE = 16;
	static const int REQ_BUF_SIZE = 1024;

	thread_safe_FIFO_queue<message<io_request> > queue;
	message<io_request> msg_buf[MSG_BUF_SIZE];
	io_request req_buf[REQ_BUF_SIZE];
	underlying_io_proxy *io;
	size_t num_reqs;
public:
	underlying_io_thread(const std::vector<disk_io_thread *> &remotes,
			file_mapper *mapper, int node_id): thread(
			std::string("underlying_io_thread-") + itoa(node_id),
			node_id), queue("underlying_io_queue", node_id,
			MSG_BUF_SIZE, INT_MAX) {
		io = new underlying_io_proxy(remotes, mapper, this);
		num_reqs = 0;
		this->start();
	}

	void run() {
		while (!queue.is_empty()) {
			int ret = queue.fetch(msg_buf, MSG_BUF_SIZE);
			for (int i = 0; i < ret; i++) {
				num_reqs += msg_buf[i].get_num_objs();
				while (!msg_buf[i].is_empty()) {
					int num_reqs = msg_buf[i].get_next_objs(req_buf,
							REQ_BUF_SIZE);
					io->access(req_buf, num_reqs);
				}
			}
		}
		io->wait4complete(0);
	}

	void add_reqs(message<io_request> &requests) {
		int ret = queue.add(&requests, 1);
		assert(ret == 1);
		this->activate();
	}

	io_interface *get_io() const {
		return io;
	}

	void print_stat() {
		printf("%s issued %ld reqs to I/O threads\n",
				this->get_thread_name().c_str(), num_reqs);
	}
};

class req_stealer;
class flush_timer_task: public timer_task
{
	static const int FLUSH_TIMEOUT = 300000;
	req_stealer *io;
public:
	flush_timer_task(req_stealer *io): timer_task(FLUSH_TIMEOUT) {
		this->io = io;
	}

	void run();
};

class req_stealer: public io_interface
{
	underlying_io_thread *t;
	message<io_request> req_buf;
	slab_allocator *alloc;
	periodic_timer *flush_timer;

	atomic_integer access_guard;
public:
	req_stealer(underlying_io_thread *t, thread *curr_thread): io_interface(
			curr_thread) {
		flush_timer = new periodic_timer(curr_thread, new flush_timer_task(this));
		flush_timer->start();
		this->t = t;
		alloc = new slab_allocator(
				"io_msg_allocator-" + itoa(curr_thread->get_node_id()),
				PAGE_SIZE, PAGE_SIZE * 128, INT_MAX, curr_thread->get_node_id());
		message<io_request> tmp(alloc, false);
		req_buf = tmp;
	}

	virtual int get_file_id() const {
		// This function shouldn't be called.
		assert(0);
		return -1;
	}

	void access(io_request *reqs, int num, io_status *status = NULL) {
		access_guard.inc(1);
		int num_added = 0;
		for (int i = 0; i < num; i++) {
			// the requests are created by global_cached_io and user_data
			// isn't used by global_cached_io, so we should be safe here.
			assert(reqs[i].get_user_data() == NULL);
			reqs[i].set_user_data((void *) reqs[i].get_io());
			// We want I/O threads to notify the IO proxy of request completion
			// first.
			reqs[i].set_io(t->get_io());
		}
		while (num_added < num) {
			int ret = req_buf.add(reqs + num_added, num - num_added);
			num_added += ret;
			if (ret == 0) {
				t->add_reqs(req_buf);
				message<io_request> tmp(alloc, false);
				req_buf = tmp;
			}
		}
		access_guard.inc(1);
	}

	void flush_requests_timeout() {
		// The timeout signal is sent to the thread that uses the I/O
		// instance. When the access guard is even, the thread isn't using
		// req_buf, so we are safe to use it here.
		if (access_guard.get() % 2 == 0) {
			t->add_reqs(req_buf);
			message<io_request> tmp(alloc, false);
			req_buf = tmp;
		}
	}

	virtual void flush_requests() {
	}

	void print_state() {
		printf("req_stealer %d has %d buffered reqs in thread %d, timer %d is enabled? %d\n",
				get_io_idx(), req_buf.get_num_objs(), get_thread()->get_tid(),
				flush_timer->get_id(), flush_timer->is_enabled());
	}
};

void flush_timer_task::run()
{
	io->flush_requests_timeout();
}

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
	// This IO just dispatches any number of requests sent to it.
	// So we don't care how many requests are pending.
	set_max_num_pending_ios(INT_MAX);

	processing_thread_id = 0;
	pthread_key_create(&replier_key, NULL);

	for (int i = 0; i < MSG_BUF_SIZE; i++) {
		message<io_request> tmp(local_group->msg_allocator, false);
		local_msgs[i] = tmp;
	}
}

/* process the requests sent to this thread */
int node_cached_io::process_requests()
{
	if (processing_thread_id == 0)
		processing_thread_id = pthread_self();
	assert(processing_thread_id == pthread_self());
	int num_processed = 0;
	int num_access = 0;
	int num_received = 0;
	while (!request_queue->is_empty()) {
		int num = request_queue->fetch(tmp_msg_buf, MSG_BUF_SIZE);

		// We need to copy the message buffers to the local memory.
		for (int i = 0; i < num; i++) {
			tmp_msg_buf[i].copy_to(local_msgs[i]);
		}
		for (int i = 0; i < num; i++) {
			if (local_msgs[i].get_num_objs() > NUMA_REQ_BUF_SIZE - num_received) {
				global_cached_io::access(local_reqs, num_received);
				num_access++;
				num_processed += num_received;
				num_received = 0;
			}
			int ret = local_msgs[i].get_next_objs(local_reqs + num_received,
					NUMA_REQ_BUF_SIZE - num_received);
			num_received += ret;
			assert(local_msgs[i].is_empty());
		}
	}
	global_cached_io::access(local_reqs, num_received);
	num_access++;
	num_processed += num_received;
	if (num_processed > 0) {
		flush_requests();
		process_all_completed_requests();
	}
	else
		process_all_completed_requests();

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
	stack_array<io_reply> local_reply_buf(num);
	int num_remote = 0;
	int num_local = 0;

	for (int i = 0; i < num; i++) {
		page_id_t pg_id(reqs[i]->get_file_id(), reqs[i]->get_offset());
		int idx = cache_conf->page2cache(pg_id);
		if (idx != get_node_id())
			local_reply_buf[num_remote++] = io_reply(reqs[i], true, 0);
		else {
			assert(reqs[i]->get_io() == this);
			local_reqs[num_local++] = reqs[i];
		}
	}

	// The reply must be sent to the thread on a different node.
	int num_sent = reply_sender->send(local_reply_buf.data(), num_remote);
	// We use blocking queues here, so the send must succeed.
	assert(num_sent == num_remote);
	if (reply_queue->get_num_entries() > 0)
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
	struct thread_group *group;

	process_request_thread(struct thread_group *group): thread(
			std::string("process_request_thread-") + itoa(group->id), group->id) {
		this->io = NULL;
		this->group = group;
	}
public:
	static process_request_thread *create(struct thread_group *group) {
		void *addr = numa_alloc_onnode(sizeof(process_request_thread),
				group->id);
		return new(addr) process_request_thread(group);
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

	void init() {
		io = node_cached_io::create(
				new req_stealer(group->underlying_thread, this), group);
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

	void print_state() {
		io->print_state();
	}
};

void process_request_thread::run()
{
	io->process_requests();
}

void thread_group::print_state()
{
	for (unsigned i = 0; i < process_request_threads.size(); i++) {
		process_request_threads[i]->print_state();
	}
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
	curr_idx = (curr_idx + 1) % senders.size();
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
		const std::vector<disk_io_thread *> &io_threads,
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
		group.underlying_thread = new underlying_io_thread(io_threads,
				mapper, node_id);
		group.cache = cache_conf->create_cache_on_node(node_id,
				MAX_NUM_FLUSHES_PER_FILE * num_ssds / config->get_num_caches());
		group.msg_allocator = new slab_allocator(std::string("msg_allocator-")
				+ itoa(node_id), NUMA_MSG_SIZE,
				NUMA_MSG_SIZE * 128, INT_MAX, node_id);
		groups.insert(std::pair<int, struct thread_group>(node_id, group));

		struct thread_group *groupp = &groups[node_id];
		// Create processing threads.
		for (int i = 0; i < params.get_numa_num_process_threads(); i++) {
			process_request_thread *t = process_request_thread::create(groupp);
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

void part_io_process_table::print_stat(int num_user_threads)
{
	int tot_gcached_ios = num_user_threads + get_num_threads();
	std::set<int> node_ids = get_node_ids();

	printf("\n");
	for (std::set<int>::const_iterator it = node_ids.begin();
			it != node_ids.end(); it++) {
		int node_id = *it;
		const struct thread_group *group = get_thread_group(node_id);
		printf("cache on node %d\n", group->id);
		group->cache->print_stat();
	}
	printf("\n");
	for (std::set<int>::const_iterator it = node_ids.begin();
			it != node_ids.end(); it++) {
		int node_id = *it;
		const struct thread_group *group = get_thread_group(node_id);
		for (size_t i = 0; i < group->process_request_threads.size(); i++) {
			process_request_thread *thread = (process_request_thread *) group
				->process_request_threads[i];
			thread->get_io()->print_stat(tot_gcached_ios);
		}
	}
	printf("\n");
	for (std::set<int>::const_iterator it = node_ids.begin();
			it != node_ids.end(); it++) {
		int node_id = *it;
		const struct thread_group *group = get_thread_group(node_id);
		int tot_group_hits = 0;
		for (size_t i = 0; i < group->process_request_threads.size(); i++) {
			process_request_thread *thread = (process_request_thread *) group
				->process_request_threads[i];
			printf("group %d thread %ld gets %ld requests\n", group->id, i,
					thread->get_num_requests());
			tot_group_hits += thread->get_cache_hits();
		}
		printf("group %d gets %d hits\n", group->id, tot_group_hits);
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

class debug_process_request_threads: public debug_task
{
	part_io_process_table *table;
public:
	debug_process_request_threads(part_io_process_table *table) {
		this->table = table;
	}

	void run() {
		table->print_state();
	}
};

part_io_process_table *part_global_cached_io::open_file(
		const std::vector<disk_io_thread *> &io_threads,
		file_mapper *mapper, const cache_config *config)
{
	part_io_process_table *table = new part_io_process_table(io_threads,
			mapper, config);
	debug.register_task(new debug_process_request_threads(table));
	return table;
}

part_global_cached_io::part_global_cached_io(io_interface *underlying,
		part_io_process_table *table): io_interface(underlying->get_thread())
{
	int node_id = underlying->get_node_id();
	this->local_group = table->get_thread_group(node_id);

	this->underlying = new global_cached_io(underlying->get_thread(),
			new req_stealer(local_group->underlying_thread, underlying->get_thread()),
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
	reply_queue = msg_queue<io_reply>::create(node_id, "reply_queue", 
			// We don't want that adding replies is blocked, so we allow
			// the reply queue to be expanded to an arbitrary size.
			NUMA_REPLY_MSG_QUEUE_SIZE, INT_MAX / sizeof(io_reply), false);

	reply_sender = thread_safe_msg_sender<io_reply>::create(
			0, local_group->msg_allocator, reply_queue);

	req_senders = table->create_req_senders(node_id);

	for (int i = 0; i < MSG_BUF_SIZE; i++) {
		message<io_reply> tmp(local_group->msg_allocator, false);
		local_reply_msgs[i] = tmp;
	}
}

part_global_cached_io::~part_global_cached_io()
{
	msg_queue<io_reply>::destroy(reply_queue);
	thread_safe_msg_sender<io_reply>::destroy(reply_sender);
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
	io_request *reqs[1] = {&reply->get_request()};
	final_cb->invoke(reqs, 1);
	return ret;
}

int part_global_cached_io::process_replies()
{
	int num_processed = 0;
	message<io_reply> tmp_msgs[MSG_BUF_SIZE];
	// The reply buffer size should be able to contain more replies than
	// a reply message.
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
		page_id_t pg_id(requests[i].get_file_id(), requests[i].get_offset());
		int idx = cache_conf->page2cache(pg_id);
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
		reply_sender->flush_all();
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
	tot_remote_reads += remote_reads;
	seen_threads++;
	int tot_gcached_ios = nthreads + global_table->get_num_threads();
	underlying->print_stat(tot_gcached_ios);

	if (seen_threads == nthreads) {
		printf("part_global_cached_io: in total, there are %ld requests sent to the remote nodes\n",
				tot_remote_reads);
		global_table->print_stat(nthreads);
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
		page_id_t pg_id(get_file_id(), ROUND_PAGE(offset));
		page_id_t old_id;
		// We only preload data to the local cache.
		if (cache_conf->page2cache(pg_id) != get_node_id())
			continue;
		thread_safe_page *p = (thread_safe_page *) local_group->cache->search(
					pg_id, old_id);
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

void part_global_cached_io::print_state()
{
	printf("part global cached io %d has %d pending reqs, %ld remote pending reqs, %d local pending reqs %d replies in the queue\n",
			get_io_idx(), num_pending_ios(), sent_requests - processed_replies,
			underlying->num_pending_ios(), reply_queue->get_num_objs());
	for (std::tr1::unordered_map<int, group_request_sender *>::const_iterator it
			= req_senders.begin(); it != req_senders.end(); it++) {
		group_request_sender *sender = it->second;
		printf("\treq sender has %d buffered reqs\n", sender->get_num_remaining());
	}
	printf("\treply sender has %d buffered reqs\n", reply_sender->get_num_remaining());
	underlying->print_state();
}
