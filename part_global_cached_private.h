#ifndef __PART_GLOBAL_CACHED_PRIVATE_H__
#define __PART_GLOBAL_CACHED_PRIVATE_H__

#include <errno.h>

#include <tr1/unordered_map>

#include "parameters.h"
#include "messaging.h"
#include "garbage_collection.h"
#include "global_cached_private.h"
#include "thread.h"
#include "file_mapper.h"

class part_global_cached_io;

struct thread_group
{
	int id;
	int nthreads;
	part_global_cached_io **ios;
	page_cache *cache;
	std::vector<thread *> process_request_threads;
	
	blocking_FIFO_queue<io_request> *request_queue;
	blocking_FIFO_queue<io_reply> *reply_queue;

	thread *reply_processor;
	std::tr1::unordered_map<int, thread_safe_msg_sender<io_reply> *> reply_senders;
};

/**
 * This provides interface for application threads issue IO requests
 * in the NUMA architecture. If the requests are sent to disks connected
 * to the local processor, they are processed directly. No message passing
 * is needed.
 */
class part_global_cached_io: public global_cached_io
{
	static std::tr1::unordered_map<int, struct thread_group> groups;
	/* this mutex just for helping initialize cache. */
	static pthread_mutex_t init_mutex;
	/* indicates the number of threads that finish initialization. */
	static int num_finish_init;
	/* used for thread initialization. */
	static pthread_mutex_t wait_mutex;
	static pthread_cond_t cond;

	/*
	 * These counts are used for stopping protocols.
	 */
	static atomic_integer num_finish_issuing_threads;
	static atomic_integer num_finished_threads;

	int group_idx;
	struct thread_group *local_group;

	/* the size of the cache associated to the thread. */
	long cache_size;
	int cache_type;

	/*
	 * there is a sender for each node.
	 */
	// group id <-> msg sender
	std::tr1::unordered_map<int, request_sender *> req_senders;

	file_mapper *mapper;
	int hash_req(io_request *req)
	{
		int idx = mapper->map2file(req->get_offset() / PAGE_SIZE);
		return mapper->get_file_node_id(idx);
	}

	io_interface *underlying;

	atomic_integer processed_requests;
	atomic_integer processed_replies;

	long remote_reads;

	int thread_id;

	// It's the callback from the user.
	callback *final_cb;

	int process_replies(int num);
	int distribute_reqs(io_request *requests, int num);

public:
	/* get the location of a thread in the group. */
	static inline int thread_idx(int thread_id, int num_groups) {
		int remaining = nthreads % num_groups;
		int group_size = nthreads / num_groups;
		if (thread_id <= remaining * (group_size + 1))
			return thread_id % (group_size + 1);
		else
			return (thread_id - remaining * (group_size + 1)) % group_size;
	}

	~part_global_cached_io() {
		delete mapper;
		// TODO delete all senders
	}

	int init();

	part_global_cached_io(int num_groups, io_interface *underlying,
			int idx, long cache_size, int cache_type, file_mapper *mapper);

	virtual page_cache *get_global_cache() {
		return groups[group_idx].cache;
	}

	virtual bool set_callback(callback *cb) {
		this->final_cb = cb;
		return true;
	}

	virtual callback *get_callback() {
		return final_cb;
	}

	int reply(io_request *requests, io_reply *replies, int num);

	int process_requests(int max_nreqs);

	int process_reply(io_reply *reply);

	ssize_t access(io_request *requests, int num);
	ssize_t access(char *, off_t, ssize_t, int) {
		return -1;
	}

	void cleanup();

	int get_group_id() {
		return group_idx;
	}

	bool support_aio() {
		return true;
	}
	friend class node_cached_io;

#ifdef STATISTICS
	virtual void print_stat();
#endif
};

#endif
