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

#define REPLY_BUF_SIZE 1000
#define REQ_BUF_SIZE 1000

class part_global_cached_io;

struct thread_group
{
	int id;
	int nthreads;
	part_global_cached_io **ios;
	page_cache *cache;
	std::vector<thread *> process_request_threads;
	
	blocking_FIFO_queue<io_request> *request_queue;

	std::tr1::unordered_map<part_global_cached_io *,
		thread_safe_msg_sender<io_reply> *> reply_senders;
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
	static int nthreads;

	/*
	 * These counts are used for stopping protocols.
	 */
	static atomic_integer num_finish_issuing_threads;
	static atomic_integer num_finished_threads;

	int group_idx;
	struct thread_group *local_group;

	const cache_config *cache_conf;
	blocking_FIFO_queue<io_reply> *reply_queue;
	io_reply local_reply_buf[REPLY_BUF_SIZE];
	io_request local_req_buf[REQ_BUF_SIZE];

	/*
	 * there is a sender for each node.
	 */
	// group id <-> msg sender
	std::tr1::unordered_map<int, request_sender *> req_senders;

	io_interface *underlying;

	long processed_requests;
	long sent_requests;
	long processed_replies;

	long remote_reads;

	int thread_id;

	// It's the callback from the user.
	callback *final_cb;

	int process_replies();
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

	static void set_num_threads(const int num) {
		nthreads = num;
	}

	~part_global_cached_io() {
		// TODO delete all senders
	}

	int init();

	part_global_cached_io(int num_groups, io_interface *underlying,
			int idx, const cache_config *config);

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

	void access(io_request *requests, int num, io_status *status);
	io_status access(char *, off_t, ssize_t, int) {
		return IO_UNSUPPORTED;
	}

	void cleanup();
	int preload(off_t start, long size);

	int get_group_id() {
		return group_idx;
	}

	bool support_aio() {
		return true;
	}

	blocking_FIFO_queue<io_reply> *get_reply_queue() {
		return reply_queue;
	}

	friend class node_cached_io;

#ifdef STATISTICS
	virtual void print_stat();
#endif
};

#endif
