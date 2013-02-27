#ifndef __PART_GLOBAL_CACHED_PRIVATE_H__
#define __PART_GLOBAL_CACHED_PRIVATE_H__

#include <errno.h>

#include "parameters.h"
#include "messaging.h"
#include "garbage_collection.h"
#include "global_cached_private.h"

class part_global_cached_io;

struct thread_group
{
	int id;
	int nthreads;
	part_global_cached_io **ios;
	page_cache *cache;
};

class part_global_cached_io: public global_cached_io
{
	static thread_group *groups;
	/* this mutex just for helping initialize cache. */
	static pthread_mutex_t init_mutex;
	/* indicates the number of threads that finish initialization. */
	static int num_finish_init;
	/* used for thread initialization. */
	static pthread_mutex_t wait_mutex;
	static pthread_cond_t cond;

	int num_groups;
	int group_idx;

	/* the size of the cache associated to the thread. */
	long cache_size;
	int cache_type;
	
	thread_safe_FIFO_queue<io_request> *request_queue;
	thread_safe_FIFO_queue<io_reply> *reply_queue;

	volatile int finished_threads;

	/*
	 * there is a sender for each node.
	 */
	msg_sender<io_request> **req_senders;
	msg_sender<io_reply> **reply_senders;

	int hash_req(io_request *req)
	{
		return req->get_offset() % num_groups;
	}

	long processed_requests;

	long remote_reads;

	int thread_id;

	callback *cb;

public:
	static inline int group_id(int thread_id, int num_groups) {
		int remaining = nthreads % num_groups;
		int group_size = nthreads / num_groups;
		if (thread_id <= remaining * (group_size + 1))
			return thread_id / (group_size + 1);
		else
			return (thread_id - remaining * (group_size + 1)) / group_size + remaining;
	}

	/* get the location of a thread in the group. */
	static inline int thread_idx(int thread_id, int num_groups) {
		int remaining = nthreads % num_groups;
		int group_size = nthreads / num_groups;
		if (thread_id <= remaining * (group_size + 1))
			return thread_id % (group_size + 1);
		else
			return (thread_id - remaining * (group_size + 1)) % group_size;
	}

	part_global_cached_io *id2thread(int thread_id) {
		return groups[group_id(thread_id, num_groups)].ios[thread_idx(thread_id, num_groups)];
	}

	~part_global_cached_io() {
		// TODO delete all senders
		delete request_queue;
		delete reply_queue;
	}

	int init();

	part_global_cached_io(int num_groups, io_interface *underlying,
			int idx, long cache_size, int cache_type);

	virtual page_cache *get_global_cache() {
		return groups[group_idx].cache;
	}

	virtual bool set_callback(callback *cb);

	virtual callback *get_callback() {
		return cb;
	}

	int reply(io_request *requests, io_reply *replies, int num);

	int distribute_reqs(io_request *requests, int num);

	int process_requests(int max_nreqs);

	int process_replies(int max_nreplies);

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

#ifdef STATISTICS
	virtual void print_stat() {
		global_cached_io::print_stat();
		static long tot_remote_reads = 0;
		static int seen_threads = 0;
		tot_remote_reads += remote_reads;
		seen_threads++;
		if (seen_threads == nthreads) {
			printf("there are %ld requests sent to the remote nodes\n",
					tot_remote_reads);
		}
	}
#endif
};

#endif
