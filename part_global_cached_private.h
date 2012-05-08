#ifndef __PART_GLOBAL_CACHED_PRIVATE_H__
#define __PART_GLOBAL_CACHED_PRIVATE_H__

#include <errno.h>

#define BUF_SIZE 1000
#define REQ_QUEUE_SIZE 100000
#define REPLY_QUEUE_SIZE REQ_QUEUE_SIZE

#include "messaging.h"
#include "garbage_collection.h"
#include "global_cached_private.h"

class part_global_cached_private;

struct thread_group
{
	int id;
	int nthreads;
	part_global_cached_private **threads;
	page_cache *cache;
};

class part_global_cached_private: public global_cached_private
{
	memory_manager *manager;
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
	
	bulk_queue<io_request> *request_queue;
	bulk_queue<io_reply> *reply_queue;

	volatile int finished_threads;

	/*
	 * there is a sender for each node.
	 */
	msg_sender<io_request> **req_senders;
	msg_sender<io_reply> **reply_senders;

	int hash_req(io_request *req)
	{
		/* the size of each group needs to access */
		long size = get_size() / num_groups;
		return req->get_offset() / size;
	}

	long processed_requests;

	long remote_reads;

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

	part_global_cached_private *id2thread(int thread_id) {
		return groups[group_id(thread_id, num_groups)].threads[thread_idx(thread_id, num_groups)];
	}

	~part_global_cached_private() {
		// TODO delete all senders
		delete request_queue;
		delete reply_queue;
	}

	int thread_init();

	part_global_cached_private(int num_groups, read_private *underlying,
			int idx, long cache_size, int entry_size, int cache_type,
			memory_manager *manager);

	virtual page_cache *get_global_cache() {
		return groups[group_idx].cache;
	}

	int reply(io_request *requests, io_reply *replies, int num);

	void distribute_reqs(io_request *requests, int num);

	int process_requests(int max_nreqs);

	int process_replies(int max_nreplies);

	int process_reply(io_reply *reply);

	ssize_t access(io_request *requests, int num, int access_method);

	void cleanup();

	int get_group_id() {
		return group_idx;
	}

	bool support_bulk() {
		return true;
	}

#ifdef STATISTICS
	virtual void print_stat() {
		global_cached_private::print_stat();
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
