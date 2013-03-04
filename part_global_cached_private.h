#ifndef __PART_GLOBAL_CACHED_PRIVATE_H__
#define __PART_GLOBAL_CACHED_PRIVATE_H__

#include <errno.h>

#include <tr1/unordered_map>

#include "access_mapper.h"
#include "parameters.h"
#include "messaging.h"
#include "garbage_collection.h"
#include "global_cached_private.h"
#include "thread.h"

class part_global_cached_io;

struct thread_group
{
	int id;
	int nthreads;
	part_global_cached_io **ios;
	page_cache *cache;
	std::vector<thread *> process_request_threads;
	std::vector<thread *> process_reply_threads;
	
	blocking_FIFO_queue<io_request> *request_queue;
	blocking_FIFO_queue<io_reply> *reply_queue;
};

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

	int num_groups;
	int group_idx;
	struct thread_group *local_group;

	/* the size of the cache associated to the thread. */
	long cache_size;
	int cache_type;

	volatile int finished_threads;

	/*
	 * there is a sender for each node.
	 */
	std::tr1::unordered_map<int, msg_sender<io_request> *> req_senders;
	std::tr1::unordered_map<int, msg_sender<io_reply> *> reply_senders;

	access_mapper *mapper;
	int hash_req(io_request *req)
	{
		return mapper->map(req->get_offset());
	}

	long processed_requests;

	long remote_reads;

	int thread_id;

	// It's the callback from the user.
	callback *final_cb;
	// It's the callback used by global_cached_io.
	callback *my_cb;

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
		// TODO delete all senders
	}

	int init();

	part_global_cached_io(int num_groups, io_interface *underlying,
			int idx, long cache_size, int cache_type, access_mapper *mapper);

	virtual page_cache *get_global_cache() {
		return groups[group_idx].cache;
	}

	virtual bool set_callback(callback *cb) {
		this->final_cb = cb;
		return true;
	}

	virtual callback *get_callback() {
		return my_cb;
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
