#ifndef __PART_GLOBAL_CACHED_PRIVATE_H__
#define __PART_GLOBAL_CACHED_PRIVATE_H__

#include <errno.h>

#define BUF_SIZE 1000
#define REQ_QUEUE_SIZE 100000
#define REPLY_QUEUE_SIZE REQ_QUEUE_SIZE

#include "garbage_collection.h"
#include "global_cached_private.h"

class io_reply
{
	char *buf;
	off_t offset;
	ssize_t size: 32;
	int success: 1;
	int status: 16;
	int access_method: 1;
	void init(char *buf, off_t off, ssize_t size, int success,
			int status, int access_method) {
		this->buf = buf;
		this->offset = off;
		this->size = size;
		this->success = success;
		this->status = status;
		this->access_method = access_method;
	}
public:
	io_reply() {
		init(NULL, 0, 0, 0, 0, READ);
	}

	io_reply(io_request *req, int success, int status) {
		init(req->get_buf(), req->get_offset(), req->get_size(),
					success, status, req->get_access_method());
	}

	int get_status() {
		return status;
	}

	bool is_success() {
		return success;
	}

	char *get_buf() {
		return buf;
	}

	off_t get_offset() {
		return offset;
	}

	ssize_t get_size() {
		return size;
	}

	int get_access_method() {
		return access_method;
	}
};

class part_global_cached_private;

struct thread_group
{
	int id;
	int nthreads;
	part_global_cached_private **threads;
	page_cache *cache;
};

inline int min(int v1, int v2)
{
	return v1 > v2 ? v2 : v1;
}

/**
 * this is a thread-safe FIFO queue.
 * It supports bulk operations.
 */
template<class T>
class bulk_queue
{
	T *buf;
	volatile int size;
	int start;
	int num_entries;
	pthread_spinlock_t _lock;
public:
	bulk_queue(int size) {
		buf = new T[size];
		this->size = size;
		start = 0;
		num_entries = 0;
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

	~bulk_queue() {
		pthread_spin_destroy(&_lock);
		delete [] buf;
	}

	int fetch(T *entries, int num);

	int add(T *entries, int num);

	int get_num_entries() {
		return num_entries;
	}

	bool is_full() {
		return num_entries == size;
	}

	bool is_empty() {
		return num_entries == 0;
	}
};

class part_global_cached_private: public global_cached_private
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
	
	bulk_queue<io_request> *request_queue;
	bulk_queue<io_reply> *reply_queue;

	volatile int finished_threads;

	/* 
	 * we have request and reply buffer to distribute requests and replies
	 * in case we can't send them all.
	 */
	/* there is a request buffer for each group it sends to */
	io_request **thread_reqs;
	int *nreqs;
	/* there is a reply buffer for each thread it sends to */
	io_reply **thread_replies;
	int *nreplies;

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
		for (int i = 0; i < num_groups; i++) {
			numa_free(thread_reqs[i], sizeof(io_request) * BUF_SIZE);
			numa_free(thread_replies[i], sizeof(io_reply) * BUF_SIZE);
		}
		numa_free(thread_reqs, sizeof(thread_reqs[0]) * num_groups);
		numa_free(thread_replies, sizeof(thread_replies[0]) * num_groups);
		delete request_queue;
		delete reply_queue;
	}

	int thread_init();

	part_global_cached_private(int num_groups, const char *names[],
			int num, long size, int idx, long cache_size, int entry_size,
			int cache_type);

	virtual page_cache *get_global_cache() {
		return groups[group_idx].cache;
	}

	io_request *send(int node_id, io_request *reqs, int num);

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
