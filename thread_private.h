#ifndef __THREAD_PRIVATE_H__
#define __THREAD_PRIVATE_H__

#include "rand_buf.h"
#include "garbage_collection.h"
#include "messaging.h"

#define NUM_PAGES (4096 * nthreads)

class callback
{
public:
	virtual int invoke(io_request *) = 0;
};

/* this data structure stores the thread-private info. */
class thread_private
{
	int entry_size;
public:
#ifdef USE_PROCESS
	pid_t id;
#else
	pthread_t id;
#endif
	/* the location in the thread descriptor array. */
	int idx;
	workload_gen *gen;
	ssize_t read_bytes;
	struct timeval start_time;
	struct timeval end_time;
	rand_buf *buf;
	callback *cb;

#ifdef STATISTICS
	int cache_hits;
#endif

	virtual ssize_t access(char *, off_t, ssize_t, int) = 0;
	virtual int thread_init() {
		buf = new rand_buf(NUM_PAGES / (nthreads / NUM_NODES) * PAGE_SIZE,
				get_entry_size());
		return 0;
	}

	/* by default, the base class doesn't support bulk operations. */
	virtual bool support_bulk() {
		return false;
	}

	virtual void cleanup() {
	}

	virtual ssize_t access(io_request *requests, int num, int access_method) {
		return -1;
	}

	virtual int get_group_id() {
		return 0;
	}

	thread_private(int idx, int entry_size) {
		this->idx = idx;
		read_bytes = 0;
		this->entry_size = entry_size;
		buf = NULL;
#ifdef STATISTICS
		cache_hits = 0;
#endif
	}

	int get_entry_size() {
		return entry_size;
	}

#ifdef STATISTICS
	virtual void print_stat() {
		static int tot_hits = 0;
		static int seen_threads = 0;
		seen_threads++;
		tot_hits += cache_hits;
		printf("read %ld bytes, start at %f seconds, takes %f seconds\n",
				read_bytes, time_diff(global_start, start_time),
				time_diff(start_time, end_time));
		if (seen_threads == nthreads)
			printf("there are %d cache hits\n", tot_hits);
	}
#endif
};

#endif
