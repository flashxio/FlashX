#ifndef __THREAD_PRIVATE_H__
#define __THREAD_PRIVATE_H__

#include "rand_buf.h"

#define NUM_PAGES 16384

extern int nthreads;
extern struct timeval global_start;

inline float time_diff(struct timeval time1, struct timeval time2)
{
	return time2.tv_sec - time1.tv_sec
			+ ((float)(time2.tv_usec - time1.tv_usec))/1000000;
}

/* this data structure stores the thread-private info. */
class thread_private
{
public:
#ifdef USE_PROCESS
	pid_t id;
#else
	pthread_t id;
#endif
	/* the location in the thread descriptor array. */
	int idx;
	rand_buf buf;
	workload_gen *gen;
	ssize_t read_bytes;
	struct timeval start_time;
	struct timeval end_time;

#ifdef STATISTICS
	int cache_hits;
#endif

	virtual ssize_t access(char *, off_t, ssize_t, int) = 0;
	virtual int thread_init() = 0;

	thread_private(int idx, int entry_size): buf(NUM_PAGES / nthreads * PAGE_SIZE, entry_size) {
		this->idx = idx;
		read_bytes = 0;
#ifdef STATISTICS
		cache_hits = 0;
#endif
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
