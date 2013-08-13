#ifndef __THREAD_PRIVATE_H__
#define __THREAD_PRIVATE_H__

#include "workload.h"
#include "rand_buf.h"
#include "garbage_collection.h"
#include "messaging.h"
#include "io_interface.h"

extern int buf_type;

enum {
	SINGLE_LARGE_BUF,
	SINGLE_SMALL_BUF,
	MULTI_BUF,
};

class cleanup_callback;

/* this data structure stores the thread-private info. */
class thread_private
{
#ifdef USE_PROCESS
	pid_t id;
#else
	pthread_t id;
#endif
	/* the location in the thread descriptor array. */
	int idx;
	int node_id;
	workload_gen *gen;
	rand_buf *buf;
	io_interface *io;
	file_io_factory *factory;

	ssize_t read_bytes;
	long num_accesses;

	/* compute average number of pending IOs. */
	long tot_num_pending;
	long num_sampling;
	int max_num_pending;

	cleanup_callback *cb;
	
	struct timeval start_time, end_time;

#ifdef STATISTICS
public:
	atomic_integer num_completes;
	atomic_integer num_pending;
#endif

public:
	thread_private() {
		cb = NULL;
		num_sampling = 0;
		tot_num_pending = 0;
		max_num_pending = 0;
	}

	~thread_private() {
		if (buf)
			delete buf;
	}

	int thread_init();

	int get_node_id() {
		return node_id;
	}

	thread_private(int node_id, int idx, int entry_size, file_io_factory *factory,
			workload_gen *gen) {
		this->cb = NULL;
		this->node_id = node_id;
		this->idx = idx;
		buf = NULL;
		this->gen = gen;
		this->io = NULL;
		this->factory = factory;
		read_bytes = 0;
		num_accesses = 0;
		num_sampling = 0;
		tot_num_pending = 0;
		max_num_pending = 0;
	}

	int attach2cpu();

	int run();

	ssize_t get_read_bytes();

	int start_thread();

	int wait_thread_end();

	void print_stat() {
#ifdef STATISTICS
		extern int nthreads;
		io->print_stat(nthreads);
		int avg_num_pending = 0;
		if (num_sampling > 0)
			avg_num_pending = tot_num_pending / num_sampling;
		printf("access %ld bytes in %ld accesses (%d completes), avg pending: %d, max pending: %d\n",
				get_read_bytes(), num_accesses, num_completes.get(),
				avg_num_pending, max_num_pending);
#endif
		extern struct timeval global_start;
		printf("thread %d: start at %f seconds, takes %f seconds, access %ld bytes in %ld accesses\n", idx,
				time_diff(global_start, start_time), time_diff(start_time, end_time),
				get_read_bytes(), num_accesses);
	}
};

#endif
