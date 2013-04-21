#ifndef __THREAD_PRIVATE_H__
#define __THREAD_PRIVATE_H__

#include "rand_buf.h"
#include "garbage_collection.h"
#include "messaging.h"

extern int buf_type;

enum {
	SINGLE_LARGE_BUF,
	SINGLE_SMALL_BUF,
	MULTI_BUF,
};

class callback
{
public:
	virtual int invoke(io_request *reqs[], int num) = 0;
};

/**
 * The interface for all IO classes.
 */
class io_interface
{
	int node_id;
public:
	io_interface(int node_id) {
		this->node_id = node_id;
	}

	virtual ~io_interface() { }

	/* When a thread begins, this method will be called. */
	virtual int init() {
		return 0;
	}

	/**
	 * set the callback if the class supports the asynchronous fashion.
	 * If the class doesn't support async IO, return false.
	 */
	virtual bool set_callback(callback *cb) {
		return false;
	}

	virtual callback *get_callback() {
		return NULL;
	}

	virtual bool support_aio() {
		return false;
	}

	/**
	 * The total size accessible with this IO interface.
	 */
	virtual ssize_t get_size() const {
		return 0;
	}

	/**
	 * The size of data on the local node.
	 */
	virtual ssize_t get_local_size() const {
		return 0;
	}

	virtual void cleanup() {
	}

	/**
	 * The asynchronous IO interface
	 */
	virtual ssize_t access(io_request *requests, int num) {
		return -1;
	}
	virtual void flush_requests() {
	}

	/**
	 * The synchronous IO interface.
	 */
	virtual ssize_t access(char *, off_t, ssize_t, int) {
		return -1;
	}

	virtual void print_stat() {
	}

	virtual io_interface *clone() const {
		return NULL;
	}

	int get_node_id() const {
		return node_id;
	}
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
	workload_gen *gen;
	rand_buf *buf;
	io_interface *io;

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

	virtual int thread_init();

	virtual int get_group_id() {
		return 0;
	}

	thread_private(int idx, int entry_size, io_interface *io) {
		this->cb = NULL;
		this->idx = idx;
		buf = NULL;
		this->gen = NULL;
		this->io = io;
		read_bytes = 0;
		num_accesses = 0;
		num_sampling = 0;
		tot_num_pending = 0;
		max_num_pending = 0;
	}

	io_interface *get_io() {
		return io;
	}

	int attach2cpu();

	int get_idx() {
		return idx;
	}

	void set_workload(workload_gen *gen) {
		this->gen = gen;
	}

	int run();

	ssize_t get_read_bytes();

	int start_thread();

	int wait_thread_end();

	virtual void print_stat() {
#ifdef STATISTICS
		io->print_stat();
		int avg_num_pending = 0;
		if (num_sampling > 0)
			avg_num_pending = tot_num_pending / num_sampling;
		printf("access %ld bytes in %ld accesses (%d completes), avg pending: %d, max pending: %d\n",
				get_read_bytes(), num_accesses, num_completes.get(),
				avg_num_pending, max_num_pending);
#endif
		printf("thread %d: start at %f seconds, takes %f seconds, access %ld bytes in %ld accesses\n", idx,
				time_diff(global_start, start_time), time_diff(start_time, end_time),
				get_read_bytes(), num_accesses);
	}
};

#endif
