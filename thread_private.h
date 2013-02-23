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
	virtual int invoke(io_request *) = 0;
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

	virtual ssize_t get_size() {
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
	int node_id;

	ssize_t read_bytes;
	long num_accesses;

	cleanup_callback *cb;
	
	struct timeval start_time, end_time;

public:
	thread_private(int node_id) {
		cb = NULL;
		this->node_id = node_id;
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

#ifdef STATISTICS
	virtual void print_stat() {
		io->print_stat();
		printf("access %ld bytes in %ld accesses, start at %f seconds, takes %f seconds\n",
				get_read_bytes(), num_accesses, time_diff(global_start, start_time),
				time_diff(start_time, end_time));
	}
#endif
};

#endif
