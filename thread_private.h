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

/**
 * The interface for all IO classes.
 */
class io_interface
{
public:
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
};

class cleanup_callback;

/* this data structure stores the thread-private info. */
class thread_private
{
	int entry_size;
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

	cleanup_callback *cb;
	
	struct timeval start_time, end_time;

public:
	thread_private() {
		cb = NULL;
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
		this->idx = idx;
		this->entry_size = entry_size;
		buf = NULL;
		this->gen = NULL;
		this->io = io;
		read_bytes = 0;
	}

	io_interface *get_io() {
		return io;
	}

	int attach2cpu();

	int get_entry_size() {
		return entry_size;
	}

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
		printf("read %ld bytes, start at %f seconds, takes %f seconds\n",
				get_read_bytes(), time_diff(global_start, start_time),
				time_diff(start_time, end_time));
	}
#endif
};

#endif
