#ifndef __IO_INTERFACE_H__
#define __IO_INTERFACE_H__

#include <stdlib.h>

#include "exception.h"

class io_request;

class callback
{
public:
	virtual int invoke(io_request *reqs[], int num) = 0;
};

enum io_status
{
	IO_OK,
	IO_PENDING = -1,
	IO_FAIL = -2,
	IO_UNSUPPORTED = -3,
};

/**
 * The interface for all IO classes.
 */
class io_interface
{
	int node_id;
	// This is an index for locating this IO object in a global table.
	int io_idx;
public:
	io_interface(int node_id) {
		this->node_id = node_id;
		this->io_idx = -1;
	}

	virtual ~io_interface() { }

	/* When a thread begins, this method will be called. */
	virtual int init() {
		return 0;
	}

	void set_io_idx(int idx) {
		io_idx = idx;
	}

	int get_io_idx() const {
		return io_idx;
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

	/**
	 * The main interface to send requests.
	 */
	virtual void access(io_request *requests, int num, io_status *status = NULL) {
		throw unsupported_exception();
	}
	/**
	 * When requests are passed to the access method, an IO layer may buffer
	 * the requests. This method guarantees that all requests are flushed to
	 * the underlying devices.
	 */
	virtual void flush_requests() {
	}
	/**
	 * This method waits for all requests currently being sent by the access
	 * method to complete.
	 */
	virtual void wait4complete() {
	}

	/**
	 * The synchronous IO interface.
	 */
	virtual io_status access(char *, off_t, ssize_t, int) {
		return IO_UNSUPPORTED;
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

void register_io(io_interface *io);
io_interface *get_io(int idx);

#endif
