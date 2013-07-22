#ifndef __IO_INTERFACE_H__
#define __IO_INTERFACE_H__

#include <stdlib.h>

#include <vector>

#include "exception.h"
#include "common.h"

class io_request;

class callback
{
public:
	virtual int invoke(io_request *reqs[], int num) = 0;
};

class io_status
{
	long status: 8;
	long priv_data: 56;
public:
	io_status() {
		status = 0;
		priv_data = 0;
	}

	io_status(int status) {
		this->status = status;
		priv_data = 0;
	}

	void set_priv_data(long data) {
		priv_data = data;
	}

	long get_priv_data() const {
		return priv_data;
	}

	io_status &operator=(int status) {
		this->status = status;
		return *this;
	}

	bool operator==(int status) {
		return this->status == status;
	}
};

enum
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
	 * This method waits for at least the specified number of requests currently
	 * being sent by the access method to complete.
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

io_interface *allocate_io();
void release_io(io_interface *io);
io_interface *get_io(int idx);
int get_num_ios();

enum {
	READ_ACCESS,
	DIRECT_ACCESS,
#ifdef ENABLE_AIO
	AIO_ACCESS,
#endif
	REMOTE_ACCESS,
	GLOBAL_CACHE_ACCESS,
	PART_GLOBAL_ACCESS,
};

class RAID_config;
class cache_config;

std::vector<io_interface *> create_ios(const RAID_config &raid_conf,
		cache_config *cache_conf, const std::vector<int> &node_id_array,
		int nthreads, int access_option, long size, bool preload);

#endif
