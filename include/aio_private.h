#ifndef __AIO_PRIVATE_H__
#define __AIO_PRIVATE_H__

#include <deque>
#include <tr1/unordered_map>

#include "wpaio.h"
#include "io_interface.h"
#include "thread.h"
#include "container.h"
#include "io_request.h"

void aio_callback(io_context_t, struct iocb*, void *, long, long);

struct thread_callback_s;

class buffered_io;
class logical_file_partition;
class callback_allocator;

class async_io: public io_interface
{
	int buf_idx;
	struct aio_ctx *ctx;
	callback *cb;
	const int AIO_DEPTH;
	callback_allocator *cb_allocator;

	int num_iowait;
	int num_completed_reqs;

	// file id <-> buffered io
	std::tr1::unordered_map<int, buffered_io *> open_files;
	buffered_io *default_io;

	struct iocb *construct_req(io_request &io_req, callback_t cb_func);
public:
	/**
	 * @aio_depth_per_file
	 * @node_id: the NUMA node where the disks to be read are connected to.
	 */
	async_io(const logical_file_partition &partition,
			int aio_depth_per_file, int node_id);

	virtual ~async_io();

	virtual io_status access(char *, off_t, ssize_t, int) {
		return IO_UNSUPPORTED;
	}
	virtual void access(io_request *requests, int num, io_status *status = NULL);

	bool set_callback(callback *cb) {
		this->cb = cb;
		return true;
	}

	callback *get_callback() {
		return cb;
	}

	bool support_aio() {
		return true;
	}

	int get_file_id() const;

	virtual void cleanup();

	void return_cb(thread_callback_s *tcbs[], int num);

	int num_available_IO_slots() const {
		return max_io_slot(ctx);
	}

	virtual int num_pending_ios() const {
		return AIO_DEPTH - max_io_slot(ctx);
	}

	virtual void notify_completion(io_request *reqs[], int num);
	int wait4complete(int num) {
		return io_wait(ctx, NULL, num);
	}
	virtual int get_max_num_pending_ios() const {
		return AIO_DEPTH;
	}
	virtual void set_max_num_pending_ios(int max) {
	}

	int get_num_iowait() const {
		return num_iowait;
	}

	int get_num_completed_reqs() const {
		return num_completed_reqs;
	}

	virtual void flush_requests();

	// These two interfaces allow users to open and close more files.
	
	/*
	 * It opens a virtual file.
	 * Actually, it opens physical files on the underlying filesystems
	 * within the partition of the virtual file, managed by the IO interface.
	 */
	int open_file(const logical_file_partition &partition);
	int close_file(int file_id);
};

void init_aio(std::vector<int> node_ids);

#endif
