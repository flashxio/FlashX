#ifndef __AIO_PRIVATE_H__
#define __AIO_PRIVATE_H__

#include <deque>
#include <tr1/unordered_map>

#include "wpaio.h"
#include "read_private.h"
#include "thread.h"
#include "container.h"

#ifdef ENABLE_AIO

void aio_callback(io_context_t, struct iocb*, void *, long, long);

struct thread_callback_s;
class aio_complete_thread;

class aio_complete_queue
{
	blocking_FIFO_queue<thread_callback_s *> queue;
public:
	aio_complete_queue(): queue("aio_completes", AIO_DEPTH_PER_FILE,
			// This max size seems enough to keep all requests.
			AIO_DEPTH_PER_FILE * 200) {
	}

	blocking_FIFO_queue<thread_callback_s *> *get_queue() {
		return &queue;
	}

	int process(int max_num, bool blocking);
};

class aio_complete_sender: public simple_msg_sender<thread_callback_s *>
{
public:
	aio_complete_sender(aio_complete_queue *queue): simple_msg_sender<thread_callback_s *>(
			queue->get_queue(), AIO_DEPTH_PER_FILE) {
	}
};

class async_io: public buffered_io
{
	int buf_idx;
	struct aio_ctx *ctx;
	callback *cb;
	const int AIO_DEPTH;
	// This is for allocating memory in the case that the memory given
	// by the user isn't in the local NUMA node.
	slab_allocator allocator;
	obj_allocator<thread_callback_s> cb_allocator;
	std::tr1::unordered_map<int, aio_complete_sender *> complete_senders;
	std::tr1::unordered_map<int, fifo_queue<thread_callback_s *> *> remote_tcbs;

	int num_iowait;
	int num_completed_reqs;
	int num_local_alloc;

	struct iocb *construct_req(io_request &io_req, callback_t cb_func);
public:
	/**
	 * @size: the size of data to be accessed in all files
	 * @aio_depth_per_file
	 * @node_id: the NUMA node where the disks to be read are connected to.
	 */
	async_io(const logical_file_partition &partition,
			const std::tr1::unordered_map<int, aio_complete_thread *> &complete_threads,
			long size, int aio_depth_per_file, int node_id);

	virtual ~async_io();

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method) {
		return -1;
	}

	ssize_t access(io_request *requests, int num);

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

	virtual void cleanup();

	void return_cb(thread_callback_s *tcbs[], int num);

	int num_available_IO_slots() const {
		return max_io_slot(ctx);
	}

	int num_pending_IOs() const {
		return AIO_DEPTH - max_io_slot(ctx);
	}

	void wait4complete() {
		io_wait(ctx, NULL, 1);
	}

	int get_num_iowait() const {
		return num_iowait;
	}

	int get_num_completed_reqs() const {
		return num_completed_reqs;
	}

	int get_num_local_alloc() const {
		return num_local_alloc;
	}
};

class aio_complete_thread: public thread
{
	aio_complete_queue queue;
	int num_completed_reqs;
public:
	aio_complete_thread(int node_id): thread("aio_complete_thread",
			node_id) {
		num_completed_reqs = 0;
		start();
	}

	void run();

	int get_num_completed_reqs() const {
		return num_completed_reqs;
	}

	aio_complete_queue *get_queue() {
		return &queue;
	}
};

#endif

#endif
