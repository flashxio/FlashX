#ifndef __AIO_PRIVATE_H__
#define __AIO_PRIVATE_H__

#include <deque>
#include <tr1/unordered_map>

#include "wpaio.h"
#include "read_private.h"
#include "thread.h"
#include "container.h"
#include "messaging.h"
#include "slab_allocator.h"

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

class async_io;
class callback_allocator;

struct thread_callback_s
{
	struct io_callback_s cb;
	async_io *aio;
	callback *aio_callback;
	callback_allocator *cb_allocator;
	io_request req;
};

/**
 * This slab allocator makes sure all requests in the callback structure
 * are extended requests.
 */
class callback_allocator: public obj_allocator<thread_callback_s>
{
	class callback_initiator: public obj_initiator<thread_callback_s>
	{
	public:
		void init(thread_callback_s *cb) {
			cb->req.init();
		}
	};
public:
	callback_allocator(long increase_size,
			long max_size = MAX_SIZE): obj_allocator<thread_callback_s>(
				increase_size, max_size, new callback_initiator()) {
	}

	virtual int alloc_objs(thread_callback_s **cbs, int num) {
		int ret = obj_allocator<thread_callback_s>::alloc_objs(cbs, num);
		// Make sure all requests are extended requests.
		for (int i = 0; i < ret; i++) {
			if (!cbs[i]->req.is_extended_req()) {
				io_request tmp(true);
				cbs[i]->req = tmp;
			}
		}
		return ret;
	}

	virtual thread_callback_s *alloc_obj() {
		thread_callback_s *cb = obj_allocator<thread_callback_s>::alloc_obj();
		if (!cb->req.is_extended_req()) {
			io_request tmp(true);
			cb->req = tmp;
		}
		return cb;
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
	callback_allocator cb_allocator;
	std::tr1::unordered_map<int, aio_complete_sender *> complete_senders;
	std::tr1::unordered_map<int, fifo_queue<thread_callback_s *> *> remote_tcbs;

	int num_iowait;
	int num_completed_reqs;
	int num_local_alloc;

	struct iocb *construct_req(io_request &io_req, callback_t cb_func);
public:
	/**
	 * @aio_depth_per_file
	 * @node_id: the NUMA node where the disks to be read are connected to.
	 */
	async_io(const logical_file_partition &partition,
			const std::tr1::unordered_map<int, aio_complete_thread *> &complete_threads,
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
