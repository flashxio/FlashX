#ifndef __AIO_PRIVATE_H__
#define __AIO_PRIVATE_H__

#include <deque>

#include "wpaio.h"
#include "read_private.h"

#ifdef ENABLE_AIO

void aio_callback(io_context_t, struct iocb*,
		struct io_callback_s *, long, long);

class async_io;
struct thread_callback_s
{
	struct io_callback_s cb;
	async_io *aio;
	/*
	 * the thread that initiates the request.
	 * it is needed when the request is sent with AIO.
	 */
	io_interface *initiator;
	void *priv;
};

class async_io: public buffered_io
{
	int buf_idx;
	struct aio_ctx *ctx;
	std::deque<thread_callback_s *> cbs;
	callback *cb;

	struct iocb *construct_req(char *buf, off_t offset, ssize_t size,
			int access_method, callback_t callback,
			io_interface *initiator, void *priv);
public:
	/**
	 * @names: the names of files to be accessed
	 * @num: the number of files
	 * @size: the size of data to be accessed in all files
	 */
	async_io(const char *names[], int num, long size);

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

	void return_cb(thread_callback_s *tcb) {
		cbs.push_back(tcb);
		io_callback_s *cb = (io_callback_s *) tcb;

		// TODO it only support READ right now.
		if (this->cb) {
			io_request req(cb->buf, cb->offset, cb->size, READ,
					tcb->initiator, tcb->priv);
			this->cb->invoke(&req);
		}
	}
};
#endif

#endif
