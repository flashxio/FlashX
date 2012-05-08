#ifndef __AIO_PRIVATE_H__
#define __AIO_PRIVATE_H__

#include <deque>

#include "wpaio.h"
#include "read_private.h"

#ifdef ENABLE_AIO

class aio_private;
void aio_callback(io_context_t, struct iocb*,
		struct io_callback_s *, long, long);

struct thread_callback_s
{
	struct io_callback_s cb;
	aio_private *thread;
	/*
	 * the thread that initiates the request.
	 * it is needed when the request is sent with AIO.
	 */
	thread_private *initiator;
	void *priv;
};

class aio_private: public read_private
{
	int buf_idx;
	struct aio_ctx *ctx;
	std::deque<thread_callback_s *> cbs;

public:
	/**
	 * @names: the names of files to be accessed
	 * @num: the number of files
	 * @size: the size of data to be accessed in all files
	 * @idx: the thread index
	 * @entry_size: the size of an entry to be accessed.
	 */
	aio_private(const char *names[], int num, long size, int idx,
			int entry_size);

	~aio_private();

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method);
	ssize_t access(io_request *requests, int num, int access_method);
	struct iocb *construct_req(char *buf, off_t offset, ssize_t size,
			int access_method, callback_t callback,
			thread_private *initiator, void *priv);

	bool support_bulk() {
		return true;
	}

	virtual void cleanup();

	void return_cb(thread_callback_s *tcb) {
		cbs.push_back(tcb);
		io_callback_s *cb = (io_callback_s *) tcb;
		read_bytes += cb->size;

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
