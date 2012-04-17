#ifndef __AIO_PRIVATE_H__
#define __AIO_PRIVATE_H__

#include <deque>

#include "wpaio.h"
#include "read_private.h"

#ifdef ENABLE_AIO

#define AIO_DEPTH 128

class aio_private;
void aio_callback(io_context_t, struct iocb*,
		struct io_callback_s *, long, long);

struct thread_callback_s
{
	struct io_callback_s cb;
	aio_private *thread;
};

class aio_private: public read_private
{
	char *pages;
	int buf_idx;
	struct aio_ctx *ctx;
	std::deque<thread_callback_s *> cbs;

public:
	aio_private(const char *names[], int num, long size, int idx,
			int entry_size);

	~aio_private();

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method);
	ssize_t access(io_request *requests, int num, int access_method);
	struct iocb *construct_req(char *buf, off_t offset,
			ssize_t size, int access_method, callback_t callback);

	bool support_bulk() {
		return false;
//		return true;
	}

	void return_cb(thread_callback_s *cb) {
		cbs.push_back(cb);
	}

	void return_cb1(thread_callback_s *tcb) {
		cbs.push_back(tcb);
		io_callback_s *cb = (io_callback_s *) tcb;
		buf->free_entry(cb->buf);
	}
};
#endif

#endif
