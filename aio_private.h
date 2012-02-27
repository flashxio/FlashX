#ifndef __AIO_PRIVATE_H__
#define __AIO_PRIVATE_H__

#include "read_private.h"

#ifdef ENABLE_AIO

#define AIO_DEPTH 32

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
			int entry_size): read_private(names, num, size, idx,
			entry_size, O_DIRECT | O_RDWR) {
		printf("aio is used\n");
		pages = (char *) valloc(PAGE_SIZE * 4096);
		buf_idx = 0;
		ctx = create_aio_ctx(AIO_DEPTH);
		for (int i = 0; i < AIO_DEPTH * 5; i++) {
			cbs.push_back(new thread_callback_s());
		}
	}

	~aio_private() {
		int slot = max_io_slot(ctx);

		while (slot < AIO_DEPTH) {
			io_wait(ctx, NULL);
			slot = max_io_slot(ctx);
		}
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method) {
		struct iocb *req;
		int slot = max_io_slot(ctx);

		assert(access_method == READ);
		if (slot == 0) {
			io_wait(ctx, NULL);
		}

		if (cbs.empty()) {
			fprintf(stderr, "no callback object left\n");
			return -1;
		}

		thread_callback_s *tcb = cbs.front();
		io_callback_s *cb = (io_callback_s *) tcb;
		cbs.pop_front();
		cb->buf = buf;
		cb->offset = offset;
		cb->size = size;
		cb->func = aio_callback;
		tcb->thread = this;

		/* for simplicity, I assume all request sizes are smaller than a page size */
		assert(size <= PAGE_SIZE);
		if (ROUND_PAGE(offset) == offset
				&& (long) buf == ROUND_PAGE(buf)
				&& size == PAGE_SIZE) {
			req = make_io_request(ctx, get_fd(), PAGE_SIZE, offset,
					buf, A_READ, cb);
		}
		else {
			buf_idx++;
			if (buf_idx == 4096)
				buf_idx = 0;
			char *page = pages + buf_idx * PAGE_SIZE;
			req = make_io_request(ctx, get_fd(), PAGE_SIZE, ROUND_PAGE(offset),
					page, A_READ, cb);
		}
		submit_io_request(ctx, &req, 1);
		return 0;
	}

	void return_cb(thread_callback_s *cb) {
		cbs.push_back(cb);
	}
};

void aio_callback(io_context_t ctx, struct iocb* iocb,
		struct io_callback_s *cb, long res, long res2) {
	thread_callback_s *tcb = (thread_callback_s *) cb;

	memcpy(cb->buf, ((char *) iocb->u.c.buf)
			+ (cb->offset - ROUND_PAGE(cb->offset)), cb->size);
	if(*(unsigned long *) cb->buf != cb->offset / sizeof(long))
		printf("%ld %ld\n", *(unsigned long *) cb->buf, cb->offset / sizeof(long));
	assert(*(unsigned long *) cb->buf == cb->offset / sizeof(long));
	tcb->thread->return_cb(tcb);
	tcb->thread->read_bytes += cb->size;
}
#endif

#endif
