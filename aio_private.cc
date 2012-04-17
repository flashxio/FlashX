#include "aio_private.h"

void aio_callback(io_context_t ctx, struct iocb* iocb,
		struct io_callback_s *cb, long res, long res2) {
	thread_callback_s *tcb = (thread_callback_s *) cb;

	memcpy(cb->buf, ((char *) iocb->u.c.buf)
			+ (cb->offset - ROUND_PAGE(cb->offset)), cb->size);
//	if(*(unsigned long *) cb->buf != cb->offset / sizeof(long))
//		printf("%ld %ld\n", *(unsigned long *) cb->buf, cb->offset / sizeof(long));
//	assert(*(unsigned long *) cb->buf == cb->offset / sizeof(long));
	tcb->thread->return_cb(tcb);
	tcb->thread->read_bytes += cb->size;
}

void aio_callback1(io_context_t ctx, struct iocb* iocb,
		struct io_callback_s *cb, long res, long res2) {
	thread_callback_s *tcb = (thread_callback_s *) cb;

	memcpy(cb->buf, ((char *) iocb->u.c.buf)
			+ (cb->offset - ROUND_PAGE(cb->offset)), cb->size);
//	if(*(unsigned long *) cb->buf != cb->offset / sizeof(long))
//		printf("%ld %ld\n", *(unsigned long *) cb->buf, cb->offset / sizeof(long));
//	assert(*(unsigned long *) cb->buf == cb->offset / sizeof(long));
	tcb->thread->return_cb1(tcb);
	tcb->thread->read_bytes += cb->size;
}

aio_private::aio_private(const char *names[], int num, long size,
		int idx, int entry_size): read_private(names, num, size, idx,
			entry_size, O_DIRECT | O_RDWR)
{
	printf("aio is used\n");
	pages = (char *) valloc(PAGE_SIZE * 4096);
	buf_idx = 0;
	ctx = create_aio_ctx(AIO_DEPTH);
	for (int i = 0; i < AIO_DEPTH * 5; i++) {
		cbs.push_back(new thread_callback_s());
	}
}

aio_private::~aio_private()
{
	int slot = max_io_slot(ctx);

	while (slot < AIO_DEPTH) {
		io_wait(ctx, NULL);
		slot = max_io_slot(ctx);
	}
}

struct iocb *aio_private::construct_req(char *buf, off_t offset,
		ssize_t size, int access_method, callback_t cb_func)
{
	struct iocb *req;

	if (cbs.empty()) {
		fprintf(stderr, "no callback object left\n");
		return NULL;
	}

	thread_callback_s *tcb = cbs.front();
	io_callback_s *cb = (io_callback_s *) tcb;
	cbs.pop_front();
	cb->buf = buf;
	cb->offset = offset;
	cb->size = size;
	cb->func = cb_func;
	tcb->thread = this;

	/* for simplicity, I assume all request sizes are smaller than a page size */
	assert(size <= PAGE_SIZE);
	if (ROUND_PAGE(offset) == offset
			&& (long) buf == ROUND_PAGE(buf)
			&& size == PAGE_SIZE) {
		req = make_io_request(ctx, get_fd(offset), PAGE_SIZE, offset,
				buf, A_READ, cb);
	}
	else {
		buf_idx++;
		if (buf_idx == 4096)
			buf_idx = 0;
		char *page = pages + buf_idx * PAGE_SIZE;
		req = make_io_request(ctx, get_fd(offset), PAGE_SIZE, ROUND_PAGE(offset),
				page, A_READ, cb);
	}
	return req;
}

ssize_t aio_private::access(char *buf, off_t offset,
		ssize_t size, int access_method) {
	struct iocb *req;
	int slot = max_io_slot(ctx);

	assert(access_method == READ);
	if (slot == 0) {
		io_wait(ctx, NULL);
	}

	req = construct_req(buf, offset, size, access_method, aio_callback);
	submit_io_request(ctx, &req, 1);
	return 0;
}

ssize_t aio_private::access(io_request *requests, int num, int access_method)
{
	ssize_t ret = 0;
//	printf("send %d requests\n", num);
	while (num > 0) {
		int slot = max_io_slot(ctx);
		struct iocb *reqs[slot];
		int min = slot > num ? num : slot;
		for (int i = 0; i < min; i++) {
			reqs[i] = construct_req(requests->get_buf(),
					requests->get_offset(), requests->get_size(),
					access_method, aio_callback1);
			requests++;
			ret += requests->get_size();
		}
//		printf("submit %d requests\n", min);
		submit_io_request(ctx, reqs, min);
		num -= min;
		min = min / 2;
		if (min == 0)
			min = 1;
		int done = io_wait(ctx, NULL, min);
//		printf("finish %d requests\n", done);
	}
//	printf("read %ld bytes\n", ret);
	return ret;
}
