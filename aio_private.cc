#include "aio_private.h"

#define AIO_DEPTH (40 * num_open_files() / nthreads)

#define EVEN_DISTRIBUTE

const int MAX_BUF_REQS = 1024 * 3;

/* 
 * each file gets the same number of outstanding requests.
 */
#ifdef EVEN_DISTRIBUTE
#define MAX_OUTSTANDING_NREQS (AIO_DEPTH / num_open_files())
#define ALLOW_DROP
#endif

void aio_callback(io_context_t ctx, struct iocb* iocb,
		struct io_callback_s *cb, long res, long res2) {
	thread_callback_s *tcb = (thread_callback_s *) cb;

	tcb->thread->return_cb(tcb);
}

aio_private::aio_private(const char *names[], int num, long size,
		int idx, int entry_size): read_private(names, num, size, idx,
			entry_size, O_DIRECT | O_RDONLY)
{
	printf("aio is used\n");
	buf_idx = 0;
	ctx = create_aio_ctx(AIO_DEPTH);
	for (int i = 0; i < AIO_DEPTH * 5; i++) {
		cbs.push_back(new thread_callback_s());
	}
}

void aio_private::cleanup()
{
	int slot = max_io_slot(ctx);

	while (slot < AIO_DEPTH) {
		io_wait(ctx, NULL);
		slot = max_io_slot(ctx);
	}
	read_private::cleanup();
}

aio_private::~aio_private()
{
	cleanup();
}

struct iocb *aio_private::construct_req(char *buf, off_t offset,
		ssize_t size, int access_method, callback_t cb_func,
		thread_private *initiator, void *priv)
{
	struct iocb *req = NULL;

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
	tcb->priv = priv;

	assert(size >= MIN_BLOCK_SIZE);
	assert(size % MIN_BLOCK_SIZE == 0);
	assert(offset % MIN_BLOCK_SIZE == 0);
	assert((long) buf % MIN_BLOCK_SIZE == 0);
	req = make_io_request(ctx, get_fd(offset), size, offset, buf, A_READ, cb);
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

	/* the initiator of the request must be itself. */
	req = construct_req(buf, offset, size, access_method, aio_callback,
			this, NULL);
	submit_io_request(ctx, &req, 1);
	return 0;
}

ssize_t aio_private::access(io_request *requests, int num, int access_method)
{
	ssize_t ret = 0;

	while (num > 0) {
		int slot = max_io_slot(ctx);
		if (slot == 0) {
			io_wait(ctx, NULL, 10);
			slot = max_io_slot(ctx);
		}
		struct iocb *reqs[slot];
		int min = slot > num ? num : slot;
		for (int i = 0; i < min; i++) {
			reqs[i] = construct_req(requests->get_buf(),
					requests->get_offset(), requests->get_size(),
					requests->get_access_method(), aio_callback, 
					/*
					 * The initiator thread is normally itself.
					 * However, if the request is received from 
					 * another thread with message passing, 
					 * it will be another thread.
					 */
					requests->get_thread(), requests->get_priv());
			ret += requests->get_size();
			requests++;
		}
		submit_io_request(ctx, reqs, min);
		num -= min;
	}
	return ret;
}
