#include "aio_private.h"

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
		void *cb, long res, long res2) {
	thread_callback_s *tcb = (thread_callback_s *) cb;

	tcb->aio->return_cb(tcb);
}

async_io::async_io(const char *names[], int num,
		long size, int aio_depth_per_file, int node_id): buffered_io(names,
			num, size, node_id, O_DIRECT | O_RDWR), AIO_DEPTH(aio_depth_per_file * num)
{
	printf("aio is used\n");
	buf_idx = 0;
	ctx = create_aio_ctx(AIO_DEPTH);
	for (int i = 0; i < AIO_DEPTH * 5; i++) {
		cbs.push_back(new thread_callback_s());
	}
	cb = NULL;
}

void async_io::cleanup()
{
	int slot = max_io_slot(ctx);

	while (slot < AIO_DEPTH) {
		io_wait(ctx, NULL);
		slot = max_io_slot(ctx);
	}
	buffered_io::cleanup();
}

async_io::~async_io()
{
	cleanup();
}

struct iocb *async_io::construct_req(io_request &io_req, callback_t cb_func)
{
	if (cbs.empty()) {
		fprintf(stderr, "no callback object left\n");
		return NULL;
	}

	thread_callback_s *tcb = cbs.front();
	io_callback_s *cb = (io_callback_s *) tcb;
	cbs.pop_front();

	cb->func = cb_func;
	tcb->req = io_req;
	tcb->aio = this;

	assert(tcb->req.get_size() >= MIN_BLOCK_SIZE);
	assert(tcb->req.get_size() % MIN_BLOCK_SIZE == 0);
	assert(tcb->req.get_offset() % MIN_BLOCK_SIZE == 0);
	assert((long) tcb->req.get_buf() % MIN_BLOCK_SIZE == 0);
	int io_type = tcb->req.get_access_method() == READ ? A_READ : A_WRITE;
	if (tcb->req.get_num_bufs() == 1)
		return make_io_request(ctx, get_fd(tcb->req.get_offset()),
				tcb->req.get_size(), tcb->req.get_offset(), tcb->req.get_buf(),
				io_type, cb);
	else {
		int num_bufs = tcb->req.get_num_bufs();
		for (int i = 0; i < num_bufs; i++) {
			assert((long) tcb->req.get_buf(i) % MIN_BLOCK_SIZE == 0);
			assert(tcb->req.get_buf_size(i) % MIN_BLOCK_SIZE == 0);
		}
		return make_io_request(ctx, get_fd(tcb->req.get_offset()),
				/* 
				 * iocb only contains a pointer to the io vector.
				 * the space for the IO vector is stored
				 * in the callback structure.
				 */
				tcb->req.get_vec(), num_bufs, tcb->req.get_offset(),
				io_type, cb);
	}
}

ssize_t async_io::access(io_request *requests, int num)
{
	ssize_t ret = 0;

	while (num > 0) {
		int slot = max_io_slot(ctx);
		if (slot == 0) {
			/*
			 * To achieve the best performance, we need to submit requests
			 * as long as there is a slot available.
			 */
			io_wait(ctx, NULL, 1);
			slot = max_io_slot(ctx);
		}
		struct iocb *reqs[slot];
		int min = slot > num ? num : slot;
		for (int i = 0; i < min; i++) {
			ret += requests->get_size();
			reqs[i] = construct_req(*requests, aio_callback);
			requests++;
		}
		submit_io_request(ctx, reqs, min);
		num -= min;
	}
	return ret;
}
