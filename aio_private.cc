#include "aio_private.h"

#define AIO_DEPTH 128

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
	reqs_array = new std::deque<io_request>[this->num_open_files()];
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

void aio_private::buffer_reqs(io_request *requests, int num)
{
	for (int i = 0; i < num; i++) {
		int fd_idx = get_fd_idx(requests->get_offset());
		reqs_array[fd_idx].push_back(*requests);
		requests++;
	}
}

ssize_t aio_private::process_reqs(io_request *requests, int num)
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
					requests->get_access_method(), aio_callback1);
			requests++;
			ret += requests->get_size();
		}
		submit_io_request(ctx, reqs, min);
		num -= min;
		min = min / 2;
		if (min == 0)
			min = 1;
	}
	return ret;
}

const int MAX_REQ_SIZE = 128;

/**
 * Random workload may cause some load imbalance at some particular moments
 * when the requests are distributed to multiple files.
 * If there is only one file to read, simply submit requests with AIO.
 * Otherwise, before submitting requests to AIO, I need to rebalance them first.
 * The way to rebalance requests is to have a queue for each file, and buffer
 * all requests in queues according to the files where they are submitted to.
 * Fetch requests from queues in a round-robin fashion, so it's guaranteed that
 * all requests will be distributed to files evenly.
 *
 * If num is 0, the invocation processes up to num_open_files() reuqests.
 * This is good, so the number of requests in each queue won't be extremely
 * large if the distribution is highly skewed.
 */
ssize_t aio_private::access(io_request *requests, int num, int access_method)
{
	if (num_open_files() == 1)
		return process_reqs(requests, num);

	ssize_t ret = 0;
	int one_empty = false;
	buffer_reqs(requests, num);
	io_request req_buf[MAX_REQ_SIZE];
	int req_buf_size = 0;

	/* if a request queue for a file is empty, stop processing. */
	while (!one_empty) {
		for (int i = 0; i < num_open_files(); i++) {
			if (reqs_array[i].empty()) {
				one_empty = true;
				continue;
			}
			req_buf[req_buf_size++] = reqs_array[i].front();
			reqs_array[i].pop_front();
			/*
			 * if the temporary request buffer is full,
			 * process all requests.
			 */
			if (req_buf_size == MAX_REQ_SIZE) {
				ret += process_reqs(req_buf, req_buf_size);
				req_buf_size = 0;
			}
		}
	}
	/* process the remaining requests in the buffer. */
	if (req_buf_size > 0)
		ret += process_reqs(req_buf, req_buf_size);
	return ret;
}

/* process remaining requests in the queues. */
void aio_private::cleanup()
{
	read_private::cleanup();

	for (int i = 0; i < num_open_files(); i++) {
		for (unsigned int j = 0; j < reqs_array[i].size(); j++) {
			io_request req = reqs_array[i][j];
			access(req.get_buf(), req.get_offset(),
					req.get_size(), req.get_access_method());
		}
	}
}
