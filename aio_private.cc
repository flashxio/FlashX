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
#else
#define MAX_OUTSTANDING_NREQS (AIO_DEPTH)
#endif

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
	buf_idx = 0;
	ctx = create_aio_ctx(AIO_DEPTH);
	for (int i = 0; i < AIO_DEPTH * 5; i++) {
		cbs.push_back(new thread_callback_s());
	}
	reqs_array = new std::deque<io_request>[this->num_open_files()];
	outstanding_nreqs = new int[this->num_open_files()];
	memset(outstanding_nreqs, 0, sizeof(int) * this->num_open_files());
}

aio_private::~aio_private()
{
	int slot = max_io_slot(ctx);

	while (slot < AIO_DEPTH) {
		io_wait(ctx, NULL);
		slot = max_io_slot(ctx);
	}
	delete [] reqs_array;
	delete [] outstanding_nreqs;
}

struct iocb *aio_private::construct_req(char *buf, off_t offset,
		ssize_t size, int access_method, callback_t cb_func)
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

	/* for simplicity, I assume all request sizes are smaller than a page size */
	assert(size <= PAGE_SIZE);
	if (ROUND_PAGE(offset) == offset
			&& (long) buf == ROUND_PAGE(buf)
			&& size == PAGE_SIZE) {
		req = make_io_request(ctx, get_fd(offset), PAGE_SIZE, offset,
				buf, A_READ, cb);
	}
	else {
		assert(size == PAGE_SIZE);
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
			ret += requests->get_size();
			requests++;
		}
		submit_io_request(ctx, reqs, min);
		num -= min;
	}
	return ret;
}

const int MAX_REQ_SIZE = 128;

/**
 * some SSDs can process requests faster than others. 
 * The idea here is to send more requests to faster SSDs
 * and send fewer requests to slower SSDs.
 * To do so, I need to track the number of outstanding requests
 * to each SSD.
 */
ssize_t aio_private::access(io_request *requests, int num, int access_method)
{
	if (num_open_files() == 1)
		return process_reqs(requests, num);

	/* first put all requests in the queues. */
	buffer_reqs(requests, num);

	ssize_t ret = 0;
	io_request req_buf[MAX_REQ_SIZE];
	int req_buf_size = 0;
	int remaining = 0;	// the number of remaining requests
	int busy = 0;

	while(true) {
		busy = 0;
		remaining = 0;
		for (int i = 0; i < num_open_files(); i++) {
			int available = MAX_OUTSTANDING_NREQS - outstanding_nreqs[i];
			for (int j = 0; j < available && !reqs_array[i].empty(); j++) {
				req_buf[req_buf_size++] = reqs_array[i].front();
				reqs_array[i].pop_front();
				outstanding_nreqs[i]++;
				/*
				 * if the temporary request buffer is full,
				 * process all requests.
				 */
				if (req_buf_size == MAX_REQ_SIZE) {
					ret += process_reqs(req_buf, req_buf_size);
					req_buf_size = 0;
				}
			}
			/*
			 * if there are more requests than we can send to a file,
			 * we consider it as busy.
			 */
			if (outstanding_nreqs[i] == MAX_OUTSTANDING_NREQS
					&& !reqs_array[i].empty()) {
				busy++;
			}
			remaining += reqs_array[i].size();
		}

		/* process the remaining requests in the buffer. */
		if (req_buf_size > 0) {
			ret += process_reqs(req_buf, req_buf_size);
			req_buf_size = 0;
		}

		/*
		 * if all files are busy or there are too many requests
		 * buffered, we should wait and then process the remaining
		 * requests again.
		 */
		if (busy == num_open_files()
#ifndef ALLOW_DROP
				|| remaining > MAX_BUF_REQS
#endif
				)
			io_wait(ctx, NULL, 10);
		else
			break;
	}

#ifdef ALLOW_DROP
	if (remaining > MAX_BUF_REQS) {
		unsigned int max = reqs_array[0].size();
		int max_idx = 0;
		for (int i = 1; i < num_open_files(); i++)
			if (max < reqs_array[i].size()) {
				max = reqs_array[i].size();
				max_idx = i;
			}
		/* 
		 * there are too many requests for the file accumulated.
		 * it may be caused by the slow underlying block device.
		 * drop half of these requests, so we can send more
		 * requests to other files.
		 */
		int num_drops = max / 2;
		for (int i = 0; i < num_drops; i++) {
			io_request req = reqs_array[max_idx].front();
			reqs_array[max_idx].pop_front();
			drop_req(&req);
		}
//		printf("drop %d requests to file %d, max: %d, remaining: %d\n",
//				num_drops, max_idx, max, remaining);
	}
#endif
	
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
