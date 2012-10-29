#include "balanced_aio_private.h"

const int MAX_REQ_SIZE = 128;

void aio_private::buffer_reqs(io_request *requests, int num)
{
	for (int i = 0; i < num; i++) {
		int fd_idx = get_fd_idx(requests->get_offset());
		reqs_array[fd_idx].push_back(*requests);
		requests++;
	}
}

/**
 * some SSDs can process requests faster than others. 
 * The idea here is to send more requests to faster SSDs
 * and send fewer requests to slower SSDs.
 * To do so, I need to track the number of outstanding requests
 * to each SSD.
 */
ssize_t balanced_aio_private::access(io_request *requests,
		int num, int access_method)
{
	if (num_open_files() == 1)
		return aio_private::access(requests, num, access_method);

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
					ret += aio_private::access(req_buf,
							req_buf_size, access_method);
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
			ret += aio_private::access(req_buf, req_buf_size, access_method);
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
void balanced_aio_private::cleanup()
{
	for (int i = 0; i < num_open_files(); i++) {
		for (unsigned int j = 0; j < reqs_array[i].size(); j++) {
			io_request req = reqs_array[i][j];
			access(req.get_buf(), req.get_offset(),
					req.get_size(), req.get_access_method());
		}
	}

	aio_private::cleanup();
}
