#include <algorithm>

#include "virt_aio_ctx.h"

const int MIN_DELAY = 100;
const int MAX_RAND_DELAY = 100;		// in microseconds

struct timeval add2timeval(const struct timeval time, long added_delay)
{
	struct timeval res = time;
	res.tv_usec += added_delay;
	res.tv_sec += res.tv_usec / 1000000;
	res.tv_usec = res.tv_usec % 1000000;
	return res;
}

struct comp_issued_request
{
	bool operator() (const struct req_entry &req1,
			const struct req_entry &req2) {
		return time_diff_us(req1.issue_time, req2.issue_time) > 0;
	}
} issued_req_comparator;

void virt_aio_ctx::submit_io_request(struct iocb* ioq[], int num)
{
	struct req_entry entries[num];
	struct timeval curr;

	gettimeofday(&curr, NULL);
	for (int i = 0; i < num; i++) {
		entries[i].req = ioq[i];

		// We introduce some random delay in each request.
		// The delay is in microseconds.
		int rand_delay = random() % MAX_RAND_DELAY;
		entries[i].issue_time = add2timeval(curr, rand_delay + MIN_DELAY);
	}
	std::sort(entries, entries + num, issued_req_comparator);
	assert(time_diff_us(entries[0].issue_time,
				entries[num - 1].issue_time) >= 0);

	int num_existing = pending_reqs.get_num_entries();
	struct req_entry origs[num_existing];
	int ret = pending_reqs.fetch(origs, num_existing);
	assert(ret == num_existing);

	int tot = num_existing + num;
	struct req_entry merge_buf[tot];
	assert(tot <= max_aio);
	std::merge(entries, entries + num, origs, origs + num_existing, merge_buf,
			issued_req_comparator);
	assert(time_diff_us(merge_buf[0].issue_time,
				merge_buf[tot - 1].issue_time) >= 0);
	pending_reqs.add(merge_buf, tot);
}

int virt_aio_ctx::io_wait(struct timespec* to, int num)
{
	struct req_entry entries[num];
	int ret = pending_reqs.fetch(entries, num);

	// If there aren't any requests we can wait for, return immediately.
	if (ret == 0)
		return 0;

	// We have to wait until the specified number of requests are completed.
	struct timeval curr;
	gettimeofday(&curr, NULL);
	do {
		long sleep_time = time_diff_us(curr, entries[ret - 1].issue_time);
		if (sleep_time <= 0)
			break;
		struct timespec req = {0, sleep_time * 1000};
		struct timespec rem;
		nanosleep(&req, &rem);
		gettimeofday(&curr, NULL);
	} while (time_diff_us(curr, entries[ret - 1].issue_time) > 0);

	// Notify the application of the completion of the requests.
	// TODO I should fill the requests with right data.
	struct iocb *iocbs[ret];
	long res[ret];
	long res2[ret];
	io_callback_s *cbs[ret];
	callback_t cb_func = NULL;
	for (int i = 0; i < ret; i++)
	{
		cbs[i] = (io_callback_s *) entries[i].req->data;
		if (cb_func == NULL)
			cb_func = cbs[i]->func;
		assert(cb_func == cbs[i]->func);
		iocbs[i] = entries[i].req;
		if (params.is_verify_content()) {
			if (iocbs[i]->aio_lio_opcode == IO_CMD_PREADV) {
				off_t offset = iocbs[i]->u.c.offset;
				int num_vecs = iocbs[i]->u.c.nbytes;
				struct iovec *iov = (struct iovec *) iocbs[i]->u.c.buf;
				int fd = iocbs[i]->aio_fildes;
				for (int j = 0; j < num_vecs; j++) {
					data->create_data(fd, (char *) iov[j].iov_base,
							iov[j].iov_len, offset);
					offset += iov[j].iov_len;
				}
			}
			else if (iocbs[i]->aio_lio_opcode == IO_CMD_PREAD) {
				data->create_data(iocbs[i]->aio_fildes,
						(char *) iocbs[i]->u.c.buf, iocbs[i]->u.c.nbytes,
						iocbs[i]->u.c.offset);
			}
			else if (iocbs[i]->aio_lio_opcode == IO_CMD_PWRITEV) {
				off_t offset = iocbs[i]->u.c.offset;
				int num_vecs = iocbs[i]->u.c.nbytes;
				struct iovec *iov = (struct iovec *) iocbs[i]->u.c.buf;
				int fd = iocbs[i]->aio_fildes;
				for (int j = 0; j < num_vecs; j++) {
					assert(data->verify_data(fd, (char *) iov[j].iov_base,
							iov[j].iov_len, offset));
					offset += iov[j].iov_len;
				}
			}
			else if (iocbs[i]->aio_lio_opcode == IO_CMD_PWRITE) {
				assert(data->verify_data(iocbs[i]->aio_fildes,
						(char *) iocbs[i]->u.c.buf, iocbs[i]->u.c.nbytes,
						iocbs[i]->u.c.offset));
			}
		}
		res[i] = 0;
		res2[i] = 0;
	}

	cb_func(NULL, iocbs, (void **) cbs, res, res2, ret);
	destroy_io_requests(iocbs, ret);
	return ret;
}

int virt_aio_ctx::max_io_slot()
{
	return max_aio - pending_reqs.get_num_entries();
}
