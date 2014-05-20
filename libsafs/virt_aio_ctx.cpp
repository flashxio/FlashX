/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <algorithm>

#include "virt_aio_ctx.h"

const int MIN_READ_DELAY = 100;
const int MAX_RAND_READ_DELAY = 100;		// in microseconds
const int MIN_WRITE_DELAY = 200;
const int MAX_RAND_WRITE_DELAY = 500;		// in microseconds

struct timeval add2timeval(const struct timeval time, long added_delay)
{
	struct timeval res = time;
	res.tv_usec += added_delay;
	res.tv_sec += res.tv_usec / 1000000;
	res.tv_usec = res.tv_usec % 1000000;
	return res;
}

class naive_ssd_perf_model: public ssd_perf_model
{
public:
	virtual long get_read_delay(off_t off, size_t size);
	virtual long get_write_delay(off_t off, size_t size);
};

long naive_ssd_perf_model::get_read_delay(off_t off, size_t size)
{
	// We introduce some random delay in each request.
	// The delay is in microseconds.
	int rand_delay = random() % MAX_RAND_READ_DELAY;
	int num_pages = size / PAGE_SIZE;
	if (num_pages == 0)
		num_pages = 1;
	return rand_delay + MIN_READ_DELAY * num_pages;
}

long naive_ssd_perf_model::get_write_delay(off_t off, size_t size)
{
	int rand_delay = random() % MAX_RAND_WRITE_DELAY;
	int num_pages = size / PAGE_SIZE;
	if (num_pages == 0)
		num_pages = 1;
	return rand_delay + MIN_WRITE_DELAY * num_pages;
}

virt_aio_ctx::virt_aio_ctx(virt_data *data, int node_id,
		int max_aio): aio_ctx(node_id, max_aio), pending_reqs(node_id, max_aio)
{
	this->max_aio = max_aio;
	this->data = data;
	this->model = new naive_ssd_perf_model();

	read_bytes = 0;
	write_bytes = 0;
	read_bytes_ps = 0;
	write_bytes_ps = 0;
	memset(&prev_print_time, 0, sizeof(prev_print_time));
}

struct comp_issued_request
{
	bool operator() (const struct req_entry &req1,
			const struct req_entry &req2) {
		return time_diff_us(req1.issue_time, req2.issue_time) > 0;
	}
} issued_req_comparator;

size_t get_size(struct iocb *req)
{
	if (req->aio_lio_opcode == IO_CMD_PREAD
			|| req->aio_lio_opcode == IO_CMD_PWRITE) {
		return req->u.c.nbytes;
	}
	else {
		size_t size = 0;
		int num_vecs = req->u.c.nbytes;
		struct iovec *iov = (struct iovec *) req->u.c.buf;
		for (int i = 0; i < num_vecs; i++) {
			size += iov[i].iov_len;
		}
		return size;
	}
}

void virt_aio_ctx::submit_io_request(struct iocb* ioq[], int num)
{
	struct req_entry entries[num];
	struct timeval curr;

	gettimeofday(&curr, NULL);
	for (int i = 0; i < num; i++) {
		entries[i].req = ioq[i];

		off_t off = ioq[i]->u.c.offset;
		if (ioq[i]->aio_lio_opcode == IO_CMD_PREAD
				|| ioq[i]->aio_lio_opcode == IO_CMD_PREADV) {
			long delay = model->get_read_delay(off, get_size(ioq[i]));
			entries[i].issue_time = add2timeval(curr, delay);
		}
		else {
			long delay = model->get_write_delay(off, get_size(ioq[i]));
			entries[i].issue_time = add2timeval(curr, delay);
		}
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

	long time_diff = time_diff_us(prev_print_time, curr);
	if (time_diff >= params.get_vaio_print_freq()) {
		printf("read %.2fMB/s, write %.2fMB/s\n",
				((double) read_bytes_ps) / 1024 / 1024 / (((double) time_diff) / 1000000),
				((double) write_bytes_ps) / 1024 / 1024 / (((double) time_diff) / 1000000));
		read_bytes_ps = 0;
		write_bytes_ps = 0;
		prev_print_time = curr;
	}

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
		if (iocbs[i]->aio_lio_opcode == IO_CMD_PREADV) {
			off_t offset = iocbs[i]->u.c.offset;
			int num_vecs = iocbs[i]->u.c.nbytes;
			struct iovec *iov = (struct iovec *) iocbs[i]->u.c.buf;
			int fd = iocbs[i]->aio_fildes;
			for (int j = 0; j < num_vecs; j++) {
				if (params.is_verify_content()) {
					data->create_data(fd, (char *) iov[j].iov_base,
							iov[j].iov_len, offset);
					offset += iov[j].iov_len;
				}
				read_bytes_ps += iov[j].iov_len;
				read_bytes += iov[j].iov_len;
			}
		}
		else if (iocbs[i]->aio_lio_opcode == IO_CMD_PREAD) {
			if (params.is_verify_content()) {
				data->create_data(iocbs[i]->aio_fildes,
						(char *) iocbs[i]->u.c.buf, iocbs[i]->u.c.nbytes,
						iocbs[i]->u.c.offset);
			}
			read_bytes_ps += iocbs[i]->u.c.nbytes;
			read_bytes += iocbs[i]->u.c.nbytes;
		}
		else if (iocbs[i]->aio_lio_opcode == IO_CMD_PWRITEV) {
			off_t offset = iocbs[i]->u.c.offset;
			int num_vecs = iocbs[i]->u.c.nbytes;
			struct iovec *iov = (struct iovec *) iocbs[i]->u.c.buf;
			int fd = iocbs[i]->aio_fildes;
			for (int j = 0; j < num_vecs; j++) {
				if (params.is_verify_content()) {
					assert(data->verify_data(fd, (char *) iov[j].iov_base,
								iov[j].iov_len, offset));
					offset += iov[j].iov_len;
				}
				write_bytes_ps += iov[j].iov_len;
				write_bytes += iov[j].iov_len;
			}
		}
		else if (iocbs[i]->aio_lio_opcode == IO_CMD_PWRITE) {
			if (params.is_verify_content()) {
				assert(data->verify_data(iocbs[i]->aio_fildes,
							(char *) iocbs[i]->u.c.buf, iocbs[i]->u.c.nbytes,
							iocbs[i]->u.c.offset));
			}
			write_bytes_ps += iocbs[i]->u.c.nbytes;
			write_bytes += iocbs[i]->u.c.nbytes;
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
