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

#include <stdlib.h>

#include "slab_allocator.h"

#include "direct_comp_access.h"
#include "cache.h"
#include "remote_access.h"
#include "comp_io_scheduler.h"

namespace safs
{

// 8MB.
static const size_t MAX_PEND_COMP_SIZE = 8 * 1024 * 1024;

/*
 * This class represents the data requested by user compute.
 * The data is stored in the memory buffer read from the disks by remote I/O.
 */
class direct_byte_array: public page_byte_array
{
	// The offset of the data requested by a user compute on the disks.
	off_t req_off;
	size_t size;
	// This buffer contains the data requested by a user compute.
	char *buf;

	void assign(direct_byte_array &arr) {
		this->req_off = arr.req_off;
		this->size = arr.size;
		this->buf = arr.buf;
		arr.buf = NULL;
	}

	direct_byte_array(direct_byte_array &arr) {
		assign(arr);
	}

	direct_byte_array &operator=(direct_byte_array &arr) {
		assign(arr);
		return *this;
	}

	// The offset of the data in the memory buffer on the disks.
	// This offset isn't aligned to the page size.
	off_t get_buf_offset() const {
		return ROUND(req_off, MIN_BLOCK_SIZE);
	}
public:
	direct_byte_array(byte_array_allocator &alloc): page_byte_array(alloc) {
		req_off = 0;
		size = 0;
		buf = NULL;
	}

	// When the memory buffer is passed to the byte array, the array has
	// the ownership of the buffer. i.e., the array has the responsibility
	// of deallocating the buffer.
	direct_byte_array(off_t req_off, size_t size, char *buf,
			byte_array_allocator &alloc): page_byte_array(alloc) {
		this->req_off = req_off;
		this->size = size;
		this->buf = buf;
	}

	~direct_byte_array() {
		if (buf)
			free(buf);
	}

	virtual off_t get_offset() const {
		return req_off;
	}

	virtual off_t get_offset_in_first_page() const {
		return req_off - get_buf_offset();
	}

	virtual const char *get_page(int pg_idx) const {
		return buf + pg_idx * PAGE_SIZE;
	}

	virtual size_t get_size() const {
		return size;
	}

	void lock() {
		ABORT_MSG("lock isn't implemented");
	}

	void unlock() {
		ABORT_MSG("unlock isn't implemented");
	}

	page_byte_array *clone() {
		direct_byte_array *arr = (direct_byte_array *) get_allocator().alloc();
		*arr = *this;
		return arr;
	}
};

class direct_byte_array_allocator: public byte_array_allocator
{
	class array_initiator: public obj_initiator<direct_byte_array>
	{
		direct_byte_array_allocator *alloc;
	public:
		array_initiator(direct_byte_array_allocator *alloc) {
			this->alloc = alloc;
		}

		virtual void init(direct_byte_array *obj) {
			new (obj) direct_byte_array(*alloc);
		}
	};

	class array_destructor: public obj_destructor<direct_byte_array>
	{
	public:
		void destroy(direct_byte_array *obj) {
			obj->~direct_byte_array();
		}
	};

	obj_allocator<direct_byte_array> allocator;
public:
	direct_byte_array_allocator(thread *t): allocator(
			"byte-array-allocator", t->get_node_id(), false, 1024 * 1024,
			params.get_max_obj_alloc_size(),
			obj_initiator<direct_byte_array>::ptr(new array_initiator(this)),
			obj_destructor<direct_byte_array>::ptr(new array_destructor())) {
	}

	virtual page_byte_array *alloc() {
		return allocator.alloc_obj();
	}

	virtual void free(page_byte_array *arr) {
		allocator.free((direct_byte_array *) arr);
	}
};

/*
 * When issuing an I/O request to remote I/O, the offset and the size of
 * the request need to be aligned to BLOCK SIZE. We use this data structure
 * to keep the original offset and size of the request.
 */
class orig_comp_request
{
	off_t off;
	size_t size;
	user_compute *compute;
public:
	orig_comp_request(off_t off, size_t size, user_compute *compute) {
		this->off = off;
		this->size = size;
		this->compute = compute;
	}

	off_t get_offset() const {
		return off;
	}

	size_t get_size() const {
		return size;
	}

	user_compute *get_compute() const {
		return compute;
	}
};

class comp_callback: public callback
{
	direct_comp_io &io;
public:
	comp_callback(direct_comp_io &_io): io(_io) {
	}

	virtual int invoke(io_request *reqs[], int num) {
		for (int i = 0; i < num; i++)
			io.complete_req(*reqs[i]);
		return 0;
	}
};

direct_comp_io::direct_comp_io(std::shared_ptr<remote_io> io): io_interface(
		io->get_thread()), req_buf(io->get_node_id(), 10240, true)
{
	comp_sched = comp_io_scheduler::ptr(new default_comp_io_scheduler(
				this->get_node_id()));
	comp_sched->set_io(this);
	underlying = io;
	underlying->set_callback(callback::ptr(new comp_callback(*this)));
	alloc_mem_size = 0;
	num_issued_areqs = 0;
	num_completed_areqs = 0;
	num_disk_bytes = 0;
	num_req_bytes = 0;
	arr_alloc = std::unique_ptr<direct_byte_array_allocator>(
			new direct_byte_array_allocator(io->get_thread()));
}

direct_comp_io::~direct_comp_io()
{
	cleanup();
}

void direct_comp_io::complete_req(const io_request &req)
{
	num_completed_areqs++;
	// When an I/O request is complete, we need to invoke the user task.
	orig_comp_request *orig = (orig_comp_request *) req.get_user_data();
	assert(orig);
	assert(orig->get_offset() + orig->get_size()
			<= req.get_offset() + req.get_size());
	// We pass the ownership of the memory buffer to the byte array.
	// The byte array will free the buffer once it's done.
	direct_byte_array arr(orig->get_offset(), orig->get_size(), req.get_buf(),
			*arr_alloc);
	orig->get_compute()->run(arr);
	comp_sched->post_comp_process(orig->get_compute());
	delete orig;

	alloc_mem_size -= req.get_size();
}

int direct_comp_io::get_file_id() const
{
	return underlying->get_file_id();
}

void direct_comp_io::cleanup()
{
	// wait4complete may generate more requests because of user compute
	// tasks. We have to make sure all requests are completed.
	while (num_pending_ios() > 0 || !comp_sched->is_empty())
		wait4complete(num_pending_ios());
	underlying->cleanup();
}

static void conv_comp_to_basic(io_request &req)
{
	// The data requested by the user compute may not be aligned to
	// BLOCK_SIZE, we need to align the offset and size.
	off_t align_off = ROUND(req.get_offset(), MIN_BLOCK_SIZE);
	off_t align_off_end = ROUNDUP(req.get_offset() + req.get_size(),
			MIN_BLOCK_SIZE);
	size_t req_size = align_off_end - align_off;
	char *buf = NULL;
	int ret = posix_memalign((void **) &buf, MIN_BLOCK_SIZE, req_size);
	assert(ret == 0);
	orig_comp_request *orig = new orig_comp_request(req.get_offset(),
			req.get_size(), req.get_compute());
	data_loc_t loc(req.get_file_id(), align_off);
	req = io_request(buf, loc, req_size, req.get_access_method(),
			req.get_io(), req.get_node_id());
	req.set_user_data(orig);
}

void direct_comp_io::access(io_request *requests, int num, io_status *status)
{
	num_issued_areqs += num;
	int i;
	for (i = 0; i < num; i++) {
		assert(requests[i].get_user_data() == NULL);
		// It has to be a user-compute request, and we need to convert it into
		// a basic I/O request.
		assert(requests[i].get_req_type() == io_request::USER_COMPUTE);
		user_compute *compute = requests[i].get_compute();
		compute->inc_ref();

		if (alloc_mem_size <= MAX_PEND_COMP_SIZE) {
			// We count the number of bytes accessed by I/O requests.
			num_req_bytes = requests[i].get_size();
			// This convert the original request to the request that accesses
			// the aligned area on the disks.
			conv_comp_to_basic(requests[i]);
			num_disk_bytes = requests[i].get_size();
			alloc_mem_size += requests[i].get_size();
		}
		else
			break;
	}
	// If there are requests that didn't get converted (because there are
	// too many pending requests), we need to buffer the rest of requests.
	if (i < num)
		req_buf.add(&requests[i], num - i);
	// We issue the requests that have been converted to the underlying I/O.
	underlying->access(requests, i, status);
}

void direct_comp_io::process_buf_reqs()
{
	const int BUF_SIZE = 8;
	io_request local_buf[BUF_SIZE];
	int num_reqs = 0;
	while (!req_buf.is_empty() && alloc_mem_size <= MAX_PEND_COMP_SIZE) {
		io_request req = req_buf.pop_front();
		assert(req.get_user_data() == NULL);
		assert(req.get_req_type() == io_request::USER_COMPUTE);

		conv_comp_to_basic(req);
		alloc_mem_size += req.get_size();
		// The IO instance in the request gets notified from the I/O thread
		// when the IO request is completed. direct_comp_io gets notification
		// from remote_io, so we should set io in the request to
		// the underlying IO instance.
		req.set_io(underlying.get());
		local_buf[num_reqs++] = req;

		if (num_reqs == BUF_SIZE) {
			underlying->access(local_buf, num_reqs, NULL);
			num_reqs = 0;
		}
	}

	if (num_reqs > 0) {
		underlying->access(local_buf, num_reqs, NULL);
	}
}

void direct_comp_io::process_incomplete_computes()
{
	if (alloc_mem_size >= MAX_PEND_COMP_SIZE)
		return;
	size_t num_new = comp_sched->get_requests(req_buf, req_buf.get_size());
	num_issued_areqs += num_new;
	comp_sched->gc_computes();
	if (!req_buf.is_empty())
		process_buf_reqs();
}

void direct_comp_io::flush_requests()
{
	process_buf_reqs();
	if (!req_buf.is_full())
		process_incomplete_computes();
	underlying->flush_requests();
}

int direct_comp_io::wait4complete(int num_to_complete)
{
	// Issue buffered requests to the underlying IO.
	flush_requests();
	size_t prev_completed_areqs = num_completed_areqs;
	underlying->wait4complete(num_to_complete);
	// We can now issue more buffered requests.
	flush_requests();
	return num_completed_areqs - prev_completed_areqs;
}

int direct_comp_io::num_pending_ios() const
{
	assert(num_issued_areqs >= num_completed_areqs);
	return num_issued_areqs - num_completed_areqs;
}

}
