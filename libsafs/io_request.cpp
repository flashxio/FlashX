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

#include "io_request.h"
#include "cache.h"
#include "io_interface.h"

namespace safs
{

void *io_buf::get_buf() const
{
	if (is_page)
		return u.p->get_data();
	else
		return u.buf;
}

void io_buf::init(thread_safe_page *p)
{
	assert(p->get_ref() > 0);
	u.p = p;
	size = PAGE_SIZE;
	is_page = 1;
}

void io_request::init(char *buf, const data_loc_t &loc, ssize_t size,
		int access_method, io_interface *io, int node_id)
{
	assert(loc.get_offset() <= MAX_FILE_SIZE);
	this->file_id = loc.get_file_id();
	this->offset = loc.get_offset();
	this->io_addr = (long) io;
	if (is_extended_req()) {
		if (buf)
			add_buf(buf, size);
	}
	else if (payload_type == BASIC_REQ) {
		set_int_buf_size(size);
		this->payload.buf_addr = buf;
	}
	else {
		set_int_buf_size(size);
		this->payload.compute = NULL;
	}
	this->access_method = access_method & 0x1;
	// by default, a request is of high priority.
	assert(node_id <= MAX_NODE_ID);
	this->node_id = node_id;
}

int io_request::get_overlap_size(thread_safe_page *pg) const
{
	off_t start = max(pg->get_offset(), this->get_offset());
	off_t end = min(pg->get_offset() + PAGE_SIZE,
			this->get_offset() + this->get_size());
	return end - start;
}

file_id_t io_request::get_file_id() const
{
	return file_id;
}

void io_req_extension::add_io_buf(const io_buf &buf)
{
	if (num_bufs >= vec_capacity) {
		if (vec_pointer == embedded_vecs) {
			vec_capacity = MIN_NUM_ALLOC_IOVECS;
			vec_pointer = new io_buf[vec_capacity];
			memcpy(vec_pointer, embedded_vecs,
					sizeof(embedded_vecs[0]) * NUM_EMBEDDED_IOVECS);
		}
		else {
			vec_capacity *= 2;
			io_buf *tmp = new io_buf[vec_capacity];
			memcpy(tmp, vec_pointer,
					sizeof(vec_pointer[0]) * vec_capacity / 2);
			delete [] vec_pointer;
			vec_pointer = tmp;
		}
	}
	assert(num_bufs < vec_capacity);
	vec_pointer[num_bufs] = buf;
	num_bufs++;
}

void io_req_extension::add_buf(char *buf, int size, bool is_page)
{
	io_buf tmp;
	tmp.init(buf, size, is_page);
	add_io_buf(tmp);
}

void io_req_extension::add_buf_front(char *buf, int size, bool is_page)
{
	if (num_bufs >= vec_capacity) {
		if (vec_pointer == embedded_vecs) {
			vec_capacity = MIN_NUM_ALLOC_IOVECS;
			vec_pointer = new io_buf[vec_capacity];
			memcpy(vec_pointer + 1, embedded_vecs,
					sizeof(embedded_vecs[0]) * NUM_EMBEDDED_IOVECS);
		}
		else {
			vec_capacity *= 2;
			io_buf *tmp = new io_buf[vec_capacity];
			memcpy(tmp + 1, vec_pointer,
					sizeof(vec_pointer[0]) * vec_capacity / 2);
			delete [] vec_pointer;
			vec_pointer = tmp;
		}
	}
	else {
		memmove(vec_pointer + 1, vec_pointer,
				sizeof(vec_pointer[0]) * num_bufs);
	}
	assert(num_bufs < vec_capacity);
	vec_pointer[0].init((void *) buf, size, is_page);
	num_bufs++;
}

int user_compute::fetch_requests(io_interface *io, user_comp_req_queue &reqs,
		int max_fetch)
{
	int num_issues = 0;
	while (has_requests() && num_issues < max_fetch) {
		request_range range = get_next_request();
		user_compute *compute = range.get_compute();
		compute->inc_ref();
		io_request req(compute, range.get_loc(), range.get_size(),
				range.get_access_method(), io, io->get_node_id());
		reqs.push_back(req);
		num_issues++;
	}
	return num_issues;
}

bool user_compute::fetch_request(io_interface *io, io_request &req)
{
	if (!has_requests())
		return false;

	request_range range = get_next_request();
	user_compute *compute = range.get_compute();
	compute->inc_ref();
	req = io_request(compute, range.get_loc(), range.get_size(),
			range.get_access_method(), io, io->get_node_id());
	return true;
}

}
