/*
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

#include "in_mem_io.h"
#include "slab_allocator.h"

namespace safs
{

class in_mem_byte_array: public page_byte_array
{
	off_t off;
	size_t size;
	const char *pages;

	void assign(in_mem_byte_array &arr) {
		this->off = arr.off;
		this->size = arr.size;
		this->pages = arr.pages;
	}

	in_mem_byte_array(in_mem_byte_array &arr) {
		assign(arr);
	}

	in_mem_byte_array &operator=(in_mem_byte_array &arr) {
		assign(arr);
		return *this;
	}
public:
	in_mem_byte_array(byte_array_allocator &alloc): page_byte_array(alloc) {
		off = 0;
		size = 0;
		pages = NULL;
	}

	in_mem_byte_array(const io_request &req, char *pages,
			byte_array_allocator &alloc): page_byte_array(alloc) {
		this->off = req.get_offset();
		this->size = req.get_size();
		this->pages = pages;
	}

	virtual off_t get_offset() const {
		return off;
	}

	virtual off_t get_offset_in_first_page() const {
		return off % PAGE_SIZE;
	}

	virtual const char *get_page(int pg_idx) const {
		return pages + pg_idx * PAGE_SIZE;
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
		in_mem_byte_array *arr = (in_mem_byte_array *) get_allocator().alloc();
		*arr = *this;
		return arr;
	}
};

class in_mem_byte_array_allocator: public byte_array_allocator
{
	class array_initiator: public obj_initiator<in_mem_byte_array>
	{
		in_mem_byte_array_allocator *alloc;
	public:
		array_initiator(in_mem_byte_array_allocator *alloc) {
			this->alloc = alloc;
		}

		virtual void init(in_mem_byte_array *obj) {
			new (obj) in_mem_byte_array(*alloc);
		}
	};

	class array_destructor: public obj_destructor<in_mem_byte_array>
	{
	public:
		void destroy(in_mem_byte_array *obj) {
			obj->~in_mem_byte_array();
		}
	};

	obj_allocator<in_mem_byte_array> allocator;
public:
	in_mem_byte_array_allocator(thread *t): allocator(
			"byte-array-allocator", t->get_node_id(), false, 1024 * 1024,
			params.get_max_obj_alloc_size(),
			obj_initiator<in_mem_byte_array>::ptr(new array_initiator(this)),
			obj_destructor<in_mem_byte_array>::ptr(new array_destructor())) {
	}

	virtual page_byte_array *alloc() {
		return allocator.alloc_obj();
	}

	virtual void free(page_byte_array *arr) {
		allocator.free((in_mem_byte_array *) arr);
	}
};

void in_mem_io::process_req(const io_request &req)
{
	assert(req.get_req_type() == io_request::USER_COMPUTE);
	in_mem_byte_array byte_arr(req,
			data.get() + ROUND_PAGE(req.get_offset()), *array_allocator);
	user_compute *compute = req.get_compute();
	compute->run(byte_arr);
	comp_io_sched->post_comp_process(compute);
}

void in_mem_io::process_computes()
{
	while (true) {
		comp_io_sched->get_requests(req_buf, req_buf.get_size());
		if (req_buf.is_empty())
			break;

		while (!req_buf.is_empty()) {
			io_request new_req = req_buf.pop_front();
			process_req(new_req);
		}
	}
}

void in_mem_io::access(io_request *requests, int num, io_status *)
{
	for (int i = 0; i < num; i++) {
		io_request &req = requests[i];
		if (req.get_req_type() == io_request::USER_COMPUTE) {
			// Let's possess a reference to the user compute first. process_req()
			// will release the reference when the user compute is completed.
			req.get_compute()->inc_ref();
			process_req(req);
		}
		else {
			assert(req.get_req_type() == io_request::BASIC_REQ);
			if (req.get_access_method() == READ)
				memcpy(req.get_buf(), data.get() + req.get_offset(), req.get_size());
			else
				memcpy(data.get() + req.get_offset(), req.get_buf(), req.get_size());
			io_request *reqs[1];
			reqs[0] = &req;
			if (this->have_callback())
				this->get_callback().invoke(reqs, 1);
		}
	}
	process_computes();
	comp_io_sched->gc_computes();
}

in_mem_io::in_mem_io(std::shared_ptr<char> data, int file_id,
		thread *t): io_interface(t), req_buf(get_node_id(), 1024)
{
	this->data = data;
	this->file_id = file_id;
	array_allocator = std::unique_ptr<byte_array_allocator>(
			new in_mem_byte_array_allocator(t));
	comp_io_sched = comp_io_scheduler::ptr(
			new default_comp_io_scheduler(get_node_id()));
	comp_io_sched->set_io(this);
}

}
