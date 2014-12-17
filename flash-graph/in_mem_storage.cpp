/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashGraph.
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
#include <malloc.h>

#include <boost/format.hpp>

#include "log.h"
#include "safs_file.h"
#include "cache.h"
#include "slab_allocator.h"
#include "native_file.h"
#include "comp_io_scheduler.h"
#include "io_interface.h"

#include "in_mem_storage.h"
#include "graph_file_header.h"

using namespace safs;

namespace fg
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

class in_mem_io: public io_interface
{
	const in_mem_graph &graph;
	int file_id;
	fifo_queue<io_request> req_buf;
	comp_io_scheduler::ptr comp_io_sched;
	std::unique_ptr<byte_array_allocator> array_allocator;

	callback::ptr cb;

	void process_req(const io_request &req);
	void process_computes();
public:
	in_mem_io(const in_mem_graph &_graph, int file_id,
			thread *t): io_interface(t), graph(_graph), req_buf(
				get_node_id(), 1024) {
		this->file_id = file_id;
		array_allocator = std::unique_ptr<byte_array_allocator>(
				new in_mem_byte_array_allocator(t));
		comp_io_sched = comp_io_scheduler::ptr(
				new default_comp_io_scheduler(get_node_id()));
		comp_io_sched->set_io(this);
	}

	virtual int get_file_id() const {
		return file_id;
	}

	virtual bool support_aio() {
		return true;
	}

	virtual bool set_callback(callback::ptr cb) {
		this->cb = cb;
		return true;
	}

	virtual bool have_callback() const {
		return cb != NULL;
	}

	virtual callback &get_callback() {
		return *cb;
	}

	virtual void flush_requests() { }

	virtual int num_pending_ios() const {
		return 0;
	}

	virtual io_status access(char *buf, off_t off, ssize_t size,
			int access_method) {
		assert(access_method == READ);
		memcpy(buf, graph.graph_data + off, size);
		return IO_OK;
	}

	virtual void access(io_request *requests, int num, io_status *status);
	virtual int wait4complete(int) {
		return 0;
	}
};

enum
{
	IN_QUEUE,
};

void in_mem_io::process_req(const io_request &req)
{
	assert(req.get_req_type() == io_request::USER_COMPUTE);
	in_mem_byte_array byte_arr(req,
			graph.graph_data + ROUND_PAGE(req.get_offset()), *array_allocator);
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
			memcpy(req.get_buf(), graph.graph_data + req.get_offset(), req.get_size());
			io_request *reqs[1];
			reqs[0] = &req;
			if (this->have_callback())
				this->get_callback().invoke(reqs, 1);
		}
	}
	process_computes();
	comp_io_sched->gc_computes();
}

class in_mem_io_factory: public file_io_factory
{
	const in_mem_graph &graph;
	int file_id;
public:
	in_mem_io_factory(const in_mem_graph &_graph, int file_id,
			const std::string file_name): file_io_factory(file_name), graph(_graph) {
		this->file_id = file_id;
	}

	virtual int get_file_id() const {
		return file_id;
	}

	virtual io_interface::ptr create_io(thread *t) {
		return io_interface::ptr(new in_mem_io(graph, file_id, t));
	}

	virtual void destroy_io(io_interface *) {
	}
};

in_mem_graph::ptr in_mem_graph::load_graph(const std::string &file_name)
{
	native_file local_f(file_name);
	if (!local_f.exist())
		throw io_exception(file_name + std::string(" doesn't exist"));

	ssize_t size = local_f.get_size();
	assert(size > 0);

	in_mem_graph::ptr graph = in_mem_graph::ptr(new in_mem_graph());
	graph->graph_size = size;
	graph->graph_data = (char *) malloc(size);
	assert(graph->graph_data);
	graph->graph_file_name = file_name;
	BOOST_LOG_TRIVIAL(info) << boost::format("load a graph of %1% bytes")
		% graph->graph_size;

	FILE *fd = fopen(file_name.c_str(), "r");
	if (fd == NULL)
		throw io_exception(std::string("can't open ") + file_name);
	if (fread(graph->graph_data, size, 1, fd) != 1)
		throw io_exception(std::string("can't read from ") + file_name);
	fclose(fd);

	graph_header *header = (graph_header *) graph->graph_data;
	header->verify();

	return graph;
}

in_mem_graph::ptr in_mem_graph::load_safs_graph(const std::string &file_name)
{
	file_io_factory::shared_ptr io_factory = ::create_io_factory(file_name,
			REMOTE_ACCESS);

	in_mem_graph::ptr graph = in_mem_graph::ptr(new in_mem_graph());
	graph->graph_size = io_factory->get_file_size();
	size_t num_pages = ROUNDUP_PAGE(graph->graph_size) / PAGE_SIZE;
	int ret = posix_memalign((void **) &graph->graph_data, PAGE_SIZE,
			num_pages * PAGE_SIZE);
	BOOST_VERIFY(ret == 0);
	graph->graph_file_name = file_name;

	BOOST_LOG_TRIVIAL(info) << boost::format("load a graph of %1% bytes")
		% graph->graph_size;
#if 0
	graph->graph_file_id = io_factory->get_file_id();
#endif
	io_interface::ptr io = io_factory->create_io(thread::get_curr_thread());
	const size_t MAX_IO_SIZE = 256 * 1024 * 1024;
	for (off_t off = 0; (size_t) off < graph->graph_size; off += MAX_IO_SIZE) {
		data_loc_t loc(io_factory->get_file_id(), off);
		size_t req_size = min(MAX_IO_SIZE, graph->graph_size - off);
		io_request req(graph->graph_data + off, loc, req_size, READ);
		io->access(&req, 1);
		io->wait4complete(1);
	}

	graph_header *header = (graph_header *) graph->graph_data;
	header->verify();

	return graph;
}

file_io_factory::shared_ptr in_mem_graph::create_io_factory() const
{
	return file_io_factory::shared_ptr(new in_mem_io_factory(*this,
				graph_file_id, graph_file_name));
}

}
