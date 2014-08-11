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

#include "safs_file.h"
#include "in_mem_storage.h"
#include "cache.h"

class in_mem_io: public io_interface
{
	const in_mem_graph &graph;
	int file_id;
public:
	in_mem_io(const in_mem_graph &_graph, int file_id,
			thread *t): io_interface(t), graph(_graph) {
		this->file_id = file_id;
	}

	virtual int get_file_id() const {
		return file_id;
	}

	virtual bool support_aio() {
		return true;
	}

	virtual void access(io_request *requests, int num, io_status *status);
	virtual void flush_requests() { }
	virtual int wait4complete(int) {
		return 0;
	}
	virtual int num_pending_ios() const {
		return 0;
	}
};

class in_mem_byte_array: public page_byte_array
{
	io_request &req;
	thread_safe_page *pages;
public:
	in_mem_byte_array(io_request &_req, thread_safe_page *pages): req(_req) {
		this->pages = pages;
	}

	virtual off_t get_offset_in_first_page() const {
		return req.get_offset() % PAGE_SIZE;
	}

	virtual thread_safe_page *get_page(int pg_idx) const {
		return &pages[pg_idx];
	}

	virtual size_t get_size() const {
		return req.get_size();
	}

	void lock() {
		assert(0);
	}

	void unlock() {
		assert(0);
	}
};

void in_mem_io::access(io_request *requests, int num, io_status *)
{
	for (int i = 0; i < num; i++) {
		io_request &req = requests[i];
		assert(req.get_req_type() == io_request::USER_COMPUTE);
		in_mem_byte_array byte_arr(req, &graph.graph_pages[req.get_offset() / PAGE_SIZE]);
		req.get_compute()->run(byte_arr);
	}
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
	file_io_factory::shared_ptr io_factory = ::create_io_factory(file_name,
			REMOTE_ACCESS);

	in_mem_graph::ptr graph = in_mem_graph::ptr(new in_mem_graph());
	graph->graph_size = io_factory->get_file_size();
	size_t num_pages = ROUNDUP_PAGE(graph->graph_size) / PAGE_SIZE;
	graph->graph_data = (char *) memalign(PAGE_SIZE, num_pages * PAGE_SIZE);
	assert(graph->graph_data);
	graph->graph_file_name = file_name;
	graph->graph_pages = new thread_safe_page[num_pages];
	assert(graph->graph_pages);

	graph->graph_file_id = io_factory->get_file_id();
	io_interface::ptr io = io_factory->create_io(thread::get_curr_thread());
	data_loc_t loc(io_factory->get_file_id(), 0);
	io_request req(graph->graph_data, loc, graph->graph_size, READ);
	io->access(&req, 1);
	io->wait4complete(1);

	for (size_t i = 0; i < num_pages; i++) {
		page_id_t pg_id(graph->graph_file_id, i * PAGE_SIZE);
		graph->graph_pages[i] = thread_safe_page(pg_id,
				graph->graph_data + i * PAGE_SIZE, 0);
	}

	return graph;
}

file_io_factory::shared_ptr in_mem_graph::create_io_factory() const
{
	return file_io_factory::shared_ptr(new in_mem_io_factory(*this,
				graph_file_id, graph_file_name));
}
