#ifndef __VERTEX_INDEX_READER__
#define __VERTEX_INDEX_READER__

/**
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

#include "FG_basic_types.h"
#include "io_interface.h"
#include "simple_KV_store.h"
#include "vertex.h"

class vertex_compute;
class directed_vertex_compute;

template<class ValueType>
class req_vertex_task
{
	vertex_id_t vid;
	vertex_compute *compute;
public:
	req_vertex_task() {
		vid = 0;
		compute = NULL;
	}

	req_vertex_task(vertex_id_t vid, vertex_compute &_compute) {
		this->vid = vid;
		this->compute = &_compute;
	}

	size_t get_idx() const {
		return vid + vertex_index::get_header_size() / sizeof(ValueType);
	}

	size_t get_num_entries() const {
		return 2;
	}

	void run(ValueType entries[], int num) {
		assert(num == 2);
#ifdef VERIFY_INDEX_READER
		extern vertex_index *test_vindex;
		vertex_index_temp<ValueType> *index = (vertex_index_temp<ValueType> *) test_vindex;
		assert(entries[0].get_off() == index->get_vertex(vid).get_off());
		if (vid < index->get_max_id())
			assert(entries[1].get_off() == index->get_vertex(vid + 1).get_off());
#endif
		in_mem_vertex_info info(vid, entries[0].get_off(),
				entries[1].get_off() - entries[0].get_off());
		compute->issue_io_request(info);
	}

	// Override the operator so it can be used in a priority queue.
	bool operator<(const req_vertex_task &task) const {
		return this->get_idx() > task.get_idx();
	}
};

class req_part_vertex_task
{
	directed_vertex_request req;
	directed_vertex_compute *compute;
public:
	req_part_vertex_task() {
		compute = NULL;
	}

	req_part_vertex_task(directed_vertex_request &req,
			directed_vertex_compute &_compute) {
		this->req = req;
		this->compute = &_compute;
	}

	size_t get_idx() const {
		return req.get_id()
			+ vertex_index::get_header_size() / sizeof(directed_vertex_entry);
	}

	size_t get_num_entries() const {
		return 2;
	}

	void run(directed_vertex_entry entries[], int num) {
		assert(num == 2);
#ifdef VERIFY_INDEX_READER
		extern vertex_index *test_vindex;
		directed_vertex_index *index = (directed_vertex_index *) test_vindex;
		assert(entries[0].get_off() == index->get_vertex(req.get_id()).get_off());
		if (req.get_id() < index->get_max_id())
			assert(entries[1].get_off() == index->get_vertex(req.get_id() + 1).get_off());
#endif
		in_mem_directed_vertex_info info(req.get_id(), entries[0].get_off(),
				entries[1].get_off() - entries[0].get_off(),
				entries[0].get_num_in_edges(), entries[0].get_num_out_edges());
		compute->issue_io_request(req, info);
	}
};

template<class ValueType>
class req_edge_task
{
	// The vertex whose edges the task requests.
	vertex_id_t vid;
	// The vertex compute of the vertex who issues the request.
	vertex_compute *compute;
public:
	req_edge_task() {
		vid = 0;
		compute = NULL;
	}

	req_edge_task(vertex_id_t vid, vertex_compute *compute) {
		this->vid = vid;
		this->compute = compute;
	}

	size_t get_idx() const {
		return vid + vertex_index::get_header_size() / sizeof(ValueType);
	}

	size_t get_num_entries() const {
		return 2;
	}

	void run(ValueType entries[], int num) {
		assert(num == 2);
#ifdef VERIFY_INDEX_READER
		extern vertex_index *test_vindex;
		vertex_index_temp<ValueType> *index = (vertex_index_temp<ValueType> *) test_vindex;
		assert(entries[0].get_off() == index->get_vertex(vid).get_off());
		if (vid < index->get_max_id())
			assert(entries[1].get_off() == index->get_vertex(vid + 1).get_off());
#endif
		compute->run_on_vertex_size(vid,
				entries[1].get_off() - entries[0].get_off());
	}

	// Override the operator so it can be used in a priority queue.
	bool operator<(const req_edge_task &task) const {
		return this->get_idx() > task.get_idx();
	}
};

class req_directed_edge_task
{
	// The vertex whose edges the task requests.
	vertex_id_t vid;
	// The vertex compute of the vertex who issues the request.
	directed_vertex_compute *compute;
public:
	req_directed_edge_task() {
		vid = 0;
		compute = NULL;
	}

	req_directed_edge_task(vertex_id_t vid, directed_vertex_compute *compute) {
		this->vid = vid;
		this->compute = compute;
	}

	size_t get_idx() const {
		return vid + vertex_index::get_header_size() / sizeof(directed_vertex_entry);
	}

	size_t get_num_entries() const {
		return 1;
	}

	void run(directed_vertex_entry entries[], int num) {
		assert(num == 1);
#ifdef VERIFY_INDEX_READER
		extern vertex_index *test_vindex;
		directed_vertex_index *index = (directed_vertex_index *) test_vindex;
		assert(entries[0].get_off() == index->get_vertex(vid).get_off());
		if (vid < index->get_max_id())
			assert(entries[1].get_off() == index->get_vertex(vid + 1).get_off());
#endif
		compute->run_on_num_edges(vid, entries[0].get_num_in_edges(),
				entries[0].get_num_out_edges());
	}

	// Override the operator so it can be used in a priority queue.
	bool operator<(const req_directed_edge_task &task) const {
		return this->get_idx() > task.get_idx();
	}
};

/*
 * This interface reads vertex index from SSDs.
 * It accepts the requests of reading vertices or partial vertices as well
 * as other vertex information in the vertex index such as the number of edges.
 * This is a per-thread data structure.
 */
class vertex_index_reader
{
public:
	typedef std::shared_ptr<vertex_index_reader> ptr;

	virtual void request_vertices(vertex_id_t ids[], size_t num,
			vertex_compute &compute) = 0;
	virtual void request_num_edges(vertex_id_t vertices[], size_t num,
			vertex_compute &compute) = 0;
	virtual void request_num_directed_edges(vertex_id_t ids[], size_t num,
			directed_vertex_compute &compute) = 0;
	virtual void request_vertices(directed_vertex_request reqs[], size_t num,
			directed_vertex_compute &compute) = 0;
	virtual void wait4complete(int num) = 0;
	virtual size_t get_num_pending_tasks() const = 0;
};

template<class ValueType>
class ext_mem_vindex_reader_impl: public vertex_index_reader
{
	io_interface::ptr io;
	typedef simple_KV_store<ValueType, req_vertex_task<ValueType> > vertex_KV_store;
	typedef simple_KV_store<ValueType, req_edge_task<ValueType> > edge_KV_store;
	typename vertex_KV_store::ptr req_vertex_store;
	typename edge_KV_store::ptr req_edge_store;

	// The starting offset of the entries in the cached index.
	size_t cached_index_start;
	// The values in the last one or two pages in the index.
	std::vector<ValueType> cached_index;

protected:
	ext_mem_vindex_reader_impl(io_interface::ptr io) {
		this->io = io;
		req_vertex_store = vertex_KV_store::create(io);
		req_edge_store = edge_KV_store::create(io);

		// Get the vertex index header.
		char *hdr_buf = new char[vertex_index::get_header_size()];
		io->access(hdr_buf, 0, vertex_index::get_header_size(), READ);
		vertex_index *index = (vertex_index *) hdr_buf;
		vsize_t num_vertices = index->get_num_vertices();
		size_t graph_size = index->get_graph_size();
		delete [] hdr_buf;

		// Keep the last one or two pages in memory.
		off_t read_start;
		size_t index_size = num_vertices * sizeof(ValueType)
			+ vertex_index::get_header_size();
		if (num_vertices * sizeof(ValueType) >= PAGE_SIZE)
			read_start = index_size - PAGE_SIZE;
		else
			read_start = vertex_index::get_header_size();
		assert(PAGE_SIZE % sizeof(ValueType) == 0);
		// The index contains a header.
		cached_index_start = (read_start
				- vertex_index::get_header_size()) / sizeof(ValueType);
		// We keep all entries in the index after `cached_index_start`
		// in the vector. We also add another entry to tell the size of
		// the graph image.
		cached_index.resize((index_size - read_start) / sizeof(ValueType) + 1);
		io->access((char *) cached_index.data(), read_start,
				(index_size - read_start), READ);
		cached_index.back() = ValueType(graph_size);
	}

	size_t get_cached_index_start() const {
		return cached_index_start;
	}

	const ValueType &get_cached_entry(int idx) const {
		return cached_index[idx];
	}
public:
	static ptr create(io_interface::ptr io) {
		return ptr(new ext_mem_vindex_reader_impl(io));
	}

	void request_vertices(vertex_id_t ids[], size_t num, vertex_compute &compute) {
		for (size_t i = 0; i < num; i++) {
			req_vertex_task<ValueType> task(ids[i], compute);
			if (ids[i] >= get_cached_index_start()) {
				int off_in_cache = ids[i] - get_cached_index_start();
				int num_entries = task.get_num_entries();
				ValueType vs[num_entries];
				for (int j = 0; j < num_entries; j++)
					vs[j] = get_cached_entry(off_in_cache + j);
				task.run(vs, num_entries);
			}
			else
				req_vertex_store->async_request(task);
		}
	}

	void request_num_edges(vertex_id_t ids[], size_t num,
			vertex_compute &compute) {
		for (size_t i = 0; i < num; i++) {
			req_edge_task<ValueType> task(ids[i], &compute);
			if (ids[i] >= get_cached_index_start()) {
				int off_in_cache = ids[i] - get_cached_index_start();
				int num_entries = task.get_num_entries();
				ValueType vs[num_entries];
				for (int j = 0; j < num_entries; j++)
					vs[j] = get_cached_entry(off_in_cache + j);
				task.run(vs, num_entries);
			}
			else
				req_edge_store->async_request(task);
		}
	}

	virtual void request_num_directed_edges(vertex_id_t ids[], size_t num,
			directed_vertex_compute &compute) {
		assert(0);
	}

	virtual void request_vertices(directed_vertex_request reqs[], size_t num,
			directed_vertex_compute &compute) {
		assert(0);
	}

	void wait4complete(int num) {
		req_vertex_store->flush_requests();
		req_edge_store->flush_requests();
		io->wait4complete(num);
	}

	size_t get_num_pending_tasks() const {
		return req_vertex_store->get_num_pending_tasks()
			+ req_edge_store->get_num_pending_tasks();
	}
};

typedef ext_mem_vindex_reader_impl<vertex_offset> ext_mem_undirected_vindex_reader;

class ext_mem_directed_vindex_reader: public ext_mem_vindex_reader_impl<directed_vertex_entry>
{
	typedef simple_KV_store<directed_vertex_entry,
			req_part_vertex_task> part_vertex_KV_store;
	part_vertex_KV_store::ptr req_part_vertex_store;
	typedef simple_KV_store<directed_vertex_entry,
			req_directed_edge_task> directed_edge_KV_store;
	directed_edge_KV_store::ptr req_directed_edge_store;

	ext_mem_directed_vindex_reader(
			io_interface::ptr io): ext_mem_vindex_reader_impl<directed_vertex_entry>(io) {
		req_part_vertex_store = part_vertex_KV_store::create(io);
		req_directed_edge_store = directed_edge_KV_store::create(io);
	}
public:
	static ptr create(io_interface::ptr io) {
		return ptr(new ext_mem_directed_vindex_reader(io));
	}

	virtual void request_vertices(directed_vertex_request reqs[], size_t num,
			directed_vertex_compute &compute) {
		for (size_t i = 0; i < num; i++) {
			req_part_vertex_task task(reqs[i], compute);
			if (reqs[i].get_id() >= get_cached_index_start()) {
				int off_in_cache = reqs[i].get_id() - get_cached_index_start();
				int num_entries = task.get_num_entries();
				directed_vertex_entry vs[num_entries];
				for (int j = 0; j < num_entries; j++)
					vs[j] = get_cached_entry(off_in_cache + j);
				task.run(vs, num_entries);
			}
			else
				req_part_vertex_store->async_request(task);
		}
	}

	void request_num_directed_edges(vertex_id_t ids[], size_t num,
			directed_vertex_compute &compute) {
		for (size_t i = 0; i < num; i++) {
			req_directed_edge_task task(ids[i], &compute);
			if (ids[i] >= get_cached_index_start()) {
				int off_in_cache = ids[i] - get_cached_index_start();
				int num_entries = task.get_num_entries();
				directed_vertex_entry vs[num_entries];
				for (int j = 0; j < num_entries; j++)
					vs[j] = get_cached_entry(off_in_cache + j);
				task.run(vs, num_entries);
			}
			else
				req_directed_edge_store->async_request(task);
		}
	}

	void wait4complete(int num) {
		req_part_vertex_store->flush_requests();
		req_directed_edge_store->flush_requests();
		ext_mem_vindex_reader_impl<directed_vertex_entry>::wait4complete(num);
	}

	size_t get_num_pending_tasks() const {
		return ext_mem_vindex_reader_impl<directed_vertex_entry>::get_num_pending_tasks()
			+ req_part_vertex_store->get_num_pending_tasks();
	}
};

template<class ValueType>
class in_mem_vindex_reader_impl: public vertex_index_reader
{
	typename vertex_index_temp<ValueType>::ptr index;

protected:
	in_mem_vindex_reader_impl(vertex_index::ptr index) {
		this->index = vertex_index_temp<ValueType>::cast(index);
	}
public:
	static ptr create(vertex_index::ptr index) {
		return ptr(new in_mem_vindex_reader_impl<ValueType>(index));
	}


	virtual void request_vertices(vertex_id_t ids[], size_t num,
			vertex_compute &compute) {
		for (size_t i = 0; i < num; i++) {
			vertex_id_t id = ids[i];
			in_mem_vertex_info info(id, index->get_vertex_off(id),
					index->get_vertex_size(id));
			compute.issue_io_request(info);
		}
	}

	virtual void request_num_edges(vertex_id_t ids[], size_t num,
			vertex_compute &compute) {
		for (size_t i = 0; i < num; i++) {
			vertex_id_t id = ids[i];
			compute.run_on_vertex_size(id, index->get_vertex_size(id));
		}
	}

	virtual void request_num_directed_edges(vertex_id_t ids[], size_t num,
			directed_vertex_compute &compute) {
		assert(0);
	}

	virtual void request_vertices(directed_vertex_request reqs[], size_t num,
			directed_vertex_compute &compute) {
		assert(0);
	}

	virtual void wait4complete(int num) {
	}

	virtual size_t get_num_pending_tasks() const {
		return 0;
	}
};

typedef in_mem_vindex_reader_impl<vertex_offset> in_mem_undirected_vindex_reader;

class in_mem_directed_vindex_reader: public in_mem_vindex_reader_impl<directed_vertex_entry>
{
	directed_vertex_index::ptr index;

	in_mem_directed_vindex_reader(
			vertex_index::ptr index): in_mem_vindex_reader_impl<directed_vertex_entry>(
				index) {
		this->index = directed_vertex_index::cast(index);
	}
public:
	static ptr create(vertex_index::ptr index) {
		return ptr(new in_mem_directed_vindex_reader(index));
	}

	virtual void request_num_directed_edges(vertex_id_t ids[], size_t num,
			directed_vertex_compute &compute) {
		for (size_t i = 0; i < num; i++) {
			vertex_id_t id = ids[i];
			compute.run_on_num_edges(id, index->get_num_in_edges(id),
					index->get_num_out_edges(id));
		}
	}

	virtual void request_vertices(directed_vertex_request reqs[], size_t num,
			directed_vertex_compute &compute) {
		for (size_t i = 0; i < num; i++) {
			vertex_id_t id = reqs[i].get_id();
			in_mem_directed_vertex_info info(id, index->get_vertex_off(id),
					index->get_vertex_size(id), index->get_num_in_edges(id),
					index->get_num_out_edges(id));
			compute.issue_io_request(reqs[i], info);
		}
	}
};

#endif
