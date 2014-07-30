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

#include "io_interface.h"
#include "simple_KV_store.h"
#include "vertex_compute.h"
#include "vertex_index_reader.h"

template<class ValueType>
class req_vertex_task
{
	vertex_id_t start_vid;
	int num_vertices;
	index_compute *compute;
public:
	req_vertex_task() {
		compute = NULL;
	}

	req_vertex_task(vertex_id_t id, int num_vertices, index_compute *_compute) {
		this->start_vid = id;
		this->num_vertices = num_vertices;
		this->compute = _compute;
	}

	size_t get_idx() const {
		return start_vid
			+ vertex_index::get_header_size() / sizeof(ValueType);
	}

	size_t get_num_entries() const {
		return num_vertices + 1;
	}

	bool merge(const req_vertex_task<ValueType> &task) {
		return false;
	}

	void run(page_byte_array::seq_const_iterator<ValueType> &seq_it) {
		page_index_iterator_impl<ValueType> it(seq_it);
		bool ret = compute->run(compute->get_first_vertex(), it);
		if (ret)
			compute->get_allocator().free(compute);
	}

	// Override the operator so it can be used in a priority queue.
	bool operator<(const req_vertex_task &task) const {
		return this->get_idx() > task.get_idx();
	}
};

template<class ValueType>
class ext_mem_vindex_reader_impl: public vertex_index_reader
{
	io_interface::ptr io;
	typedef simple_KV_store<ValueType, req_vertex_task<ValueType> > vertex_KV_store;
	typename vertex_KV_store::ptr req_vertex_store;

	// The starting offset of the entries in the cached index.
	size_t cached_index_start;
	// The values in the last one or two pages in the index.
	std::vector<ValueType> cached_index;

protected:
	embedded_array<ValueType> value_buf;

	ext_mem_vindex_reader_impl(io_interface::ptr io) {
		this->io = io;
		req_vertex_store = vertex_KV_store::create(io);

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

	virtual void request_index(index_compute *compute) {
		index_compute::id_range_t range = compute->get_range();
		if (range.first >= cached_index_start) {
			off_t start_off = range.first - cached_index_start;
			off_t end_off = range.second - cached_index_start;
			array_index_iterator_impl<ValueType> it(cached_index.data() + start_off,
					cached_index.data() + end_off + 1);
			bool ret = compute->run(compute->get_first_vertex(), it);
			assert(ret);
			compute->get_allocator().free(compute);
		}
		else if (range.second > cached_index_start) {
			off_t end_off = range.second - cached_index_start;
			array_index_iterator_impl<ValueType> it(cached_index.data(),
					cached_index.data() + end_off + 1);
			compute->run(cached_index_start, it);
			req_vertex_task<ValueType> task(compute->get_first_vertex(),
					cached_index_start - compute->get_first_vertex(), compute);
			req_vertex_store->async_request(task);
		}
		else {
			req_vertex_task<ValueType> task(compute->get_first_vertex(),
					compute->get_last_vertex() - compute->get_first_vertex() + 1,
					compute);
			req_vertex_store->async_request(task);
		}
	}

	void wait4complete(int num) {
		req_vertex_store->flush_requests();
		if (get_num_pending_tasks() > 0)
			io->wait4complete(num);
	}

	size_t get_num_pending_tasks() const {
		return req_vertex_store->get_num_pending_tasks();
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


	virtual void request_index(index_compute *compute) {
		index_compute::id_range_t range = compute->get_range();
		if (range.second < index->get_num_vertices()) {
			array_index_iterator_impl<ValueType> it(index->get_data() + range.first,
					// We need an additional entry.
					index->get_data() + range.second + 1);
			bool ret = compute->run(compute->get_first_vertex(), it);
			assert(ret);
		}
		else {
			array_index_iterator_impl<ValueType> it(index->get_data() + range.first,
					index->get_data() + range.second);
			compute->run(compute->get_first_vertex(), it);

			// For the last vertex,
			ValueType vs[2];
			vs[0] = index->get_data()[index->get_num_vertices() - 1];
			vs[1] = ValueType(index->get_graph_size());
			it = array_index_iterator_impl<ValueType>(vs, &vs[2]);
			bool ret = compute->run(index->get_num_vertices() - 1, it);
			assert(ret);
		}
		compute->get_allocator().free(compute);
	}

	virtual void wait4complete(int num) {
	}

	virtual size_t get_num_pending_tasks() const {
		return 0;
	}
};

vertex_index_reader::ptr vertex_index_reader::create(vertex_index::ptr index,
		bool directed)
{
	if (directed)
		return in_mem_vindex_reader_impl<directed_vertex_entry>::create(index);
	else
		return in_mem_vindex_reader_impl<vertex_offset>::create(index);
}

vertex_index_reader::ptr vertex_index_reader::create(io_interface::ptr io,
		bool directed)
{
	if (directed)
		return ext_mem_vindex_reader_impl<directed_vertex_entry>::create(io);
	else
		return ext_mem_vindex_reader_impl<vertex_offset>::create(io);
}

bool req_vertex_compute::run(vertex_id_t vid, index_iterator &it)
{
	while (it.has_next()) {
		num_gets++;
		in_mem_vertex_info info(vid, it.get_curr_off(),
				it.get_curr_vertex_size());
		get_compute(vid)->issue_io_request(info);
		vid++;
		it.move_next();
	}
	return num_gets == get_num_vertices();
}

bool genrq_vertex_compute::run(vertex_id_t start_vid, index_iterator &it)
{
	int num_vertices = get_num_vertices();
	for (int i = get_start_idx(start_vid); i < num_vertices; i++) {
		vertex_id_t id = get_id(i);
		vertex_compute *compute = get_compute(i);
		assert(id >= start_vid);
		if (!it.move_to(id - start_vid)) {
			break;
		}
		in_mem_vertex_info info(id, it.get_curr_off(),
				it.get_curr_vertex_size());
		compute->issue_io_request(info);
		num_gets++;
	}
	return num_gets == get_num_vertices();
}

bool req_part_vertex_compute::run(vertex_id_t vid, index_iterator &it)
{
	while (it.has_next()) {
		num_gets++;

		in_mem_directed_vertex_info info(vid, it.get_curr_off(),
				it.get_curr_vertex_size(), it.get_curr_num_in_edges(),
				it.get_curr_num_out_edges());
		directed_vertex_request req(vid, type);
		directed_vertex_compute *dcompute
			= (directed_vertex_compute *) get_compute(vid);
		dcompute->issue_io_request(req, info);
		vid++;
		it.move_next();
	}
	return num_gets == get_num_vertices();
}

bool genrq_part_vertex_compute::run(vertex_id_t start_vid, index_iterator &it)
{
	int num_vertices = get_num_vertices();
	for (int i = get_start_idx(start_vid); i < num_vertices; i++) {
		vertex_id_t id = get_id(i);
		directed_vertex_compute *compute
			= (directed_vertex_compute *) get_compute(i);
		assert(id >= start_vid);
		if (!it.move_to(id - start_vid)) {
			break;
		}

		in_mem_directed_vertex_info info(id, it.get_curr_off(),
				it.get_curr_vertex_size(), it.get_curr_num_in_edges(),
				it.get_curr_num_out_edges());
		directed_vertex_request req(id, type);
		compute->issue_io_request(req, info);
		num_gets++;
	}
	return num_gets == get_num_vertices();
}

bool req_edge_compute::run(vertex_id_t vid, index_iterator &it)
{
	while (it.has_next()) {
		num_gets++;
		get_compute(vid)->run_on_vertex_size(vid, it.get_curr_vertex_size());
		vid++;
		it.move_next();
	}
	return num_gets == get_num_vertices();
}

bool genrq_edge_compute::run(vertex_id_t vid, index_iterator &it)
{
	assert(0);
}

bool req_directed_edge_compute::run(vertex_id_t vid, index_iterator &it)
{
	while (it.has_next()) {
		num_gets++;
		directed_vertex_compute *dcompute
			= (directed_vertex_compute *) get_compute(vid);
		dcompute->run_on_num_edges(vid, it.get_curr_num_in_edges(),
				it.get_curr_num_out_edges());
		vid++;
		it.move_next();
	}
	return num_gets == get_num_vertices();
}

bool genrq_directed_edge_compute::run(vertex_id_t start_vid, index_iterator &it)
{
	int num_vertices = get_num_vertices();
	for (int i = get_start_idx(start_vid); i < num_vertices; i++) {
		vertex_id_t id = get_id(i);
		directed_vertex_compute *compute
			= (directed_vertex_compute *) get_compute(i);
		assert(id >= start_vid);
		if (!it.move_to(id - start_vid)) {
			break;
		}

		compute->run_on_num_edges(id, it.get_curr_num_in_edges(),
				it.get_curr_num_out_edges());
		num_gets++;
	}
	return num_gets == get_num_vertices();
}

template<class RequestType, class ComputeType, class AllocatorType>
void process_requests(std::vector<RequestType> &reqs,
		AllocatorType &alloc, vertex_index_reader &index_reader)
{
	assert(!reqs.empty());
	ComputeType *compute = (ComputeType *) alloc.alloc();
	const RequestType &req = reqs[0];
	compute->init(req.first, req.second);
	for (size_t i = 1; i < reqs.size(); i++) {
		const RequestType &req = reqs[i];

		if (!compute->add_vertex(req.first, req.second)) {
			index_reader.request_index(compute);
			compute = (ComputeType *) alloc.alloc();
			compute->init(req.first, req.second);
		}
	}
	assert(!compute->empty());
	index_reader.request_index(compute);
	reqs.clear();
}

void simple_index_reader::flush_computes()
{
	if (!vertex_comps.empty()) {
		if (!std::is_sorted(vertex_comps.begin(), vertex_comps.end(),
					id_compute_less()))
			std::sort(vertex_comps.begin(), vertex_comps.end(), id_compute_less());
		if (is_dense(vertex_comps))
			process_requests<id_compute_t, req_vertex_compute,
				index_comp_allocator_impl<req_vertex_compute> >(vertex_comps,
					req_vertex_comp_alloc, *index_reader);
		else
			process_requests<id_compute_t, genrq_vertex_compute,
				general_index_comp_allocator_impl<genrq_vertex_compute> >(
						vertex_comps, genrq_vertex_comp_alloc, *index_reader);
	}

	for (int type = edge_type::IN_EDGE; type < edge_type::NUM_TYPES; type++) {
		if (!part_vertex_comps[type].empty()) {
			if (!std::is_sorted(part_vertex_comps[type].begin(),
						part_vertex_comps[type].end(), directed_compute_less()))
				std::sort(part_vertex_comps[type].begin(),
						part_vertex_comps[type].end(), directed_compute_less());
			if (is_dense(part_vertex_comps[type]))
				process_requests<directed_compute_t, req_part_vertex_compute,
					index_comp_allocator_impl<req_part_vertex_compute> >(
						part_vertex_comps[type], req_part_vertex_comp_alloc,
						*index_reader);
			else
				process_requests<directed_compute_t, genrq_part_vertex_compute,
					general_index_comp_allocator_impl<genrq_part_vertex_compute> >(
						part_vertex_comps[type], genrq_part_vertex_comp_alloc,
						*index_reader);
		}
	}

	if (!edge_comps.empty()) {
		if (!std::is_sorted(edge_comps.begin(), edge_comps.end(),
					id_compute_less()))
			std::sort(edge_comps.begin(), edge_comps.end(), id_compute_less());
		process_requests<id_compute_t, req_edge_compute>(edge_comps,
				req_edge_comp_alloc, *index_reader);
	}

	if (!directed_edge_comps.empty()) {
		if (!std::is_sorted(directed_edge_comps.begin(), directed_edge_comps.end(),
					id_compute_less())) {
			std::sort(directed_edge_comps.begin(), directed_edge_comps.end(),
					id_compute_less());
			assert(std::is_sorted(directed_edge_comps.begin(), directed_edge_comps.end(),
					id_compute_less()));
		}
		if (is_dense(directed_edge_comps))
			process_requests<id_compute_t, req_directed_edge_compute,
				index_comp_allocator_impl<req_directed_edge_compute> >(
					directed_edge_comps, req_directed_edge_comp_alloc,
					*index_reader);
		else
			process_requests<id_compute_t, genrq_directed_edge_compute,
				general_index_comp_allocator_impl<genrq_directed_edge_compute> >(
					directed_edge_comps, genrq_directed_edge_comp_alloc,
					*index_reader);
	}
}
