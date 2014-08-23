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
#include "worker_thread.h"

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

#if 0
	// The starting offset of the entries in the cached index.
	size_t cached_index_start;
	// The values in the last one or two pages in the index.
	std::vector<ValueType> cached_index;
#endif

protected:
	ext_mem_vindex_reader_impl(io_interface::ptr io) {
		this->io = io;
		req_vertex_store = vertex_KV_store::create(io);

#if 0
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
#endif
	}

#if 0
	size_t get_cached_index_start() const {
		return cached_index_start;
	}

	const ValueType &get_cached_entry(int idx) const {
		return cached_index[idx];
	}
#endif
public:
	static ptr create(io_interface::ptr io) {
		return ptr(new ext_mem_vindex_reader_impl(io));
	}

	virtual void request_index(index_compute *compute) {
#if 0
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
#endif
			req_vertex_task<ValueType> task(compute->get_first_vertex(),
					compute->get_last_vertex() - compute->get_first_vertex() + 1,
					compute);
			req_vertex_store->async_request(task);
//		}
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
		id_range_t range = compute->get_range();
#if 0
		if (range.second < index->get_num_vertices()) {
#endif
			array_index_iterator_impl<ValueType> it(index->get_data() + range.first,
					// We need an additional entry.
					index->get_data() + range.second + 1);
			bool ret = compute->run(compute->get_first_vertex(), it);
			assert(ret);
#if 0
		}
		else {
			assert(range.second > range.first);
			if (range.second - range.first > 1) {
				array_index_iterator_impl<ValueType> it(
						index->get_data() + range.first,
						index->get_data() + range.second);
				compute->run(compute->get_first_vertex(), it);
			}

			// For the last vertex,
			ValueType vs[2];
			vs[0] = index->get_data()[index->get_num_vertices() - 1];
			vs[1] = ValueType(index->get_graph_size());
			array_index_iterator_impl<ValueType> it(vs, &vs[2]);
			bool ret = compute->run(index->get_num_vertices() - 1, it);
			assert(ret);
		}
#endif
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

template<class RequestType, class ComputeType, class AllocatorType, class SingleComputeType>
void process_requests(std::vector<RequestType> &reqs,
		AllocatorType &alloc, vertex_index_reader &index_reader,
		index_comp_allocator_impl<SingleComputeType> &single_alloc)
{
	assert(!reqs.empty());
	ComputeType *compute = (ComputeType *) alloc.alloc();
	const RequestType &req = reqs[0];
	compute->init(req.first, req.second);
	for (size_t i = 1; i < reqs.size(); i++) {
		const RequestType &req = reqs[i];

		if (!compute->add_vertex(req.first, req.second)) {
			if (compute->get_num_vertices() == 1) {
				SingleComputeType *single_compute
					= (SingleComputeType *) single_alloc.alloc();
				single_compute->init(*compute);
				compute->clear();
				index_reader.request_index(single_compute);
			}
			else {
				index_reader.request_index(compute);
				compute = (ComputeType *) alloc.alloc();
			}
			compute->init(req.first, req.second);
		}
	}
	assert(!compute->empty());
	if (compute->get_num_vertices() == 1) {
		SingleComputeType *single_compute
			= (SingleComputeType *) single_alloc.alloc();
		single_compute->init(*compute);
		compute->clear();
		index_reader.request_index(single_compute);
		alloc.free(compute);
	}
	else
		index_reader.request_index(compute);
	reqs.clear();
}

struct id_range_less
{
	bool operator()(const id_range_t &id1, const id_range_t &id2) {
		return id1.first < id2.first;
	}
};

void merge_vertex_requests(std::vector<id_range_t> &reqs)
{
	if (reqs.empty())
		return;

	// TODO I should test this.
	if (!std::is_sorted(reqs.begin(), reqs.end(), id_range_less()))
		std::sort(reqs.begin(), reqs.end(), id_range_less());
	id_range_t range = reqs.front();
	size_t write_back_idx = 0;
	for (size_t i = 1; i < reqs.size(); i++) {
		// Merge two ranges.
		if (range.second == reqs[i].first)
			range.second = reqs[i].second;
		else {
			// write the previous range to the vector.
			// and start a new range.
			assert(write_back_idx < i);
			reqs[write_back_idx++] = range;
			range = reqs[i];
		}
	}
	reqs[write_back_idx] = range;
	reqs.resize(write_back_idx + 1);
}

void simple_index_reader::process_self_requests(std::vector<id_range_t> &reqs,
		edge_type type)
{
	if (reqs.empty())
		return;
	merge_vertex_requests(reqs);
	// We always assume that we can merge multiple id ranges.
	sparse_self_vertex_compute *compute
		= (sparse_self_vertex_compute *) sparse_self_req_alloc->alloc();
	compute->init(reqs[0], t, type);

	for (size_t i = 1; i < reqs.size(); i++) {
		id_range_t range = reqs[i];
		// If we can merge, move to the next range.
		if (compute->add_range(range))
			continue;

		// If there are multiple ranges in the sparse index request.
		if (compute->get_num_ranges() > 1) {
			index_reader->request_index(compute);
			compute = (sparse_self_vertex_compute *) sparse_self_req_alloc->alloc();
			compute->init(range, t, type);
		}
		else {
			// If not, we prefer to query the index with the dense index
			// request.
			dense_self_vertex_compute *dense_compute
				= (dense_self_vertex_compute *) dense_self_req_alloc->alloc();
			dense_compute->init(compute->get_range(0), t, type);
			index_reader->request_index(dense_compute);
			compute->init(range, t, type);
		}
	}
	if (compute->get_num_ranges() > 1)
		index_reader->request_index(compute);
	else {
		dense_self_vertex_compute *dense_compute
			= (dense_self_vertex_compute *) dense_self_req_alloc->alloc();
		dense_compute->init(compute->get_range(0), t, type);
		index_reader->request_index(dense_compute);
		sparse_self_req_alloc->free(compute);
	}
	reqs.clear();
}

void simple_index_reader::init(worker_thread *t, bool directed)
{
	this->t = t;

	req_vertex_comp_alloc
		= new index_comp_allocator_impl<req_vertex_compute>(t);
	req_undirected_edge_comp_alloc
		= new index_comp_allocator_impl<req_undirected_edge_compute>(t);
	req_directed_edge_comp_alloc
		= new index_comp_allocator_impl<req_directed_edge_compute>(t);

	genrq_vertex_comp_alloc
		= new general_index_comp_allocator_impl<genrq_vertex_compute>(
				t, get_index_entry_size_log(directed));
	genrq_edge_comp_alloc
		= new general_index_comp_allocator_impl<genrq_edge_compute>(
				t, get_index_entry_size_log(directed));
	genrq_directed_edge_comp_alloc
		= new general_index_comp_allocator_impl<genrq_directed_edge_compute>(
				t, get_index_entry_size_log(directed));

	single_vertex_comp_alloc
		= new index_comp_allocator_impl<single_vertex_compute>(t);
	single_edge_comp_alloc
		= new index_comp_allocator_impl<single_edge_compute>(t);
	single_directed_edge_comp_alloc
		= new index_comp_allocator_impl<single_directed_edge_compute>(t);

	dense_self_req_alloc
		= new index_comp_allocator_impl<dense_self_vertex_compute>(t);
	sparse_self_req_alloc
		= new general_index_comp_allocator_impl<sparse_self_vertex_compute>(
				t, get_index_entry_size_log(directed));
}

void simple_index_reader::flush_computes()
{
	if (!vertex_comps.empty()) {
		if (!std::is_sorted(vertex_comps.begin(), vertex_comps.end(),
					id_compute_less()))
			std::sort(vertex_comps.begin(), vertex_comps.end(), id_compute_less());
		if (is_dense(vertex_comps))
			process_requests<id_compute_t, req_vertex_compute,
				index_comp_allocator_impl<req_vertex_compute>, single_vertex_compute>(
						vertex_comps, *req_vertex_comp_alloc, *index_reader,
						*single_vertex_comp_alloc);
		else
			process_requests<id_compute_t, genrq_vertex_compute,
				general_index_comp_allocator_impl<genrq_vertex_compute>,
				single_vertex_compute>(vertex_comps, *genrq_vertex_comp_alloc,
						*index_reader, *single_vertex_comp_alloc);
	}

	for (int type = edge_type::IN_EDGE; type < edge_type::NUM_TYPES; type++) {
		if (!part_vertex_comps[type].empty()) {
			if (!std::is_sorted(part_vertex_comps[type].begin(),
						part_vertex_comps[type].end(), directed_compute_less()))
				std::sort(part_vertex_comps[type].begin(),
						part_vertex_comps[type].end(), directed_compute_less());
			if (is_dense(part_vertex_comps[type]))
				process_requests<directed_compute_t, req_vertex_compute,
					index_comp_allocator_impl<req_vertex_compute>,
					single_vertex_compute>(part_vertex_comps[type],
							*req_vertex_comp_alloc, *index_reader,
							*single_vertex_comp_alloc);
			else
				process_requests<directed_compute_t, genrq_vertex_compute,
					general_index_comp_allocator_impl<genrq_vertex_compute>,
					single_vertex_compute>(part_vertex_comps[type],
							*genrq_vertex_comp_alloc, *index_reader,
							*single_vertex_comp_alloc);
		}
	}

	if (!edge_comps.empty()) {
		if (!std::is_sorted(edge_comps.begin(), edge_comps.end(),
					id_compute_less()))
			std::sort(edge_comps.begin(), edge_comps.end(), id_compute_less());
		process_requests<id_compute_t, req_undirected_edge_compute,
			index_comp_allocator_impl<req_undirected_edge_compute>, single_edge_compute>(
				edge_comps, *req_undirected_edge_comp_alloc, *index_reader,
				*single_edge_comp_alloc);
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
				index_comp_allocator_impl<req_directed_edge_compute>,
				single_directed_edge_compute>(directed_edge_comps,
						*req_directed_edge_comp_alloc, *index_reader,
						*single_directed_edge_comp_alloc);
		else
			process_requests<id_compute_t, genrq_directed_edge_compute,
				general_index_comp_allocator_impl<genrq_directed_edge_compute>,
				single_directed_edge_compute>(directed_edge_comps,
						*genrq_directed_edge_comp_alloc, *index_reader,
						*single_directed_edge_comp_alloc);
	}

	// This is to request undirected vertices.
	if (!self_undirected_reqs.empty())
		process_self_requests(self_undirected_reqs, IN_EDGE);

	// This is to request directed vertices.
	for (int i = 0; i < edge_type::NUM_TYPES; i++) {
		if (!self_part_reqs[i].empty()) {
			process_self_requests(self_part_reqs[i], (edge_type) i);
		}
	}
}

bool dense_self_vertex_compute::run(vertex_id_t start_vid, index_iterator &it)
{
	assert(start_vid == get_first_vertex());
	assert(get_num_vertices() == it.get_num_entries() - 1);

	merged_vertex_compute *compute
		= (merged_vertex_compute *) thread->get_merged_compute_allocator().alloc();
	compute->init(start_vid, get_num_vertices(), type);
	// For undirected vertices, the edge type is always IN_EDGE.
	// For directed vertices, if the edge type is BOTH_EDGES, we issue
	// two requests with one vertex compute.
	if (type == edge_type::IN_EDGE) {
		off_t first_off = it.get_curr_off();
		assert(it.move_to(get_num_vertices() - 1));
		off_t last_off = it.get_curr_off() + it.get_curr_size();
		data_loc_t loc(this->thread->get_graph().get_file_id(), first_off);
		io_request req(compute, loc, last_off - first_off, READ);
		this->thread->issue_io_request(req);
	}
	else if (type == edge_type::OUT_EDGE) {
		off_t first_off = it.get_curr_out_off();
		assert(it.move_to(get_num_vertices() - 1));
		off_t last_off = it.get_curr_out_off() + it.get_curr_out_size();
		data_loc_t loc(this->thread->get_graph().get_file_id(), first_off);
		io_request req(compute, loc, last_off - first_off, READ);
		this->thread->issue_io_request(req);
	}
	else {
		assert(type == edge_type::BOTH_EDGES);
		off_t first_in_off = it.get_curr_off();
		off_t first_out_off = it.get_curr_out_off();
		assert(it.move_to(get_num_vertices() - 1));
		off_t last_in_off = it.get_curr_off() + it.get_curr_size();
		off_t last_out_off = it.get_curr_out_off() + it.get_curr_out_size();

		data_loc_t in_loc(this->thread->get_graph().get_file_id(), first_in_off);
		io_request in_req(compute, in_loc, last_in_off - first_in_off, READ);
		this->thread->issue_io_request(in_req);

		data_loc_t out_loc(this->thread->get_graph().get_file_id(), first_out_off);
		io_request out_req(compute, out_loc, last_out_off - first_out_off, READ);
		this->thread->issue_io_request(out_req);
	}

	return true;
}

bool sparse_self_vertex_compute::run(vertex_id_t start_vid, index_iterator &it)
{
	assert(start_vid == get_first_vertex());
	assert(get_last_vertex() - get_first_vertex() + 1
			== it.get_num_entries() - 1);
	// If #ragnes is 1, we'll do it in the dense_self_vertex_compute.
	assert(num_ranges > 1);

	for (int i = 0; i < num_ranges; i++) {
		merged_vertex_compute *compute
			= (merged_vertex_compute *) thread->get_merged_compute_allocator().alloc();
		id_range_t range = ranges[i];
		vsize_t num_vertices = range.second - range.first;
		compute->init(range.first, num_vertices, type);
		off_t idx_entry_loc = range.first - start_vid;
		assert(it.move_to(idx_entry_loc));
		// For undirected vertices, the edge type is always IN_EDGE.
		// For directed vertices, if the edge type is BOTH_EDGES, we issue
		// two requests with one vertex compute.
		if (type == edge_type::IN_EDGE) {
			off_t first_off = it.get_curr_off();
			assert(it.move_to(idx_entry_loc + num_vertices - 1));
			off_t last_off = it.get_curr_off() + it.get_curr_size();
			data_loc_t loc(this->thread->get_graph().get_file_id(), first_off);
			io_request req(compute, loc, last_off - first_off, READ);
			this->thread->issue_io_request(req);
		}
		else if (type == edge_type::OUT_EDGE) {
			off_t first_off = it.get_curr_out_off();
			assert(it.move_to(idx_entry_loc + num_vertices - 1));
			off_t last_off = it.get_curr_out_off() + it.get_curr_out_size();
			data_loc_t loc(this->thread->get_graph().get_file_id(), first_off);
			io_request req(compute, loc, last_off - first_off, READ);
			this->thread->issue_io_request(req);
		}
		else {
			assert(type == edge_type::BOTH_EDGES);
			off_t first_in_off = it.get_curr_off();
			off_t first_out_off = it.get_curr_out_off();
			assert(it.move_to(idx_entry_loc + num_vertices - 1));
			off_t last_in_off = it.get_curr_off() + it.get_curr_size();
			off_t last_out_off = it.get_curr_out_off() + it.get_curr_out_size();

			data_loc_t in_loc(this->thread->get_graph().get_file_id(), first_in_off);
			io_request in_req(compute, in_loc, last_in_off - first_in_off, READ);
			this->thread->issue_io_request(in_req);

			data_loc_t out_loc(this->thread->get_graph().get_file_id(), first_out_off);
			io_request out_req(compute, out_loc, last_out_off - first_out_off, READ);
			this->thread->issue_io_request(out_req);
		}
	}

	return true;
}
