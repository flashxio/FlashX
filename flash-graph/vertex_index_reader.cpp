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

using namespace safs;

namespace fg
{

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
protected:
	ext_mem_vindex_reader_impl(io_interface::ptr io) {
		this->io = io;
		req_vertex_store = vertex_KV_store::create(io);
	}
public:
	static ptr create(io_interface::ptr io) {
		return ptr(new ext_mem_vindex_reader_impl(io));
	}

	virtual void request_index(index_compute *compute) {
		req_vertex_task<ValueType> task(compute->get_first_vertex(),
				compute->get_last_vertex() - compute->get_first_vertex() + 1,
				compute);
		req_vertex_store->async_request(task);
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

	in_mem_vindex_reader_impl(vertex_index::ptr index) {
		this->index = vertex_index_temp<ValueType>::cast(index);
	}
public:
	static ptr create(vertex_index::ptr index) {
		return ptr(new in_mem_vindex_reader_impl<ValueType>(index));
	}

	virtual void request_index(index_compute *compute) {
		id_range_t range = compute->get_range();
		array_index_iterator_impl<ValueType> it(index->get_data() + range.first,
				// We need an additional entry.
				index->get_data() + range.second + 1);
		BOOST_VERIFY(compute->run(compute->get_first_vertex(), it));
		compute->get_allocator().free(compute);
	}

	virtual void wait4complete(int num) {
	}

	virtual size_t get_num_pending_tasks() const {
		return 0;
	}
};

template<class vertex_index_type, class iterator_type>
class in_mem_cindex_reader: public vertex_index_reader
{
	typename vertex_index_type::ptr index;

protected:
	in_mem_cindex_reader(typename vertex_index_type::ptr index) {
		this->index = index;
	}
public:
	static ptr create(typename vertex_index_type::ptr index) {
		return ptr(new in_mem_cindex_reader<vertex_index_type,
				iterator_type>(index));
	}


	virtual void request_index(index_compute *compute) {
		id_range_t range = compute->get_range();
		iterator_type it(*index, range);
		BOOST_VERIFY(compute->run(compute->get_first_vertex(), it));
		compute->get_allocator().free(compute);
	}

	virtual void wait4complete(int num) {
	}

	virtual size_t get_num_pending_tasks() const {
		return 0;
	}
};

vertex_index_reader::ptr vertex_index_reader::create(
		const in_mem_query_vertex_index::ptr index, bool directed)
{
	bool compressed = index->is_compressed();
	if (!compressed && directed)
		return in_mem_vindex_reader_impl<directed_vertex_entry>::create(
				index->get_raw_index());
	else if (!compressed && !directed)
		return in_mem_vindex_reader_impl<vertex_offset>::create(
				index->get_raw_index());
	else if (compressed && directed)
		return in_mem_cindex_reader<in_mem_cdirected_vertex_index,
			   compressed_directed_index_iterator>::create(
					   in_mem_cdirected_vertex_index::cast(index));
	else
		return in_mem_cindex_reader<in_mem_cundirected_vertex_index,
			   compressed_undirected_index_iterator>::create(
					   in_mem_cundirected_vertex_index::cast(index));
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

	int entry_size_log = get_index_entry_size_log(directed);
	genrq_vertex_comp_alloc
		= new general_index_comp_allocator_impl<genrq_vertex_compute>(t,
				entry_size_log, in_mem);
	genrq_edge_comp_alloc
		= new general_index_comp_allocator_impl<genrq_edge_compute>(t,
				entry_size_log, in_mem);
	genrq_directed_edge_comp_alloc
		= new general_index_comp_allocator_impl<genrq_directed_edge_compute>(t,
				entry_size_log, in_mem);

	single_vertex_comp_alloc
		= new index_comp_allocator_impl<single_vertex_compute>(t);
	single_edge_comp_alloc
		= new index_comp_allocator_impl<single_edge_compute>(t);
	single_directed_edge_comp_alloc
		= new index_comp_allocator_impl<single_directed_edge_compute>(t);

	dense_self_req_alloc
		= new index_comp_allocator_impl<dense_self_vertex_compute>(t);
	sparse_self_req_alloc
		= new general_index_comp_allocator_impl<sparse_self_vertex_compute>(t,
				entry_size_log, in_mem);
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
	assert(get_num_vertices() == it.get_num_vertices());

	merged_vertex_compute *compute
		= (merged_vertex_compute *) thread->get_merged_compute_allocator().alloc();
	compute->init(start_vid, get_num_vertices(), type);
	// For undirected vertices, the edge type is always IN_EDGE.
	// For directed vertices, if the edge type is BOTH_EDGES, we issue
	// two requests with one vertex compute.
	if (type == edge_type::IN_EDGE) {
		off_t first_off = it.get_curr_off();
		BOOST_VERIFY(it.move_to(get_num_vertices() - 1));
		off_t last_off = it.get_curr_off() + it.get_curr_size();
		data_loc_t loc(this->thread->get_graph().get_file_id(), first_off);
		io_request req(compute, loc, last_off - first_off, READ);
		this->thread->issue_io_request(req);
	}
	else if (type == edge_type::OUT_EDGE) {
		off_t first_off = it.get_curr_out_off();
		BOOST_VERIFY(it.move_to(get_num_vertices() - 1));
		off_t last_off = it.get_curr_out_off() + it.get_curr_out_size();
		data_loc_t loc(this->thread->get_graph().get_file_id(), first_off);
		io_request req(compute, loc, last_off - first_off, READ);
		this->thread->issue_io_request(req);
	}
	else {
		assert(type == edge_type::BOTH_EDGES);
		off_t first_in_off = it.get_curr_off();
		off_t first_out_off = it.get_curr_out_off();
		BOOST_VERIFY(it.move_to(get_num_vertices() - 1));
		off_t last_in_off = it.get_curr_off() + it.get_curr_size();
		off_t last_out_off = it.get_curr_out_off() + it.get_curr_out_size();

		data_loc_t in_loc(this->thread->get_graph().get_file_id(), first_in_off);
		io_request in_req(compute, in_loc, last_in_off - first_in_off, READ);
		this->thread->issue_io_request(in_req);

		// issue_io_request doesn't really issue the I/O request to
		// the underlying I/O yet, so we don't need to worry that user
		// compute is free'd by the I/O.
		data_loc_t out_loc(this->thread->get_graph().get_file_id(), first_out_off);
		io_request out_req(compute, out_loc, last_out_off - first_out_off, READ);
		this->thread->issue_io_request(out_req);
	}

	return true;
}

typedef std::pair<off_t, off_t> off_range_t;

static sparse_vertex_compute *issue_request(worker_thread *thread,
		const off_range_t &off_range, sparse_vertex_compute *compute,
		edge_type type, bool return_compute)
{
	if (compute->get_num_ranges() == 1) {
		merged_vertex_compute *dense_compute
			= (merged_vertex_compute *) thread->get_merged_compute_allocator().alloc();
		dense_compute->init(compute->get_first_vertex(),
				compute->get_num_vertices(), type);
		data_loc_t loc(thread->get_graph().get_file_id(), off_range.first);
		io_request req(dense_compute, loc, off_range.second - off_range.first,
				READ);
		thread->issue_io_request(req);
		if (return_compute)
			return compute;
		else {
			thread->get_sparse_compute_allocator().free(compute);
			return NULL;
		}
	}
	else {
		data_loc_t loc(thread->get_graph().get_file_id(), off_range.first);
		io_request req(compute, loc, off_range.second - off_range.first, READ);
		thread->issue_io_request(req);
		if (return_compute)
			return (sparse_vertex_compute *) thread->get_sparse_compute_allocator().alloc();
		else
			return NULL;
	}
}

static sparse_vertex_compute *issue_request(worker_thread *thread,
		const off_range_t off_ranges[], sparse_vertex_compute *compute,
		bool return_compute)
{
	if (compute->get_num_ranges() == 1) {
		merged_vertex_compute *dense_compute
			= (merged_vertex_compute *) thread->get_merged_compute_allocator().alloc();
		dense_compute->init(compute->get_first_vertex(),
				compute->get_num_vertices(), BOTH_EDGES);
		data_loc_t loc(thread->get_graph().get_file_id(), off_ranges[0].first);
		io_request req(dense_compute, loc, off_ranges[0].second - off_ranges[0].first,
				READ);
		thread->issue_io_request(req);

		loc = data_loc_t(thread->get_graph().get_file_id(), off_ranges[1].first);
		req = io_request(dense_compute, loc, off_ranges[1].second - off_ranges[1].first,
				READ);
		thread->issue_io_request(req);
		if (return_compute)
			return compute;
		else {
			thread->get_sparse_compute_allocator().free(compute);
			return NULL;
		}
	}
	else {
		data_loc_t loc(thread->get_graph().get_file_id(), off_ranges[0].first);
		io_request req(compute, loc, off_ranges[0].second - off_ranges[0].first, READ);
		thread->issue_io_request(req);

		loc = data_loc_t(thread->get_graph().get_file_id(), off_ranges[1].first);
		req = io_request(compute, loc, off_ranges[1].second - off_ranges[1].first, READ);
		thread->issue_io_request(req);
		if (return_compute)
			return (sparse_vertex_compute *) thread->get_sparse_compute_allocator().alloc();
		else
			return NULL;
	}
}

static off_range_t get_in_off_range(index_iterator &it, vertex_id_t start_vid,
		const id_range_t &range)
{
	vsize_t num_vertices = range.second - range.first;
	off_t idx_entry_loc = range.first - start_vid;
	BOOST_VERIFY(it.move_to(idx_entry_loc));

	off_t first_off = it.get_curr_off();
	BOOST_VERIFY(it.move_to(idx_entry_loc + num_vertices - 1));
	off_t last_off = it.get_curr_off() + it.get_curr_size();
	return off_range_t(first_off, last_off);
}

static off_range_t get_out_off_range(index_iterator &it, vertex_id_t start_vid,
		const id_range_t &range)
{
	vsize_t num_vertices = range.second - range.first;
	off_t idx_entry_loc = range.first - start_vid;
	BOOST_VERIFY(it.move_to(idx_entry_loc));

	off_t first_off = it.get_curr_out_off();
	BOOST_VERIFY(it.move_to(idx_entry_loc + num_vertices - 1));
	off_t last_off = it.get_curr_out_off() + it.get_curr_out_size();
	return off_range_t(first_off, last_off);
}

static bool can_merge_reqs(const off_range_t &range1, const off_range_t &range2)
{
	return ROUND_PAGE(range1.second) == ROUND_PAGE(range2.first)
		|| ROUND_PAGE(range1.second) + PAGE_SIZE == ROUND_PAGE(range2.first);
}

static void merge_reqs(off_range_t &range1, const off_range_t &range2)
{
	range1.second = range2.second;
}

void sparse_self_vertex_compute::run_in_vertices(vertex_id_t start_vid,
		index_iterator &it)
{
	sparse_vertex_compute *compute
		= (sparse_vertex_compute *) thread->get_sparse_compute_allocator().alloc();
	off_range_t off_range = get_in_off_range(it, start_vid, ranges[0]);
	compute->init(ranges[0], &off_range, IN_EDGE);
	off_range_t req_range = off_range;
	for (size_t i = 1; i < num_ranges; i++) {
		off_range = get_in_off_range(it, start_vid, ranges[i]);
		if (can_merge_reqs(req_range, off_range)
				&& compute->add_range(ranges[i], &off_range)) {
			merge_reqs(req_range, off_range);
		}
		else {
			compute = issue_request(thread, req_range, compute, IN_EDGE, true);
			compute->init(ranges[i], &off_range, IN_EDGE);
			req_range = off_range;
		}
	}
	issue_request(thread, req_range, compute, IN_EDGE, false);
}

void sparse_self_vertex_compute::run_out_vertices(vertex_id_t start_vid,
		index_iterator &it)
{
	sparse_vertex_compute *compute
		= (sparse_vertex_compute *) thread->get_sparse_compute_allocator().alloc();
	off_range_t off_range = get_out_off_range(it, start_vid, ranges[0]);
	compute->init(ranges[0], &off_range, OUT_EDGE);
	off_range_t req_range = off_range;
	for (size_t i = 1; i < num_ranges; i++) {
		off_range = get_out_off_range(it, start_vid, ranges[i]);
		if (can_merge_reqs(req_range, off_range)
				&& compute->add_range(ranges[i], &off_range))
			merge_reqs(req_range, off_range);
		else {
			compute = issue_request(thread, req_range, compute, OUT_EDGE, true);
			compute->init(ranges[i], &off_range, OUT_EDGE);
			req_range = off_range;
		}
	}
	issue_request(thread, req_range, compute, OUT_EDGE, false);
}

void sparse_self_vertex_compute::run_both_vertices(vertex_id_t start_vid,
		index_iterator &it)
{
	sparse_vertex_compute *compute
		= (sparse_vertex_compute *) thread->get_sparse_compute_allocator().alloc();
	off_range_t off_ranges[2];
	off_ranges[0] = get_in_off_range(it, start_vid, ranges[0]);
	off_ranges[1] = get_out_off_range(it, start_vid, ranges[0]);
	compute->init(ranges[0], off_ranges, BOTH_EDGES);
	off_range_t req_ranges[2];
	req_ranges[0] = off_ranges[0];
	req_ranges[1] = off_ranges[1];
	for (size_t i = 1; i < num_ranges; i++) {
		off_ranges[0] = get_in_off_range(it, start_vid, ranges[i]);
		off_ranges[1] = get_out_off_range(it, start_vid, ranges[i]);
		if (can_merge_reqs(req_ranges[0], off_ranges[0])
				&& can_merge_reqs(req_ranges[1], off_ranges[1])
				&& compute->add_range(ranges[i], off_ranges)) {
			merge_reqs(req_ranges[0], off_ranges[0]);
			merge_reqs(req_ranges[1], off_ranges[1]);
		}
		else {
			compute = issue_request(thread, req_ranges, compute, true);
			compute->init(ranges[i], off_ranges, BOTH_EDGES);
			req_ranges[0] = off_ranges[0];
			req_ranges[1] = off_ranges[1];
		}
	}
	issue_request(thread, req_ranges, compute, false);
}

bool sparse_self_vertex_compute::run(vertex_id_t start_vid, index_iterator &it)
{
	assert(start_vid == get_first_vertex());
	assert(get_last_vertex() - get_first_vertex() + 1
			== it.get_num_vertices());
	// If #ragnes is 1, we'll do it in the dense_self_vertex_compute.
	assert(num_ranges > 1);

	switch(type) {
		case IN_EDGE:
			run_in_vertices(start_vid, it);
			break;
		case OUT_EDGE:
			run_out_vertices(start_vid, it);
			break;
		case BOTH_EDGES:
			run_both_vertices(start_vid, it);
			break;
		default:
			ABORT_MSG("wrong edge type");
	}

	return true;
}

}
