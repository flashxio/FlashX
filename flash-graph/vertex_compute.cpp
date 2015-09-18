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
#include "vertex_compute.h"
#include "graph_engine.h"
#include "worker_thread.h"
#include "vertex_index_reader.h"

request_range vertex_compute::get_next_request()
{
	// Get the next vertex.
	const ext_mem_vertex_info info = requested_vertices.top();
	requested_vertices.pop();
	data_loc_t loc(graph->get_file_id(), info.get_off());
	num_issued++;
	return request_range(loc, info.get_size(), READ, this);
}

void vertex_compute::start_run()
{
	issue_thread->start_run_vertex(v);
}

void vertex_compute::finish_run()
{
	bool issued_reqs = issue_thread->finish_run_vertex(v);
	// If we have completed all pending requests and this run didn't
	// issue more requests, we can be sure that the vertex has completed.
	// We need to notify the thread that initiate processing the vertex
	// of the completion of the vertex.
	if (get_num_pending() == 0) {
		// this is to double check. If there are pending requests,
		// the vertex shouldn't have issued requests in this run.
		BOOST_VERIFY(!issued_reqs);
		issue_thread->complete_vertex(v);
	}
}

void vertex_compute::request_vertices(vertex_id_t ids[], size_t num)
{
	num_requested += num;
	issue_thread->get_index_reader().request_vertices(ids, num, *this);
}

void vertex_compute::request_num_edges(vertex_id_t ids[], size_t num)
{
	num_edge_requests += num;
	issue_thread->get_index_reader().request_num_edges(ids, num, *this);
}

void vertex_compute::run_on_vertex_size(vertex_id_t id, vsize_t size)
{
	start_run();
	vsize_t num_edges = issue_thread->get_graph().cal_num_edges(size);
	vertex_header header(id, num_edges);
	issue_thread->get_vertex_program(v.is_part()).run_on_num_edges(*v, header);
	num_edge_completed++;
	finish_run();
}

void vertex_compute::issue_io_request(const ext_mem_vertex_info &info)
{
	// If the vertex compute has been issued to SAFS, SAFS will get the IO
	// request from the interface of user_compute. In this case, we only
	// need to add the I/O request to the queue.
	if (issued_to_io()) {
		requested_vertices.push(info);
	}
	else {
		// Otherwise, we need to issue the I/O request to SAFS explicitly.
		data_loc_t loc(graph->get_file_id(), info.get_off());
		io_request req(this, loc, info.get_size(), READ);
		num_issued++;
		issue_thread->issue_io_request(req);
	}
}

void vertex_compute::run(page_byte_array &array)
{
	num_complete_fetched++;
	start_run();
	page_undirected_vertex pg_v(array);
	issue_thread->get_vertex_program(v.is_part()).run(*v, pg_v);
	finish_run();
}

void directed_vertex_compute::run_on_page_vertex(page_directed_vertex &pg_v)
{
	start_run();
	issue_thread->get_vertex_program(v.is_part()).run(*v, pg_v);
	finish_run();
}

void directed_vertex_compute::run(page_byte_array &array)
{
	num_complete_fetched++;
	// If the combine map is empty, we don't need to merge
	// byte arrays.
	if (combine_map.empty()) {
		page_directed_vertex pg_v(array,
				(size_t) array.get_offset() < graph->get_in_part_size());
		run_on_page_vertex(pg_v);
		return;
	}

	vertex_id_t id = page_directed_vertex::get_id(array);
	combine_map_t::iterator it = combine_map.find(id);
	// If the vertex isn't in the combine map, we don't need to
	// merge byte arrays.
	if (it == combine_map.end()) {
		page_directed_vertex pg_v(array,
				(size_t) array.get_offset() < graph->get_in_part_size());
		run_on_page_vertex(pg_v);
		return;
	}
	else if (it->second == NULL) {
		page_byte_array *arr_copy = array.clone();
		assert(arr_copy);
		it->second = arr_copy;
	}
	else {
		page_byte_array *in_arr;
		page_byte_array *out_arr;
		if ((size_t) it->second->get_offset() < get_graph().get_in_part_size()) {
			in_arr = it->second;
			out_arr = &array;
			assert((size_t) array.get_offset() >= get_graph().get_in_part_size());
		}
		else {
			out_arr = it->second;
			in_arr = &array;
			assert((size_t) array.get_offset() < get_graph().get_in_part_size());
		}
		page_directed_vertex pg_v(*in_arr, *out_arr);
		run_on_page_vertex(pg_v);
		page_byte_array::destroy(it->second);
		combine_map.erase(it);
	}
}

void directed_vertex_compute::request_vertices(vertex_id_t ids[], size_t num)
{
	stack_array<directed_vertex_request> reqs(num);
	for (size_t i = 0; i < num; i++)
		reqs[i] = directed_vertex_request(ids[i], edge_type::BOTH_EDGES);
	request_partial_vertices(reqs.data(), num);
}

void directed_vertex_compute::request_partial_vertices(
		directed_vertex_request reqs[], size_t num)
{
	for (size_t i = 0; i < num; i++) {
		if (reqs[i].get_type() == edge_type::BOTH_EDGES)
			num_requested += 2;
		else
			num_requested++;
	}
	issue_thread->get_index_reader().request_vertices(reqs, num, *this);
}

void directed_vertex_compute::run_on_vertex_size(vertex_id_t id,
		size_t in_size, size_t out_size)
{
	start_run();
	vsize_t num_in_edges = issue_thread->get_graph().cal_num_edges(in_size);
	vsize_t num_out_edges = issue_thread->get_graph().cal_num_edges(out_size);
	directed_vertex_header header(id, num_in_edges, num_out_edges);
	issue_thread->get_vertex_program(v.is_part()).run_on_num_edges(*v, header);
	num_edge_completed++;
	finish_run();
}

void directed_vertex_compute::issue_io_request(const ext_mem_vertex_info &in_info,
		const ext_mem_vertex_info &out_info)
{
	assert(in_info.get_id() == out_info.get_id());
	if (issued_to_io()) {
		requested_vertices.push(in_info);
		requested_vertices.push(out_info);
	}
	else {
		// Otherwise, we need to issue the I/O request to SAFS explicitly.
		data_loc_t loc1(graph->get_file_id(), in_info.get_off());
		io_request req1(this, loc1, in_info.get_size(), READ);
		issue_thread->issue_io_request(req1);

		data_loc_t loc2(graph->get_file_id(), out_info.get_off());
		io_request req2(this, loc2, out_info.get_size(), READ);
		issue_thread->issue_io_request(req2);
		num_issued += 2;
	}

	combine_map.insert(combine_map_t::value_type(in_info.get_id(), NULL));
}

void directed_vertex_compute::request_num_edges(vertex_id_t ids[], size_t num)
{
	num_edge_requests += num;
	issue_thread->get_index_reader().request_num_directed_edges(ids,
			num, *this);
}

void merged_vertex_compute::start_run(compute_vertex_pointer v)
{
	issue_thread->start_run_vertex(v);
}

void merged_vertex_compute::finish_run(compute_vertex_pointer v)
{
	bool issued_reqs = issue_thread->finish_run_vertex(v);
	// TODO we have to make sure that this vertex didn't issue another vertex
	// request.
	// The vertex only issued one request, which is just processed.
	// If this run didn't issue more requests, we can be sure that
	// the vertex has completed in this iteration.
	// We need to notify the thread that initiate processing the vertex
	// of the completion of the vertex.
	if (!issued_reqs)
		issue_thread->complete_vertex(v);
}

void merged_undirected_vertex_compute::run(page_byte_array &array)
{
	off_t off = 0;
	vertex_id_t id = this->get_start_id();
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	// We don't support part vertex compute here.
	vertex_program &curr_vprog = t->get_vertex_program(false);
	for (int i = 0; i < get_num_vertices(); i++, id++) {
		sub_page_byte_array sub_arr(array, off);
		page_undirected_vertex pg_v(sub_arr);
		assert(pg_v.get_id() == id);
		compute_vertex_pointer v(&get_graph().get_vertex(pg_v.get_id()));
		start_run(v);
		curr_vprog.run(*v, pg_v);
		finish_run(v);
		off += pg_v.get_size();
	}

	complete = true;
}

void merged_directed_vertex_compute::run_on_array(page_byte_array &array)
{
	off_t off = 0;
	vertex_id_t id = this->get_start_id();
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	// We don't support part vertex compute here.
	vertex_program &curr_vprog = t->get_vertex_program(false);
	bool in_part = (size_t) array.get_offset() < get_graph().get_in_part_size();
	for (int i = 0; i < get_num_vertices(); i++, id++) {
		sub_page_byte_array sub_arr(array, off);
		page_directed_vertex pg_v(sub_arr, in_part);
		assert(pg_v.get_id() == id);
		compute_vertex_pointer v(&get_graph().get_vertex(pg_v.get_id()));
		start_run(v);
		curr_vprog.run(*v, pg_v);
		finish_run(v);
		if (in_part)
			off += pg_v.get_in_size();
		else
			off += pg_v.get_out_size();
	}
}

void merged_directed_vertex_compute::run_on_arrays(page_byte_array &in_arr,
		page_byte_array &out_arr)
{
	off_t in_off = 0;
	off_t out_off = 0;
	vertex_id_t id = this->get_start_id();
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	// We don't support part vertex compute here.
	vertex_program &curr_vprog = t->get_vertex_program(false);
	for (int i = 0; i < get_num_vertices(); i++, id++) {
		sub_page_byte_array sub_in_arr(in_arr, in_off);
		sub_page_byte_array sub_out_arr(out_arr, out_off);
		page_directed_vertex pg_v(sub_in_arr, sub_out_arr);
		assert(pg_v.get_id() == id);
		compute_vertex_pointer v(&get_graph().get_vertex(pg_v.get_id()));
		start_run(v);
		curr_vprog.run(*v, pg_v);
		finish_run(v);
		in_off += pg_v.get_in_size();
		out_off += pg_v.get_out_size();
	}
}

void merged_directed_vertex_compute::run(page_byte_array &arr)
{
	this->num_fetched_arrs++;
	assert(num_fetched_arrs <= num_required_arrs);
	if (type == BOTH_EDGES && buffered_arr) {
		page_byte_array *in_arr;
		page_byte_array *out_arr;

		if ((size_t) buffered_arr->get_offset() < get_graph().get_in_part_size()) {
			in_arr = buffered_arr;
			assert((size_t) arr.get_offset() >= get_graph().get_in_part_size());
			out_arr = &arr;
		}
		else {
			out_arr = buffered_arr;
			assert((size_t) arr.get_offset() < get_graph().get_in_part_size());
			in_arr = &arr;
		}

		run_on_arrays(*in_arr, *out_arr);
		page_byte_array::destroy(buffered_arr);
	}
	else if (type == BOTH_EDGES) {
		buffered_arr = arr.clone();
	}
	else if (type == IN_EDGE) {
		assert((size_t) arr.get_offset() < get_graph().get_in_part_size());
		run_on_array(arr);
	}
	else if (type == OUT_EDGE) {
		assert((size_t) arr.get_offset() >= get_graph().get_in_part_size());
		run_on_array(arr);
	}
	else
		ABORT_MSG("wrong type");
}

void sparse_vertex_compute::start_run(compute_vertex_pointer v)
{
	issue_thread->start_run_vertex(v);
}

void sparse_vertex_compute::finish_run(compute_vertex_pointer v)
{
	bool issued_reqs = issue_thread->finish_run_vertex(v);
	// TODO we have to make sure that this vertex didn't issue another vertex
	// request.
	// The vertex only issued one request, which is just processed.
	// If this run didn't issue more requests, we can be sure that
	// the vertex has completed in this iteration.
	// We need to notify the thread that initiate processing the vertex
	// of the completion of the vertex.
	if (!issued_reqs)
		issue_thread->complete_vertex(v);
}

void sparse_undirected_vertex_compute::run(page_byte_array &arr)
{
	assert(arr.get_offset() + arr.get_size() > (size_t) ranges[num_ranges - 1].start_off);
	vertex_program &curr_vprog = issue_thread->get_vertex_program(false);
	for (int i = 0; i < num_ranges; i++) {
		vertex_id_t id = this->ranges[i].id_range.first;
		int num_vertices
			= this->ranges[i].id_range.second - this->ranges[i].id_range.first;
		off_t off = this->ranges[i].start_off - arr.get_offset();
		for (int j = 0; j < num_vertices; j++, id++) {
			sub_page_byte_array sub_arr(arr, off);
			page_undirected_vertex pg_v(sub_arr);
			assert(pg_v.get_id() == id);
			compute_vertex_pointer v(&get_graph().get_vertex(pg_v.get_id()));
			start_run(v);
			curr_vprog.run(*v, pg_v);
			finish_run(v);
			off += pg_v.get_size();
		}
	}
	complete = true;
}

void sparse_directed_vertex_compute::run_on_array(page_byte_array &arr)
{
	assert(arr.get_offset() + arr.get_size() > (size_t) ranges[num_ranges - 1].start_off);
	vertex_program &curr_vprog = issue_thread->get_vertex_program(false);
	for (int i = 0; i < num_ranges; i++) {
		vertex_id_t id = this->ranges[i].id_range.first;
		int num_vertices
			= this->ranges[i].id_range.second - this->ranges[i].id_range.first;
		off_t off = this->ranges[i].start_off - arr.get_offset();
		// We don't support part vertex compute here.
		bool in_part = (size_t) arr.get_offset() < get_graph().get_in_part_size();
		for (int j = 0; j < num_vertices; j++, id++) {
			sub_page_byte_array sub_arr(arr, off);
			page_directed_vertex pg_v(sub_arr, in_part);
			assert(pg_v.get_id() == id);
			compute_vertex_pointer v(&get_graph().get_vertex(pg_v.get_id()));
			start_run(v);
			curr_vprog.run(*v, pg_v);
			finish_run(v);
			if (in_part)
				off += pg_v.get_in_size();
			else
				off += pg_v.get_out_size();
		}
	}
	complete = true;
}

void sparse_directed_vertex_compute::run_on_arrays(page_byte_array &in_arr,
		page_byte_array &out_arr)
{
	assert(in_arr.get_offset()
			+ in_arr.get_size() > (size_t) ranges[num_ranges - 1].start_off);
	assert(out_arr.get_offset()
			+ out_arr.get_size() > (size_t) out_start_offs[num_ranges - 1]);
	assert((size_t) num_ranges == out_start_offs.size());
	// We don't support part vertex compute here.
	vertex_program &curr_vprog = issue_thread->get_vertex_program(false);

	for (int i = 0; i < num_ranges; i++) {
		vertex_id_t id = this->ranges[i].id_range.first;
		int num_vertices
			= this->ranges[i].id_range.second - this->ranges[i].id_range.first;
		off_t in_off = this->ranges[i].start_off - in_arr.get_offset();
		off_t out_off = this->out_start_offs[i] - out_arr.get_offset();
		for (int i = 0; i < num_vertices; i++, id++) {
			sub_page_byte_array sub_in_arr(in_arr, in_off);
			sub_page_byte_array sub_out_arr(out_arr, out_off);
			page_directed_vertex pg_v(sub_in_arr, sub_out_arr);
			assert(pg_v.get_id() == id);
			compute_vertex_pointer v(&get_graph().get_vertex(pg_v.get_id()));
			start_run(v);
			curr_vprog.run(*v, pg_v);
			finish_run(v);
			in_off += pg_v.get_in_size();
			out_off += pg_v.get_out_size();
		}
	}
	complete = true;
}

void sparse_directed_vertex_compute::run(page_byte_array &arr)
{
	if (type == BOTH_EDGES && buffered_arr) {
		page_byte_array *in_arr;
		page_byte_array *out_arr;

		if ((size_t) buffered_arr->get_offset() < get_graph().get_in_part_size()) {
			in_arr = buffered_arr;
			assert((size_t) arr.get_offset() >= get_graph().get_in_part_size());
			out_arr = &arr;
		}
		else {
			out_arr = buffered_arr;
			assert((size_t) arr.get_offset() < get_graph().get_in_part_size());
			in_arr = &arr;
		}

		run_on_arrays(*in_arr, *out_arr);
		page_byte_array::destroy(buffered_arr);
	}
	else if (type == BOTH_EDGES) {
		buffered_arr = arr.clone();
	}
	else if (type == IN_EDGE) {
		assert((size_t) arr.get_offset() < get_graph().get_in_part_size());
		run_on_array(arr);
	}
	else if (type == OUT_EDGE) {
		assert((size_t) arr.get_offset() >= get_graph().get_in_part_size());
		run_on_array(arr);
	}
	else
		ABORT_MSG("wrong type");
}
