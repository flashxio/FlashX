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
#include "vertex_compute.h"
#include "graph_engine.h"
#include "worker_thread.h"

request_range vertex_compute::get_next_request()
{
	num_complete_issues++;

	// Get the next vertex.
	vertex_id_t id = requested_vertices.top();
	requested_vertices.pop();

	// Find the location of the vertex.
	const in_mem_vertex_info info = graph->get_vertex_info(id);
	data_loc_t loc(graph->get_file_id(), info.get_ext_mem_off());
	return request_range(loc, info.get_ext_mem_size(), READ, this);
}

void vertex_compute::request_vertices(vertex_id_t ids[], size_t num)
{
	// Logging the requests.
	if (graph->get_logger()) {
		std::vector<request_range> reqs;
		for (size_t i = 0; i < num; i++) {
			const in_mem_vertex_info info = graph->get_vertex_info(ids[i]);
			data_loc_t loc(graph->get_file_id(), info.get_ext_mem_off());
			request_range req(loc, info.get_ext_mem_size(), READ, this);
			reqs.push_back(req);
		}
		graph->get_logger()->log(reqs.data(), reqs.size());
	}

	for (size_t i = 0; i < num; i++)
		requested_vertices.push(ids[i]);
}

void vertex_compute::run(page_byte_array &array)
{
	ext_mem_vertex_interpreter &interpreter = graph->get_vertex_interpreter();
	stack_array<char, 64> buf(interpreter.get_vertex_size());
	const page_vertex *ext_v = interpreter.interpret(array, buf.data(),
			interpreter.get_vertex_size());
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	t->set_curr_vertex_compute(this);
	vertex_program &curr_vprog = t->get_vertex_program();
	curr_vprog.run(*v, *ext_v);
	complete_request();
	t->reset_curr_vertex_compute();
}

void vertex_compute::complete_request()
{
	num_complete_fetched++;
	// We need to notify the thread that initiate processing the vertex
	// of the completion of the vertex.
	// TODO is this a right way to complete a vertex?
	if (has_completed())
		issue_thread->complete_vertex(*v);
}

void part_directed_vertex_compute::run(page_byte_array &array)
{
	assert(comp_v);
	vsize_t num_in_edges = 0;
	vsize_t num_out_edges = 0;
	if (req.get_type() == edge_type::IN_EDGE)
		num_in_edges = array.get_size() / sizeof(vertex_id_t);
	else if (req.get_type() == edge_type::OUT_EDGE)
		num_out_edges = array.get_size() / sizeof(vertex_id_t);
	else
		assert(0);
	page_directed_vertex pg_v(req.get_id(), num_in_edges,
			num_out_edges, array);

	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	vertex_program &curr_vprog = t->get_vertex_program();
	curr_vprog.run(*comp_v, pg_v);
	assert(t->get_curr_vertex_compute() == NULL);
	num_fetched++;
	compute->complete_request();
	compute->dec_ref();
	// If the original user compute hasn't been issued to the filesystem,
	// it's possible that the reference count on it reaches 0 here.
	// We should deallocate the user compute if its ref count reaches 0.
	if (compute->get_ref() == 0) {
		assert(compute->has_completed());
		compute_allocator *alloc = compute->get_allocator();
		alloc->free(compute);
	}
}

request_range directed_vertex_compute::get_next_request()
{
	if (vertex_compute::has_requests())
		return vertex_compute::get_next_request();
	else {
		// We need to increase the number of issues, so we know when
		// the user task is completed.
		num_complete_issues++;

		directed_vertex_request req = reqs.top();
		reqs.pop();
		const in_mem_vertex_info info = get_graph().get_vertex_info(req.get_id());

		off_t start_pg = ROUND_PAGE(info.get_ext_mem_off());
		off_t end_pg = ROUNDUP_PAGE(info.get_ext_mem_off() + info.get_ext_mem_size());
		if (end_pg - start_pg <= PAGE_SIZE
				|| req.get_type() == edge_type::BOTH_EDGES) {
			data_loc_t loc(get_graph().get_file_id(), info.get_ext_mem_off());
			request_range req(loc, info.get_ext_mem_size(), READ, this);
			if (graph->get_logger())
				graph->get_logger()->log(&req, 1);
			return req;
		}
		else {
			worker_thread *t = (worker_thread *) thread::get_curr_thread();
			compute_allocator *alloc = t->get_part_compute_allocator();
			assert(alloc);
			part_directed_vertex_compute *comp
				= (part_directed_vertex_compute *) alloc->alloc();
			comp->init((compute_directed_vertex *) v, this, req);
			compute_directed_vertex &directed_v
				= (compute_directed_vertex &) get_graph().get_vertex(req.get_id());
			vsize_t num_in_edges = directed_v.get_num_in_edges();
			vsize_t num_out_edges = directed_v.get_num_out_edges();
			if (req.get_type() == edge_type::IN_EDGE) {
				data_loc_t loc(get_graph().get_file_id(), info.get_ext_mem_off()
						+ ext_mem_directed_vertex::get_header_size());
				assert(num_in_edges > 0);
				request_range req(loc, num_in_edges * sizeof(vertex_id_t),
						READ, comp);
				if (graph->get_logger())
					graph->get_logger()->log(&req, 1);
				return req;
			}
			else {
				data_loc_t loc(get_graph().get_file_id(), info.get_ext_mem_off()
						+ ext_mem_directed_vertex::get_header_size()
						+ num_in_edges * sizeof(vertex_id_t));
				assert(num_out_edges > 0);
				request_range req(loc, num_out_edges * sizeof(vertex_id_t),
						READ, comp);
				if (graph->get_logger())
					graph->get_logger()->log(&req, 1);
				return req;
			}
		}
	}
}

class empty_page_byte_array: public page_byte_array
{
public:
	virtual void lock() {
	}
	virtual void unlock() {
	}

	virtual off_t get_offset_in_first_page() const {
		return 0;
	}
	virtual thread_safe_page *get_page(int idx) const {
		return NULL;
	}
	virtual size_t get_size() const {
		return 0;
	}
};

void directed_vertex_compute::request_partial_vertices(
		directed_vertex_request reqs[], size_t num)
{
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	vertex_program &curr_vprog = t->get_vertex_program();

	assert(this->reqs.empty());
	for (size_t i = 0; i < num; i++) {
		directed_vertex_request *req = &reqs[i];
		compute_directed_vertex &info
			= (compute_directed_vertex &) get_graph().get_vertex(req->get_id());

		// If the requested edge list is empty, we can just serve the request now.
		if ((req->get_type() == edge_type::IN_EDGE
					&& info.get_num_in_edges() == 0)
				|| (req->get_type() == edge_type::OUT_EDGE
					&& info.get_num_out_edges() == 0)) {
			empty_page_byte_array array;
			page_directed_vertex pg_v(req->get_id(), 0, 0, array);
			curr_vprog.run(*v, pg_v);
			// We don't need to call complete_request() here.
		}
		else
			this->reqs.push(*req);
	}
	// We need to notify the thread that initiate processing the vertex
	// of the completion of the vertex.
	// TODO is this a right way to complete a vertex?
	if (this->reqs.empty() && has_completed())
		issue_thread->complete_vertex(*v);
}

void ts_vertex_compute::request_partial_vertices(ts_vertex_request reqs[],
		size_t num)
{
	for (size_t i = 0; i < num; i++)
		this->reqs.push(reqs[i]);
}

request_range ts_vertex_compute::get_next_request()
{
	if (vertex_compute::has_requests())
		return vertex_compute::get_next_request();
	else {
		// We need to increase the number of issues, so we know when
		// the user task is completed.
		num_complete_issues++;

		ts_vertex_request ts_req = reqs.top();
		reqs.pop();
		const in_mem_vertex_info info = get_graph().get_vertex_info(ts_req.get_id());
		data_loc_t loc(get_graph().get_file_id(), info.get_ext_mem_off());
		// There is some overhead to fetch part of a vertex, so we should
		// minize the number of vertices fetched partially.
		// If a vertex is small enough (stored on <= 3 pages), we fetch the entire
		// vertex.
		off_t start_pg = ROUND_PAGE(info.get_ext_mem_off());
		off_t end_pg = ROUNDUP_PAGE(info.get_ext_mem_off() + info.get_ext_mem_size());
		if (end_pg - start_pg <= PAGE_SIZE * 3)
			return request_range(loc, info.get_ext_mem_size(), READ, this);

		worker_thread *t = (worker_thread *) thread::get_curr_thread();
		compute_allocator *alloc = t->get_part_compute_allocator();
		assert(alloc);
		part_ts_vertex_compute *comp = (part_ts_vertex_compute *) alloc->alloc();
		comp->init(v, this, ts_req);
		// We assume the header of a ts-vertex is never larger than a page.
		return request_range(loc, PAGE_SIZE, READ, comp);
	}
}

request_range part_ts_vertex_compute::get_next_request()
{
	assert(required_vertex_header);
	assert(num_issued == 0);
	num_issued++;

	const in_mem_vertex_info info = graph->get_vertex_info(
			required_vertex_header->get_id());
	offset_pair rel_offsets = required_vertex_header->get_edge_list_offset(
			required_part.get_range());
	data_loc_t loc(graph->get_file_id(), rel_offsets.first
			+ info.get_ext_mem_off());
	return request_range(loc, rel_offsets.second - rel_offsets.first,
			READ, this);
}

void part_ts_vertex_compute::run(page_byte_array &array)
{
	assert(!has_completed());
	if (required_vertex_header == NULL) {
		ext_mem_vertex_interpreter &interpreter = graph->get_vertex_interpreter();
		char *buf = new char[interpreter.get_vertex_size()];
		required_vertex_header = (const TS_page_vertex *) interpreter.interpret(
				array, buf, interpreter.get_vertex_size());
		assert(!required_vertex_header->is_complete());
	}
	else {
		ext_mem_vertex_interpreter &interpreter = graph->get_vertex_interpreter();
		stack_array<char, 64> buf(interpreter.get_vertex_size());
		const page_vertex *ext_v = interpreter.interpret_part(
				required_vertex_header, array, buf.data(),
				interpreter.get_vertex_size());

		num_fetched++;
		assert(comp_v);

		worker_thread *t = (worker_thread *) thread::get_curr_thread();
		vertex_program &curr_vprog = t->get_vertex_program();
		curr_vprog.run(*comp_v, *ext_v);
		ts_compute->complete_request();
		ts_compute->dec_ref();
		// If the original user compute hasn't been issued to the filesystem,
		// it's possible that the reference count on it reaches 0 here.
		// We should deallocate the user compute if its ref count reaches 0.
		if (ts_compute->get_ref() == 0) {
			assert(ts_compute->has_completed());
			compute_allocator *alloc = ts_compute->get_allocator();
			alloc->free(ts_compute);
		}
		worker_thread *curr = (worker_thread *) thread::get_curr_thread();
		// Let's just assume the user doesn't issue new requests here.
		// It's easy to change it to work with the case that the user
		// wants to issue new requests.
		assert(curr->get_curr_vertex_compute() == NULL);

		char *tmp = (char *) required_vertex_header;
		delete [] tmp;
	}
}
