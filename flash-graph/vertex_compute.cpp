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
#include "vertex_index_reader.h"

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

request_range vertex_compute::get_next_request()
{
	// Get the next vertex.
	const in_mem_vertex_info info = requested_vertices.top();
	requested_vertices.pop();
	data_loc_t loc(graph->get_file_id(), info.get_ext_mem_off());
	num_issued++;
	return request_range(loc, info.get_ext_mem_size(), READ, this);
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
	size_t header_size = issue_thread->get_graph().get_vertex_header_size();
	vsize_t num_edges = (size - header_size) / sizeof(vertex_id_t);
	vertex_header header(id, num_edges);
	issue_thread->start_run_vertex(v);
	issue_thread->get_vertex_program(v.is_part()).run_on_num_edges(*v, header);
	issue_thread->finish_run_vertex(v);
	num_edge_completed++;
	if (get_num_pending() == 0)
		issue_thread->complete_vertex(v);
}

vertex_id_t vertex_compute::get_id() const
{
	return v->get_id();
}

void vertex_compute::issue_io_request(const in_mem_vertex_info &info)
{
	// If the vertex compute has been issued to SAFS, SAFS will get the IO
	// request from the interface of user_compute. In this case, we only
	// need to add the I/O request to the queue.
	if (issued_to_io()) {
		requested_vertices.push(info);
	}
	else {
		// Otherwise, we need to issue the I/O request to SAFS explicitly.
		data_loc_t loc(graph->get_file_id(), info.get_ext_mem_off());
		io_request req(this, loc, info.get_ext_mem_size(), READ);
		num_issued++;
		issue_thread->issue_io_request(req);
	}
}

void vertex_compute::run(page_byte_array &array)
{
	ext_mem_vertex_interpreter &interpreter = graph->get_vertex_interpreter();
	stack_array<char, 64> buf(interpreter.get_vertex_size());
	const page_vertex *ext_v = interpreter.interpret(array, buf.data(),
			interpreter.get_vertex_size());
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	vertex_program &curr_vprog = t->get_vertex_program(v.is_part());
	issue_thread->start_run_vertex(v);
	curr_vprog.run(*v, *ext_v);
	issue_thread->finish_run_vertex(v);
	complete_request();
}

void vertex_compute::complete_request()
{
	num_complete_fetched++;
	// We need to notify the thread that initiate processing the vertex
	// of the completion of the vertex.
	// TODO is this a right way to complete a vertex?
	if (get_num_pending() == 0)
		issue_thread->complete_vertex(v);
}

void part_directed_vertex_compute::run(page_byte_array &array)
{
	assert(comp_v.is_valid());
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
	vertex_program &curr_vprog = t->get_vertex_program(comp_v.is_part());
	t->start_run_vertex(comp_v);
	curr_vprog.run(*comp_v, pg_v);
	t->finish_run_vertex(comp_v);
	num_fetched++;
	compute->complete_request();
}

void directed_vertex_compute::run(page_byte_array &array)
{
	if (part_req.is_valid()) {
		assert(part_req.get_offset() == array.get_offset());
		assert(part_req.get_size() == array.get_size());

		vsize_t num_in_edges = 0;
		vsize_t num_out_edges = 0;
		if (part_req.get_request().get_type() == edge_type::IN_EDGE)
			num_in_edges = array.get_size() / sizeof(vertex_id_t);
		else if (part_req.get_request().get_type() == edge_type::OUT_EDGE)
			num_out_edges = array.get_size() / sizeof(vertex_id_t);
		else
			assert(0);
		page_directed_vertex pg_v(part_req.get_request().get_id(),
				num_in_edges, num_out_edges, array);

		worker_thread *t = (worker_thread *) thread::get_curr_thread();
		vertex_program &curr_vprog = t->get_vertex_program(v.is_part());
		issue_thread->start_run_vertex(v);
		curr_vprog.run(*v, pg_v);
		issue_thread->finish_run_vertex(v);
		complete_request();

		part_req.reset();
	}
	else
		vertex_compute::run(array);
}

void directed_vertex_compute::complete_empty_part(
		const directed_vertex_request &req)
{
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	assert(t == issue_thread);
	vertex_program &curr_vprog = t->get_vertex_program(v.is_part());

	empty_page_byte_array array;
	page_directed_vertex pg_v(req.get_id(), 0, 0, array);
	issue_thread->start_run_vertex(v);
	curr_vprog.run(*v, pg_v);
	issue_thread->finish_run_vertex(v);
	complete_request();
}

int directed_vertex_compute::has_requests()
{
	if (vertex_compute::has_requests())
		return true;
	else if (reqs.size() == 0)
		return false;
	else {
		while (!reqs.empty()) {
			const part_request_info &req_info = reqs.top();
			const directed_vertex_request &req = req_info.get_request();
			const in_mem_directed_vertex_info &info = req_info.get_info();
			vsize_t num_in_edges = info.get_num_in_edges();
			vsize_t num_out_edges = info.get_num_out_edges();
			// If the requested edge list is empty, we can just serve the request now.
			// and try the next request.
			if ((req.get_type() == edge_type::IN_EDGE && num_in_edges == 0)
					|| (req.get_type() == edge_type::OUT_EDGE
						&& num_out_edges == 0)) {
				num_issued++;
				complete_empty_part(req);
				reqs.pop();
			}
			else
				return true;
		}

		return false;
	}
}

request_range directed_vertex_compute::get_next_request()
{
	if (vertex_compute::has_requests())
		return vertex_compute::get_next_request();
	else {
		const part_request_info req_info = reqs.top();
		reqs.pop();
		const directed_vertex_request &req = req_info.get_request();
		const in_mem_directed_vertex_info &info = req_info.get_info();
		request_range range = generate_request(req, info);
		num_issued++;
		return range;
	}
}

request_range directed_vertex_compute::generate_request(
		const directed_vertex_request &req,
		const in_mem_directed_vertex_info &info)
{
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
		user_compute *comp;
		if (issued_to_io()) {
			part_directed_vertex_compute *part_comp
				= (part_directed_vertex_compute *) alloc->alloc();
			part_comp->init(v, this, req);
			comp = part_comp;
		}
		else
			comp = this;
		vsize_t num_in_edges = info.get_num_in_edges();
		vsize_t num_out_edges = info.get_num_out_edges();
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

void directed_vertex_compute::issue_io_request(
		const directed_vertex_request &req,
		const in_mem_directed_vertex_info &info)
{
	// This follows the same idea as in vertex_compute.
	if (issued_to_io()) {
		part_request_info req_info(req, info);
		reqs.push(req_info);
	}
	else {
		// It's possible the requested part is empty. We can notify the user
		// immediately.
		vsize_t num_in_edges = info.get_num_in_edges();
		vsize_t num_out_edges = info.get_num_out_edges();
		if ((req.get_type() == edge_type::IN_EDGE && num_in_edges == 0)
				|| (req.get_type() == edge_type::OUT_EDGE
					&& num_out_edges == 0)) {
			num_issued++;
			complete_empty_part(req);
		}
		else {
			request_range range = generate_request(req, info);
			assert(range.get_compute() == this);
			// If we request the entire vertex.
			if (info.get_ext_mem_off() == range.get_loc().get_offset())
				assert(info.get_ext_mem_size() == range.get_size());
			else {
				assert(!part_req.is_valid());
				part_req.init(req, info);
			}
			io_request io_req(range.get_compute(), range.get_loc(),
					range.get_size(), READ);
			num_issued++;
			issue_thread->issue_io_request(io_req);
		}
	}
}

void directed_vertex_compute::request_partial_vertices(
		directed_vertex_request reqs[], size_t num)
{
	num_requested += num;
	issue_thread->get_index_reader().request_vertices(reqs, num, *this);
}

void directed_vertex_compute::run_on_num_edges(vertex_id_t id,
		vsize_t num_in_edges, vsize_t num_out_edges)
{
	directed_vertex_header header(id, num_in_edges, num_out_edges);
	issue_thread->start_run_vertex(v);
	issue_thread->get_vertex_program(v.is_part()).run_on_num_edges(*v, header);
	issue_thread->finish_run_vertex(v);
	num_edge_completed++;
	if (get_num_pending() == 0)
		issue_thread->complete_vertex(v);
}

void directed_vertex_compute::request_num_edges(vertex_id_t ids[], size_t num)
{
	num_edge_requests += num;
	issue_thread->get_index_reader().request_num_directed_edges(ids,
			num, *this);
}

#if 0
void ts_vertex_compute::request_partial_vertices(ts_vertex_request reqs[],
		size_t num)
{
	num_requested += num;
	for (size_t i = 0; i < num; i++)
		this->reqs.push(reqs[i]);
}

request_range ts_vertex_compute::get_next_request()
{
	if (vertex_compute::has_requests())
		return vertex_compute::get_next_request();
	else {
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
#endif
