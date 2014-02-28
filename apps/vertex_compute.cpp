#include "vertex_compute.h"
#include "graph_engine.h"
#include "worker_thread.h"

request_range vertex_compute::get_next_request()
{
	num_complete_issues++;

	// Get the next vertex.
	vertex_id_t id = requested_vertices[fetch_idx++];

	// Find the location of the vertex.
	compute_vertex &info = graph->get_vertex(id);
	data_loc_t loc(graph->get_file_id(), info.get_ext_mem_off());
	return request_range(loc, info.get_ext_mem_size(), READ, this);
}

void vertex_compute::request_vertices(vertex_id_t ids[], int num)
{
	if (requested_vertices.empty()) {
		requested_vertices.insert(requested_vertices.end(), ids, ids + num);
		if (!std::is_sorted(requested_vertices.begin(),
					requested_vertices.end()))
			std::sort(requested_vertices.begin(), requested_vertices.end());
	}
	else {
		std::vector<vertex_id_t> merged(
				requested_vertices.size() - fetch_idx + num);
		std::merge(requested_vertices.begin() + fetch_idx, requested_vertices.end(),
				ids, ids + num, merged.begin());
		requested_vertices.swap(merged);
	}
	fetch_idx = 0;
}

void vertex_compute::run(page_byte_array &array)
{
	ext_mem_vertex_interpreter *interpreter = graph->get_vertex_interpreter();
	stack_array<char, 64> buf(interpreter->get_vertex_size());
	const page_vertex *ext_v = interpreter->interpret(array, buf.data(),
			interpreter->get_vertex_size());
	bool completed = false;
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	t->set_curr_vertex_compute(this);
	// If the algorithm doesn't need to get the full information
	// of their neighbors
	if (graph->get_required_neighbor_type() == edge_type::NONE
			// Or we haven't perform computation on the vertex yet.
			|| v == NULL) {
		v = &graph->get_vertex(ext_v->get_id());
		completed = v->run(*graph, ext_v);
	}
	else {
		num_complete_fetched++;
		completed = v->run_on_neighbors(*graph, &ext_v, 1);
	}
	// We need to notify the thread that initiate processing the vertex
	// of the completion of the vertex.
	if (completed) {
		// If the vertex has completed computation, the user compute should
		// also finish its computation.
		assert(has_completed());
		issue_thread->complete_vertex(*v);
	}
	t->reset_curr_vertex_compute();
}

class comp_ts_vertex_request
{
public:
	bool operator()(const ts_vertex_request &req1,
			const ts_vertex_request &req2) {
		return req1.get_id() < req2.get_id();
	}
};

void ts_vertex_compute::request_partial_vertices(vertex_request *reqs[], int num)
{
	for (int i = 0; i < num; i++) {
		this->reqs.push_back(*(ts_vertex_request *) reqs[i]);
	}
	if (!std::is_sorted(this->reqs.begin(), this->reqs.end(),
				comp_ts_vertex_request()))
		std::sort(this->reqs.begin(), this->reqs.end(),
				comp_ts_vertex_request());
}

request_range ts_vertex_compute::get_next_request()
{
	if (vertex_compute::has_requests())
		return vertex_compute::get_next_request();
	else {
		ts_vertex_request ts_req = reqs[fetch_idx++];
		compute_vertex &info = get_graph().get_vertex(ts_req.get_id());
		data_loc_t loc(get_graph().get_file_id(), info.get_ext_mem_off());
		// There is some overhead to fetch part of a vertex, so we should
		// minize the number of vertices fetched partially.
		// If a vertex is small enough (stored on <= 3 pages), we fetch the entire
		// vertex.
		off_t start_pg = ROUND_PAGE(info.get_ext_mem_off());
		off_t end_pg = ROUNDUP_PAGE(info.get_ext_mem_off() + info.get_ext_mem_size());
		if (end_pg - start_pg <= PAGE_SIZE * 3) {
			// We need to increase the number of issues, so we know when
			// the user task is completed.
			num_complete_issues++;
			return request_range(loc, info.get_ext_mem_size(), READ, this);
		}

		worker_thread *t = (worker_thread *) thread::get_curr_thread();
		compute_allocator *alloc = t->get_part_compute_allocator();
		assert(alloc);
		part_ts_vertex_compute *comp = (part_ts_vertex_compute *) alloc->alloc();
		comp->init(v, ts_req);
		// We assume the header of a ts-vertex is never larger than a page.
		return request_range(loc, PAGE_SIZE, READ, comp);
	}
}

request_range part_ts_vertex_compute::get_next_request()
{
	assert(required_vertex_header);
	assert(num_issued == 0);
	num_issued++;

	compute_vertex &info = graph->get_vertex(
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
	bool completed = false;
	if (required_vertex_header == NULL) {
		ext_mem_vertex_interpreter *interpreter = graph->get_vertex_interpreter();
		char *buf = new char[interpreter->get_vertex_size()];
		required_vertex_header = (const TS_page_vertex *) interpreter->interpret(
				array, buf, interpreter->get_vertex_size());
		assert(!required_vertex_header->is_complete());
	}
	else {
		ext_mem_vertex_interpreter *interpreter = graph->get_vertex_interpreter();
		stack_array<char, 64> buf(interpreter->get_vertex_size());
		const page_vertex *ext_v = interpreter->interpret_part(
				required_vertex_header, array, buf.data(),
				interpreter->get_vertex_size());

		num_fetched++;
		assert(comp_v);
		completed = comp_v->run_on_neighbors(*graph, &ext_v, 1);

		char *tmp = (char *) required_vertex_header;
		delete [] tmp;
	}
	// We need to notify the thread that initiate processing the vertex
	// of the completion of the vertex.
	if (completed)
		issue_thread->complete_vertex(*comp_v);
}
