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
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include "log.h"

#include "triangle_shared.h"

using namespace fg;

namespace {

struct directed_runtime_data_t: public runtime_data_t
{
	std::vector<vertex_id_t> in_edges;
	std::vector<vertex_id_t> out_edges;

	std::vector<vertex_id_t> selected_in_edges;
	std::vector<vertex_id_t> selected_out_edges;

	vsize_t num_edge_reqs;
	vsize_t num_tot_edge_reqs;
	vsize_t degree;

	directed_runtime_data_t(vsize_t num_exist_triangles, vsize_t num_in_edges,
			vsize_t degree): runtime_data_t(num_in_edges, num_exist_triangles) {
		this->num_tot_edge_reqs = 0;
		this->degree = degree;
		num_edge_reqs = 0;
	}

	bool is_in_edge(vertex_id_t id) const {
		return std::binary_search(in_edges.begin(), in_edges.end(), id);
	}

	bool is_out_edge(vertex_id_t id) const {
		return std::binary_search(out_edges.begin(), out_edges.end(), id);
	}

	void run_on_vertex_header(vertex_program &prog, compute_directed_vertex &v,
			const vertex_header &header);
	void run_on_itself(compute_directed_vertex &v, const page_vertex &vertex);

	size_t count_triangles(vertex_program &prog,
			compute_directed_vertex &directed_v, const page_vertex &v);

	bool collect_all_headers() const {
		return this->num_edge_reqs == this->num_tot_edge_reqs;
	}

	bool is_complete() const {
		return num_joined == num_required;
	}
};

void directed_runtime_data_t::run_on_itself(compute_directed_vertex &v,
		const page_vertex &vertex)
{
	this->in_edges.resize(vertex.get_num_edges(edge_type::IN_EDGE));
	vertex.read_edges(edge_type::IN_EDGE, this->in_edges.data(),
			this->in_edges.size());
	this->out_edges.resize(vertex.get_num_edges(edge_type::OUT_EDGE));
	vertex.read_edges(edge_type::OUT_EDGE, this->out_edges.data(),
			this->out_edges.size());

	// Request the number of edges of its neighbors.
	std::vector<vertex_id_t> edges;
	unique_merge(this->in_edges, this->out_edges, edges);
	this->num_tot_edge_reqs = edges.size();
	v.request_vertex_headers(edges.data(), edges.size());
}

void directed_runtime_data_t::run_on_vertex_header(vertex_program &prog,
		compute_directed_vertex &v, const vertex_header &header)
{
	vertex_id_t id = prog.get_vertex_id(v);
	vsize_t num_edges = header.get_num_edges();
	this->num_edge_reqs++;
	if (this->is_in_edge(header.get_id())) {
		if ((num_edges < this->degree && header.get_id() != id)
				|| (num_edges == this->degree && header.get_id() < id)) {
			this->selected_in_edges.push_back(header.get_id());
		}
	}
	if (this->is_out_edge(header.get_id())) {
		if ((num_edges < this->degree && header.get_id() != id)
				|| (num_edges == this->degree && header.get_id() < id)) {
			this->selected_out_edges.push_back(header.get_id());
		}
	}

	if (this->num_edge_reqs == this->num_tot_edge_reqs) {
		if (this->selected_in_edges.empty()
				|| this->selected_out_edges.empty()) {
			num_required = 0;
		}
		else {
			std::sort(this->selected_in_edges.begin(),
					this->selected_in_edges.end());
			std::sort(this->selected_out_edges.begin(),
					this->selected_out_edges.end());
			this->edges = this->selected_in_edges;
			this->num_required = this->selected_out_edges.size();
			std::vector<directed_vertex_request> reqs(
					this->selected_out_edges.size());
			for (size_t i = 0; i < this->selected_out_edges.size(); i++) {
				vertex_id_t id = this->selected_out_edges[i];
				reqs[i] = directed_vertex_request(id, edge_type::OUT_EDGE);
			}
			this->finalize_init();
			this->in_edges.clear();
			this->in_edges.shrink_to_fit();
			this->out_edges.clear();
			this->out_edges.shrink_to_fit();
			this->selected_in_edges.clear();
			this->selected_in_edges.shrink_to_fit();
			this->selected_out_edges.clear();
			this->selected_out_edges.shrink_to_fit();
			v.request_partial_vertices(reqs.data(), reqs.size());
		}
	}
}

size_t directed_runtime_data_t::count_triangles(vertex_program &prog,
		compute_directed_vertex &directed_v, const page_vertex &v)
{
	vertex_id_t id = prog.get_vertex_id(directed_v);
	size_t num_local_triangles = 0;
	assert(v.get_id() != id);

	if (v.get_num_edges(edge_type::OUT_EDGE) == 0)
		return 0;

	/*
	 * We search for triangles with two different ways:
	 * binary search if two adjacency lists have very different sizes,
	 * scan otherwise.
	 *
	 * when binary search for multiple neighbors, we can reduce binary search
	 * overhead by using the new end in the search range. We can further reduce
	 * overhead by searching in a reverse order (start from the largest neighbor).
	 * Since vertices of smaller ID has more neighbors, it's more likely
	 * that a neighbor is in the beginning of the adjacency list, and
	 * the search range will be narrowed faster.
	 */

	if (this->edge_set.size() > 0
			&& this->edges.size() > HASH_SEARCH_RATIO * v.get_num_edges(
				edge_type::OUT_EDGE)) {
		edge_iterator other_it = v.get_neigh_begin(edge_type::OUT_EDGE);
		edge_iterator other_end = v.get_neigh_end(edge_type::OUT_EDGE);
		for (; other_it != other_end; ++other_it) {
			vertex_id_t neigh_neighbor = *other_it;
			runtime_data_t::edge_set_t::const_iterator it
				= this->edge_set.find(neigh_neighbor);
			if (it != this->edge_set.end()) {
				if (neigh_neighbor != v.get_id() && neigh_neighbor != id) {
					num_local_triangles++;
					int idx = (*it).get_idx();
					this->triangles[idx]++;
				}
			}
		}
	}
	// If the neighbor vertex has way more edges than this vertex.
	else if (v.get_num_edges(edge_type::OUT_EDGE) / this->edges.size(
				) > BIN_SEARCH_RATIO) {
		edge_iterator other_it = v.get_neigh_begin(edge_type::OUT_EDGE);
		edge_iterator other_end = v.get_neigh_end(edge_type::OUT_EDGE);
		for (int i = this->edges.size() - 1; i >= 0; i--) {
			vertex_id_t this_neighbor = this->edges.at(i);
			// We need to skip loops.
			if (this_neighbor != v.get_id() && this_neighbor != id) {
				edge_iterator first = std::lower_bound(other_it, other_end,
						this_neighbor);
				if (first != other_end && this_neighbor == *first) {
					num_local_triangles++;
					this->triangles[i]++;
				}
				other_end = first;
			}
		}
	}
	else {
		std::vector<vertex_id_t>::const_iterator this_it = this->edges.begin();
		std::vector<int>::iterator count_it = this->triangles.begin();
		std::vector<vertex_id_t>::const_iterator this_end = this->edges.end();
		edge_seq_iterator other_it = v.get_neigh_seq_it(edge_type::OUT_EDGE, 0,
					v.get_num_edges(edge_type::OUT_EDGE));
		while (this_it != this_end && other_it.has_next()) {
			vertex_id_t this_neighbor = *this_it;
			vertex_id_t neigh_neighbor = other_it.curr();
			if (this_neighbor == neigh_neighbor) {
				// skip loop
				if (neigh_neighbor != v.get_id() && neigh_neighbor != id) {
					num_local_triangles++;
					(*count_it)++;
				}
				++this_it;
				other_it.next();
				++count_it;
			}
			else if (this_neighbor < neigh_neighbor) {
				++this_it;
				++count_it;
			}
			else
				other_it.next();
		}
	}
	return num_local_triangles;
}

class directed_triangle_vertex: public compute_directed_vertex
{
	triangle_multi_func_value local_value;

	void inc_num_triangles(size_t num) {
		if (local_value.has_num_triangles())
			local_value.inc_num_triangles(num);
		else
			local_value.get_runtime_data()->num_triangles += num;
	}
public:
	directed_triangle_vertex(vertex_id_t id): compute_directed_vertex(id) {}

	size_t get_result() const {
		return local_value.get_num_triangles();
	}

	void run(vertex_program &prog) {
		vertex_id_t id = prog.get_vertex_id(*this);
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex) {
		if (vertex.get_id() == prog.get_vertex_id(*this))
			run_on_itself(prog, vertex);
		else
			run_on_neighbor(prog, vertex);
	}

	void run_on_itself(vertex_program &prog, const page_vertex &vertex);
	void run_on_neighbor(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
		inc_num_triangles(((count_msg &) msg).get_num());
	}

	void run_on_vertex_header(vertex_program &prog, const vertex_header &header);

	void destroy_runtime() {
		directed_runtime_data_t *data
			= (directed_runtime_data_t *) local_value.get_runtime_data();
		size_t num_curr_triangles = data->num_triangles;
		delete data;
		local_value.set_num_triangles(num_curr_triangles);
	}
};

void directed_triangle_vertex::run_on_vertex_header(vertex_program &prog,
		const vertex_header &header)
{
	directed_runtime_data_t *data
		= (directed_runtime_data_t *) local_value.get_runtime_data();
	data->run_on_vertex_header(prog, *this, header);

	if (data->collect_all_headers() && data->is_complete()) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			BOOST_LOG_TRIVIAL(debug)
				<< boost::format("%1% completed vertices") % ret;
		destroy_runtime();
	}
}

void directed_triangle_vertex::run_on_itself(vertex_program &prog,
		const page_vertex &vertex)
{
	assert(!local_value.has_runtime_data());

	long ret = num_working_vertices.inc(1);
	if (ret % 100000 == 0)
		BOOST_LOG_TRIVIAL(debug)
			<< boost::format("%1% working vertices") % ret;
	// A vertex has to have in-edges and out-edges in order to form
	// a triangle. so we can simply skip the vertices that don't have
	// either of them.
	if (vertex.get_num_edges(edge_type::OUT_EDGE) == 0
			|| vertex.get_num_edges(edge_type::IN_EDGE) == 0) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			BOOST_LOG_TRIVIAL(debug)
				<< boost::format("%1% completed vertices") % ret;
		return;
	}

	directed_runtime_data_t *data = new directed_runtime_data_t(
				local_value.get_num_triangles(),
				vertex.get_num_edges(edge_type::IN_EDGE),
				vertex.get_num_edges(edge_type::BOTH_EDGES));
	local_value.set_runtime_data(data);
	data->run_on_itself(*this, vertex);
}

void directed_triangle_vertex::run_on_neighbor(vertex_program &prog,
		const page_vertex &vertex)
{
	assert(local_value.has_runtime_data());
	directed_runtime_data_t *data
		= (directed_runtime_data_t *) local_value.get_runtime_data();
	data->num_joined++;
	int ret = data->count_triangles(prog, *this, vertex);
	// If we find triangles with the neighbor, notify the neighbor
	// as well.
	if (ret > 0) {
		inc_num_triangles(ret);
		count_msg msg(ret);
		prog.send_msg(vertex.get_id(), msg);
	}

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (data->is_complete()) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			BOOST_LOG_TRIVIAL(debug)
				<< boost::format("%1% completed vertices") % ret;

		// Inform all neighbors in the in-edges.
		for (size_t i = 0; i < data->triangles.size(); i++) {
			// Inform the neighbor if they share triangles.
			if (data->triangles[i] > 0) {
				count_msg msg(data->triangles[i]);
				prog.send_msg(data->edges[i], msg);
			}
		}
		destroy_runtime();
	}
}

}

#include "save_result.h"

namespace fg
{

FG_vector<size_t>::ptr compute_directed_triangles(FG_graph::ptr fg,
		directed_triangle_type type)
{
	bool directed = fg->get_graph_header().is_directed_graph();
	if (!directed) {
		BOOST_LOG_TRIVIAL(error)
			<< "This algorithm counts triangles in a directed graph";
		return FG_vector<size_t>::ptr();
	}

	graph_index::ptr index = NUMA_graph_index<directed_triangle_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);

	BOOST_LOG_TRIVIAL(info) << "triangle counting starts";
	BOOST_LOG_TRIVIAL(info) << "prof_file: " << graph_conf.get_prof_file();
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("It takes %1% seconds to count all triangles")
		% time_diff(start, end);

	FG_vector<size_t>::ptr vec = FG_vector<size_t>::create(graph);
	graph->query_on_all(vertex_query::ptr(
				new save_query<size_t, directed_triangle_vertex>(vec)));
	return vec;
}

}
