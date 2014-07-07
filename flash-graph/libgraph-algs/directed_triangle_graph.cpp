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
#include <google/profiler.h>
#endif

#include "triangle_shared.h"

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
};

class directed_triangle_vertex: public compute_directed_vertex
{
	multi_func_value local_value;

	void inc_num_triangles(size_t num) {
		if (local_value.has_num_triangles())
			local_value.inc_num_triangles(num);
		else
			local_value.get_runtime_data()->num_triangles += num;
	}
public:
	directed_triangle_vertex(vertex_id_t id): compute_directed_vertex(id) {}

	int count_triangles(const page_vertex *v) const;

	int get_num_triangles() const {
		return local_value.get_num_triangles();
	}

	void run(vertex_program &prog) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex) {
		if (vertex.get_id() == get_id())
			run_on_itself(prog, vertex);
		else
			run_on_neighbor(prog, vertex);
	}

	void run_on_itself(vertex_program &prog, const page_vertex &vertex);
	void run_on_neighbor(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
		inc_num_triangles(((count_msg &) msg).get_num());
	}

	void run_on_num_dedges(vertex_id_t id, vsize_t num_in_edges,
			vsize_t num_out_edges);

	void destroy_runtime() {
		directed_runtime_data_t *data
			= (directed_runtime_data_t *) local_value.get_runtime_data();
		size_t num_curr_triangles = data->num_triangles;
		delete data;
		local_value.set_num_triangles(num_curr_triangles);
	}
};

void directed_triangle_vertex::run_on_num_dedges(vertex_id_t id,
		vsize_t num_in_edges, vsize_t num_out_edges)
{
	vsize_t num_edges = num_in_edges + num_out_edges;
	directed_runtime_data_t *data
		= (directed_runtime_data_t *) local_value.get_runtime_data();
	data->num_edge_reqs++;
	if (data->is_in_edge(id)) {
		if ((num_edges < data->degree && id != this->get_id())
				|| (num_edges == data->degree && id < this->get_id())) {
			data->selected_in_edges.push_back(id);
		}
	}
	if (data->is_out_edge(id)) {
		if ((num_edges < data->degree && id != this->get_id())
				|| (num_edges == data->degree && id < this->get_id())) {
			data->selected_out_edges.push_back(id);
		}
	}

	if (data->num_edge_reqs == data->num_tot_edge_reqs) {
		if (data->selected_in_edges.empty()
				|| data->selected_out_edges.empty()) {
			long ret = num_completed_vertices.inc(1);
			if (ret % 100000 == 0)
				printf("%ld completed vertices\n", ret);
			destroy_runtime();
			return;
		}
		else {
			std::sort(data->selected_in_edges.begin(),
					data->selected_in_edges.end());
			std::sort(data->selected_out_edges.begin(),
					data->selected_out_edges.end());
			data->edges = data->selected_in_edges;
			data->num_required = data->selected_out_edges.size();
			std::vector<directed_vertex_request> reqs(
					data->selected_out_edges.size());
			for (size_t i = 0; i < data->selected_out_edges.size(); i++) {
				vertex_id_t id = data->selected_out_edges[i];
				reqs[i] = directed_vertex_request(id, edge_type::OUT_EDGE);
			}
			data->finalize_init();
			data->in_edges.clear();
			data->in_edges.shrink_to_fit();
			data->out_edges.clear();
			data->out_edges.shrink_to_fit();
			data->selected_in_edges.clear();
			data->selected_in_edges.shrink_to_fit();
			data->selected_out_edges.clear();
			data->selected_out_edges.shrink_to_fit();
			request_partial_vertices(reqs.data(), reqs.size());
		}
	}
}

void directed_triangle_vertex::run_on_itself(vertex_program &prog,
		const page_vertex &vertex)
{
	assert(!local_value.has_runtime_data());

	long ret = num_working_vertices.inc(1);
	if (ret % 100000 == 0)
		printf("%ld working vertices\n", ret);
	// A vertex has to have in-edges and out-edges in order to form
	// a triangle. so we can simply skip the vertices that don't have
	// either of them.
	if (vertex.get_num_edges(edge_type::OUT_EDGE) == 0
			|| vertex.get_num_edges(edge_type::IN_EDGE) == 0) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		return;
	}

	directed_runtime_data_t *data = new directed_runtime_data_t(
				local_value.get_num_triangles(),
				vertex.get_num_edges(edge_type::IN_EDGE),
				vertex.get_num_edges(edge_type::BOTH_EDGES));
	data->in_edges.resize(vertex.get_num_edges(edge_type::IN_EDGE));
	vertex.read_edges(edge_type::IN_EDGE, data->in_edges.data(),
			data->in_edges.size());
	data->out_edges.resize(vertex.get_num_edges(edge_type::OUT_EDGE));
	vertex.read_edges(edge_type::OUT_EDGE, data->out_edges.data(),
			data->out_edges.size());
	local_value.set_runtime_data(data);

	// Request the number of edges of its neighbors.
	std::vector<vertex_id_t> edges;
	unique_merge(data->in_edges, data->out_edges, edges);
	data->num_tot_edge_reqs = edges.size();
	request_num_edges(edges.data(), edges.size());
}

void directed_triangle_vertex::run_on_neighbor(vertex_program &prog,
		const page_vertex &vertex)
{
	assert(local_value.has_runtime_data());
	runtime_data_t *data = local_value.get_runtime_data();
	data->num_joined++;
	int ret = count_triangles(&vertex);
	// If we find triangles with the neighbor, notify the neighbor
	// as well.
	if (ret > 0) {
		inc_num_triangles(ret);
		count_msg msg(ret);
		prog.send_msg(vertex.get_id(), msg);
	}

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (data->num_joined == data->num_required) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);

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

int directed_triangle_vertex::count_triangles(const page_vertex *v) const
{
	int num_local_triangles = 0;
	assert(v->get_id() != this->get_id());

	if (v->get_num_edges(edge_type::OUT_EDGE) == 0)
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

	runtime_data_t *data = local_value.get_runtime_data();
	if (data->edge_set.size() > 0
			&& data->edges.size() > HASH_SEARCH_RATIO * v->get_num_edges(
				edge_type::OUT_EDGE)) {
		page_byte_array::const_iterator<vertex_id_t> other_it
			= v->get_neigh_begin(edge_type::OUT_EDGE);
		page_byte_array::const_iterator<vertex_id_t> other_end
			= v->get_neigh_end(edge_type::OUT_EDGE);
		for (; other_it != other_end; ++other_it) {
			vertex_id_t neigh_neighbor = *other_it;
			runtime_data_t::edge_set_t::const_iterator it
				= data->edge_set.find(neigh_neighbor);
			if (it != data->edge_set.end()) {
				if (neigh_neighbor != v->get_id()
						&& neigh_neighbor != this->get_id()) {
					num_local_triangles++;
					int idx = (*it).get_idx();
					data->triangles[idx]++;
				}
			}
		}
	}
	// If the neighbor vertex has way more edges than this vertex.
	else if (v->get_num_edges(edge_type::OUT_EDGE) / data->edges.size(
				) > BIN_SEARCH_RATIO) {
		page_byte_array::const_iterator<vertex_id_t> other_it
			= v->get_neigh_begin(edge_type::OUT_EDGE);
		page_byte_array::const_iterator<vertex_id_t> other_end
			= v->get_neigh_end(edge_type::OUT_EDGE);
		for (int i = data->edges.size() - 1; i >= 0; i--) {
			vertex_id_t this_neighbor = data->edges.at(i);
			// We need to skip loops.
			if (this_neighbor != v->get_id()
					&& this_neighbor != this->get_id()) {
				page_byte_array::const_iterator<vertex_id_t> first
					= std::lower_bound(other_it, other_end, this_neighbor);
				if (first != other_end && this_neighbor == *first) {
					num_local_triangles++;
					data->triangles[i]++;
				}
				other_end = first;
			}
		}
	}
	else {
		std::vector<vertex_id_t>::const_iterator this_it = data->edges.begin();
		std::vector<int>::iterator count_it = data->triangles.begin();
		std::vector<vertex_id_t>::const_iterator this_end = data->edges.end();
		page_byte_array::seq_const_iterator<vertex_id_t> other_it
			= v->get_neigh_seq_it(edge_type::OUT_EDGE, 0,
					v->get_num_edges(edge_type::OUT_EDGE));
		while (this_it != this_end && other_it.has_next()) {
			vertex_id_t this_neighbor = *this_it;
			vertex_id_t neigh_neighbor = other_it.curr();
			if (this_neighbor == neigh_neighbor) {
				// skip loop
				if (neigh_neighbor != v->get_id()
						&& neigh_neighbor != this->get_id()) {
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

class save_ntriangles_query: public vertex_query
{
	FG_vector<size_t>::ptr vec;
public:
	save_ntriangles_query(FG_vector<size_t>::ptr vec) {
		this->vec = vec;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		directed_triangle_vertex &tv = (directed_triangle_vertex &) v;
		vec->set(tv.get_id(), tv.get_num_triangles());
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
	}

	virtual ptr clone() {
		return vertex_query::ptr(new save_ntriangles_query(vec));
	}
};

}

FG_vector<size_t>::ptr compute_directed_triangles(FG_graph::ptr fg,
		directed_triangle_type type)
{
	graph_index::ptr index = NUMA_graph_index<directed_triangle_vertex>::create(
			fg->get_index_file());
	graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
			index, fg->get_configs());

	printf("triangle counting starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
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
	printf("It takes %f seconds to count all triangles\n",
			time_diff(start, end));

	FG_vector<size_t>::ptr vec = FG_vector<size_t>::create(graph);
	graph->query_on_all(vertex_query::ptr(new save_ntriangles_query(vec)));
	return vec;
}
