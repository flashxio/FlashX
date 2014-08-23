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

#include "triangle_shared.h"

namespace {

struct undirected_runtime_data_t: public runtime_data_t
{
	vsize_t num_edge_reqs;
	vsize_t degree;

	undirected_runtime_data_t(vsize_t num_exist_triangles,
			vsize_t degree): runtime_data_t(degree, num_exist_triangles) {
		this->degree = degree;
		num_edge_reqs = 0;
	}
};

class undirected_triangle_vertex: public compute_vertex
{
	multi_func_value local_value;

	void inc_num_triangles(size_t num) {
		if (local_value.has_num_triangles())
			local_value.inc_num_triangles(num);
		else
			local_value.get_runtime_data()->num_triangles += num;
	}
public:
	undirected_triangle_vertex(vertex_id_t id): compute_vertex(id) {}

	int count_triangles(const page_vertex *v) const;

	int get_result() const {
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

	void run_on_vertex_header(vertex_program &prog, const vertex_header &header);

	void destroy_runtime() {
		undirected_runtime_data_t *data
			= (undirected_runtime_data_t *) local_value.get_runtime_data();
		size_t num_curr_triangles = data->num_triangles;
		delete data;
		local_value.set_num_triangles(num_curr_triangles);
	}
};

void undirected_triangle_vertex::run_on_vertex_header(vertex_program &prog,
		const vertex_header &header)
{
	undirected_runtime_data_t *data
		= (undirected_runtime_data_t *) local_value.get_runtime_data();
	data->num_edge_reqs++;

	if ((header.get_num_edges() < data->degree && header.get_id() != get_id())
			|| (header.get_num_edges() == data->degree
				&& header.get_id() < get_id())) {
		data->edges.push_back(header.get_id());
		data->num_required++;
	}
	if (data->num_edge_reqs == data->degree) {
		if (data->edges.empty()) {
			long ret = num_completed_vertices.inc(1);
			if (ret % 100000 == 0)
				printf("%ld completed vertices\n", ret);
			destroy_runtime();
			return;
		}
		std::sort(data->edges.begin(), data->edges.end());
		data->finalize_init();
		// We now can request the neighbors.
		request_vertices(data->edges.data(), data->edges.size());
	}
}

void undirected_triangle_vertex::run_on_itself(vertex_program &prog,
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

	std::vector<vertex_id_t> edges(vertex.get_num_edges(edge_type::IN_EDGE));
	vertex.read_edges(edge_type::IN_EDGE, edges.data(), edges.size());

	local_value.set_runtime_data(new undirected_runtime_data_t(
				local_value.get_num_triangles(),
				vertex.get_num_edges(edge_type::IN_EDGE)));
	request_vertex_headers(edges.data(), edges.size());
}

void undirected_triangle_vertex::run_on_neighbor(vertex_program &prog,
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

int undirected_triangle_vertex::count_triangles(const page_vertex *v) const
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
	page_byte_array::const_iterator<vertex_id_t> other_it
		= v->get_neigh_begin(edge_type::OUT_EDGE);
	page_byte_array::const_iterator<vertex_id_t> other_end
		= std::lower_bound(other_it, v->get_neigh_end(edge_type::OUT_EDGE),
				v->get_id());
	size_t num_v_edges = other_end - other_it;
	if (num_v_edges == 0)
		return 0;

	runtime_data_t *data = local_value.get_runtime_data();
	if (data->edge_set.size() > 0
			&& data->edges.size() > HASH_SEARCH_RATIO * num_v_edges) {
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
	else if (num_v_edges / data->edges.size() > BIN_SEARCH_RATIO) {
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
		std::vector<vertex_id_t>::const_iterator this_it = data->edges.cbegin();
		std::vector<vertex_id_t>::const_iterator this_end
			= std::lower_bound(this_it, data->edges.cend(), v->get_id());
		std::vector<int>::iterator count_it = data->triangles.begin();

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

}

#include "save_result.h"
FG_vector<size_t>::ptr compute_undirected_triangles(FG_graph::ptr fg)
{
	printf("undirected triangle counting starts\n");
	graph_index::ptr index = NUMA_graph_index<undirected_triangle_vertex>::create(
			fg->get_index_file());
	graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
			index, fg->get_configs());

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
	graph->query_on_all(vertex_query::ptr(
				new save_query<size_t, undirected_triangle_vertex>(vec)));
	return vec;
}
