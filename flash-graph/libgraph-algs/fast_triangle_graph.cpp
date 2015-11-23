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
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include "triangle_shared.h"

using namespace fg;

/*
 * This graph algorithm counts the number of triangles in a directed graph.
 * It only counts cycle triangles, in which all directed edges go to
 * the same direction.
 *     E.g A -----> B
 *         ^     /
 *         |   /
 *         | v
 *         C
 * It uses 2-D partitioning.
 */

namespace {

// The type of edges we will keep in memory.
edge_type in_mem_edge_type = OUT_EDGE;
// The type of edges from which we read neighbors.
edge_type neigh_edge_type = IN_EDGE;

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
	directed_triangle_vertex(vertex_id_t id): compute_directed_vertex(id) {
	}

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
};

runtime_data_t *construct_runtime(vertex_program &prog, const page_vertex &vertex,
		size_t num_local_edges, std::vector<vertex_id_t> &neighbors,
		vertex_id_t id_start, vertex_id_t id_end)
{
	std::vector<vertex_id_t> in_mem_edges;

	edge_seq_iterator it = vertex.get_neigh_seq_it(neigh_edge_type, 0,
				vertex.get_num_edges(neigh_edge_type));
	PAGE_FOREACH(vertex_id_t, id, it) {
		if (id < id_start || id >= id_end)
			continue;
		size_t num_local_edges1 = prog.get_num_edges(id);
		if ((num_local_edges1 < num_local_edges && id != vertex.get_id())
				|| (num_local_edges1 == num_local_edges
					&& id < vertex.get_id())) {
			neighbors.push_back(id);
		}
	} PAGE_FOREACH_END
	if (neighbors.empty()) {
#ifdef DEBUG
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			BOOST_LOG_TRIVIAL(debug)
				<< boost::format("%1% completed vertices") % ret;
#endif
		return NULL;
	}

	it = vertex.get_neigh_seq_it(in_mem_edge_type, 0,
			vertex.get_num_edges(in_mem_edge_type));
	PAGE_FOREACH(vertex_id_t, id, it) {
		size_t num_local_edges1 = prog.get_num_edges(id);
		if ((num_local_edges1 < num_local_edges && id != vertex.get_id())
				|| (num_local_edges1 == num_local_edges
					&& id < vertex.get_id())) {
			in_mem_edges.push_back(id);
		}
	} PAGE_FOREACH_END

	if (in_mem_edges.empty()) {
#ifdef DEBUG
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			BOOST_LOG_TRIVIAL(debug)
				<< boost::format("%1% completed vertices") % ret;
#endif
		return NULL;
	}

	runtime_data_t *data = new runtime_data_t(in_mem_edges.size(), 0);
	data->edges = in_mem_edges;
	data->num_required = neighbors.size();
	data->finalize_init();
	return data;
}

size_t count_triangles(runtime_data_t *data, const page_vertex &v,
		vertex_id_t this_id)
{
	size_t num_local_triangles = 0;

	if (v.get_num_edges(neigh_edge_type) == 0)
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

	if (data->edge_set.size() > 0
			&& data->edges.size() > HASH_SEARCH_RATIO * v.get_num_edges(
				neigh_edge_type)) {
		edge_iterator other_it = v.get_neigh_begin(neigh_edge_type);
		edge_iterator other_end = v.get_neigh_end(neigh_edge_type);
		for (; other_it != other_end; ++other_it) {
			vertex_id_t neigh_neighbor = *other_it;
			runtime_data_t::edge_set_t::const_iterator it
				= data->edge_set.find(neigh_neighbor);
			if (it != data->edge_set.end()) {
				if (neigh_neighbor != v.get_id()
						&& neigh_neighbor != this_id) {
					num_local_triangles++;
					int idx = (*it).get_idx();
					data->triangles[idx]++;
				}
			}
		}
	}
	// If the neighbor vertex has way more edges than this vertex.
	else if (v.get_num_edges(neigh_edge_type) / data->edges.size(
				) > BIN_SEARCH_RATIO) {
		edge_iterator other_it = v.get_neigh_begin(neigh_edge_type);
		edge_iterator other_end = v.get_neigh_end(neigh_edge_type);
		for (int i = data->edges.size() - 1; i >= 0; i--) {
			vertex_id_t this_neighbor = data->edges.at(i);
			// We need to skip loops.
			if (this_neighbor != v.get_id()
					&& this_neighbor != this_id) {
				edge_iterator first = std::lower_bound(other_it, other_end,
						this_neighbor);
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
		edge_seq_iterator other_it = v.get_neigh_seq_it(neigh_edge_type, 0,
					v.get_num_edges(neigh_edge_type));
		while (this_it != this_end && other_it.has_next()) {
			vertex_id_t this_neighbor = *this_it;
			vertex_id_t neigh_neighbor = other_it.curr();
			if (this_neighbor == neigh_neighbor) {
				// skip loop
				if (neigh_neighbor != v.get_id()
						&& neigh_neighbor != this_id) {
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

void directed_triangle_vertex::run_on_itself(vertex_program &prog,
		const page_vertex &vertex)
{
	assert(!local_value.has_runtime_data());

#ifdef DEBUG
	long ret = num_working_vertices.inc(1);
	if (ret % 100000 == 0)
		BOOST_LOG_TRIVIAL(debug)
			<< boost::format("%1% working vertices") % ret;
#endif
	// A vertex has to have in-edges and out-edges in order to form
	// a triangle. so we can simply skip the vertices that don't have
	// either of them.
	if (vertex.get_num_edges(edge_type::OUT_EDGE) == 0
			|| vertex.get_num_edges(edge_type::IN_EDGE) == 0) {
#ifdef DEBUG
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			BOOST_LOG_TRIVIAL(debug)
				<< boost::format("%1% completed vertices") % ret;
#endif
		return;
	}

	std::vector<vertex_id_t> selected_neighbors;
	runtime_data_t *data = construct_runtime(prog, vertex,
			prog.get_num_edges(vertex.get_id()), selected_neighbors, 0,
			prog.get_graph().get_max_vertex_id() + 1);
	if (data == NULL)
		return;

	// We have to set runtime data before calling request_partial_vertices.
	// It's possible that the request to a partial vertex can be completed
	// immediately and run_on_neighbor is called in request_partial_vertices.
	// TODO Maybe I should avoid that.
	data->num_triangles += local_value.get_num_triangles();
	local_value.set_runtime_data(data);
	std::vector<directed_vertex_request> reqs(selected_neighbors.size());
	for (size_t i = 0; i < selected_neighbors.size(); i++) {
		vertex_id_t id = selected_neighbors[i];
		reqs[i] = directed_vertex_request(id, neigh_edge_type);
	}
	request_partial_vertices(reqs.data(), reqs.size());
}

void directed_triangle_vertex::run_on_neighbor(vertex_program &prog,
		const page_vertex &vertex)
{
	assert(local_value.has_runtime_data());
	runtime_data_t *data = local_value.get_runtime_data();
	data->num_joined++;
	size_t ret = count_triangles(local_value.get_runtime_data(), vertex,
			prog.get_vertex_id(*this));
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
#ifdef DEBUG
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			BOOST_LOG_TRIVIAL(debug)
				<< boost::format("%1% completed vertices") % ret;
#endif

		// Inform all neighbors in the in-edges.
		for (size_t i = 0; i < data->triangles.size(); i++) {
			// Inform the neighbor if they share triangles.
			if (data->triangles[i] > 0) {
				count_msg msg(data->triangles[i]);
				prog.send_msg(data->edges[i], msg);
			}
		}
		size_t num_curr_triangles = data->num_triangles;
		delete data;
		local_value.set_num_triangles(num_curr_triangles);
	}
}

class part_directed_triangle_vertex: public part_compute_directed_vertex
{
	triangle_multi_func_value local_value;

	void inc_num_triangles(size_t num) {
		if (local_value.has_num_triangles())
			local_value.inc_num_triangles(num);
		else
			local_value.get_runtime_data()->num_triangles += num;
	}
public:
	part_directed_triangle_vertex(vertex_id_t id,
			int part_id): part_compute_directed_vertex(id, part_id) {
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
};

void part_directed_triangle_vertex::run_on_itself(vertex_program &prog,
		const page_vertex &vertex)
{
	assert(!local_value.has_runtime_data());

#ifdef DEBUG
	long ret = num_working_vertices.inc(1);
	if (ret % 100000 == 0)
		BOOST_LOG_TRIVIAL(debug)
			<< boost::format("%1% working vertices") % ret;
#endif
	// A vertex has to have in-edges and out-edges in order to form
	// a triangle. so we can simply skip the vertices that don't have
	// either of them.
	if (vertex.get_num_edges(edge_type::OUT_EDGE) == 0
			|| vertex.get_num_edges(edge_type::IN_EDGE) == 0) {
#ifdef DEBUG
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			BOOST_LOG_TRIVIAL(debug)
				<< boost::format("%1% completed vertices") % ret;
#endif
		return;
	}

	vertex_id_t max_id = prog.get_graph().get_max_vertex_id();
	vsize_t part_range
		= ceil(((double) max_id + 1) / graph_conf.get_num_vparts());
	vertex_id_t id_start = part_range * this->get_part_id();
	vertex_id_t id_end = id_start + part_range;

	std::vector<vertex_id_t> selected_neighbors;
	runtime_data_t *data = construct_runtime(prog, vertex,
			prog.get_num_edges(get_id()), selected_neighbors, id_start, id_end);
	if (data == NULL)
		return;

	// We have to set runtime data before calling request_partial_vertices.
	// It's possible that the request to a partial vertex can be completed
	// immediately and run_on_neighbor is called in request_partial_vertices.
	// TODO Maybe I should avoid that.
	assert(local_value.get_num_triangles() == 0);
	local_value.set_runtime_data(data);
	std::vector<directed_vertex_request> reqs;
	for (size_t i = 0; i < selected_neighbors.size(); i++) {
		vertex_id_t id = selected_neighbors[i];
		assert(id >= id_start && id < id_end);
		reqs.push_back(directed_vertex_request(id, neigh_edge_type));
	}
	assert(!reqs.empty());
	data->num_required = reqs.size();
	request_partial_vertices(reqs.data(), reqs.size());
}

void part_directed_triangle_vertex::run_on_neighbor(vertex_program &prog,
		const page_vertex &vertex)
{
	assert(local_value.has_runtime_data());
	runtime_data_t *data = local_value.get_runtime_data();
	data->num_joined++;
	size_t ret = count_triangles(local_value.get_runtime_data(), vertex,
			this->get_id());
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
#ifdef DEBUG
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			BOOST_LOG_TRIVIAL(debug)
				<< boost::format("%1% completed vertices") % ret;
#endif

		// Inform all neighbors in the in-edges.
		for (size_t i = 0; i < data->triangles.size(); i++) {
			// Inform the neighbor if they share triangles.
			if (data->triangles[i] > 0) {
				count_msg msg(data->triangles[i]);
				prog.send_msg(data->edges[i], msg);
			}
		}
		size_t num_curr_triangles = data->num_triangles;
		delete data;
		local_value.set_num_triangles(num_curr_triangles);
		count_msg msg(num_curr_triangles);
		prog.send_msg(get_id(), msg);
	}
}

}

#include "save_result.h"

namespace fg
{

FG_vector<size_t>::ptr compute_directed_triangles_fast(FG_graph::ptr fg,
		directed_triangle_type type)
{
	graph_index::ptr index = NUMA_graph_index<directed_triangle_vertex,
		part_directed_triangle_vertex>::create(fg->get_graph_header());
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
