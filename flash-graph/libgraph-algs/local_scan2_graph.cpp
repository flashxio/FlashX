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

#include <unordered_set>

#include "vertex.h"
#include "FG_vector.h"
#include "FGlib.h"
#include "save_result.h"

using namespace fg;

namespace 
{

enum scan2_stage_t
{
	REQ_NEIGH,
	REQ_NEIGH2,
};

class local_scan2_vertex: public compute_vertex
{
	size_t local_scan2;
	scan2_stage_t stage;
	vsize_t num_fetched;
	// The direct neighbors.
	std::unordered_set<vertex_id_t> *neighbors;
	// The vertices reached in the second hop.
	std::unordered_set<vertex_id_t> *neighbors2;
public:
	local_scan2_vertex(vertex_id_t id): compute_vertex(id) {
		neighbors2 = NULL;
		neighbors = NULL;
		local_scan2 = 0;
		num_fetched = 0;
	}

	size_t get_result() const {
		return local_scan2;
	}

	void run(vertex_program &prog) {
		vertex_id_t id = prog.get_vertex_id(*this);
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex) {
		if (vertex.get_id() == prog.get_vertex_id(*this))
			run_on_itself(prog, vertex);
		else if (stage == scan2_stage_t::REQ_NEIGH)
			run_on_neighbor(prog, vertex);
		else
			run_on_neighbor2(prog, vertex);
	}

	void run_on_itself(vertex_program &prog, const page_vertex &vertex);
	void run_on_neighbor(vertex_program &prog, const page_vertex &vertex);
	void run_on_neighbor2(vertex_program &prog, const page_vertex &vertex);
	void run_on_neighs_neigh_list(vertex_program &prog,
			const page_vertex &vertex, edge_type type);
	void run_on_neighs2_neigh_list(vertex_program &prog,
			const page_vertex &vertex, edge_type type);
};

void local_scan2_vertex::run_on_itself(vertex_program &prog,
		const page_vertex &vertex)
{
	stage = scan2_stage_t::REQ_NEIGH;
	neighbors = new std::unordered_set<vertex_id_t>();
	neighbors->insert(vertex.get_neigh_begin(IN_EDGE),
			vertex.get_neigh_end(IN_EDGE));
	neighbors->insert(vertex.get_neigh_begin(OUT_EDGE),
			vertex.get_neigh_end(OUT_EDGE));
	// We should remove the current vertex from the direct neighbor list.
	// It can happen if the graph has self loops.
	vertex_id_t id = prog.get_vertex_id(*this);
	neighbors->erase(id);

	neighbors2 = new std::unordered_set<vertex_id_t>();
	num_fetched = 0;

	// Here we request all direct neighbors.
	if (neighbors->empty()) {
		delete neighbors;
		delete neighbors2;
		neighbors = NULL;
		neighbors2 = NULL;
	}
	else {
		std::vector<vertex_id_t> neigh_vec(neighbors->begin(), neighbors->end());
		request_vertices(neigh_vec.data(), neigh_vec.size());
	}
}

/*
 * There are two tasks in this method:
 *	construct the set of vertices two hops away from the current vertex;
 *	count the number of edges connected with the direct neighbors.
 */
void local_scan2_vertex::run_on_neighs_neigh_list(vertex_program &prog,
		const page_vertex &vertex, edge_type type)
{
	vertex_id_t curr_id = prog.get_vertex_id(*this);
	auto it = vertex.get_neigh_seq_it(type);
	PAGE_FOREACH(vertex_id_t, id, it) {
		auto it1 = neighbors->find(id);
		// If the neighbor's neighbor is a direct neighbor, we need to count
		// each edge only once.
		if (it1 != neighbors->end()) {
			if (*it1 > vertex.get_id())
				local_scan2++;
		}
		// If the neighbor's neighbor isn't a direct neighbor, it's either
		// the current vertex or the vertices two hops away from the current
		// vertex. Either way, we should add it to the count directly.
		else {
			local_scan2++;
			// However, if the vertex isn't the current vertex, it must be
			// two hops away from the current vertex.
			if (*it1 != curr_id)
				neighbors2->insert(*it1);
		}
	} PAGE_FOREACH_END
}

/*
 * This method processes direct neighbors.
 */
void local_scan2_vertex::run_on_neighbor(vertex_program &prog,
		const page_vertex &vertex)
{
	num_fetched++;
	run_on_neighs_neigh_list(prog, vertex, IN_EDGE);
	run_on_neighs_neigh_list(prog, vertex, OUT_EDGE);

	// We have got all direct neighbors.
	if (num_fetched == neighbors->size()) {
		num_fetched = 0;
		stage = scan2_stage_t::REQ_NEIGH2;

		// neighbors2 contain neighbors and the current vertex.
		// We need to remove them from neighbors2.
		BOOST_FOREACH(vertex_id_t id, *neighbors)
			neighbors2->erase(id);
		vertex_id_t curr_id = prog.get_vertex_id(*this);
		neighbors2->erase(curr_id);

		if (neighbors2->empty()) {
			delete neighbors2;
			neighbors2 = NULL;
		}
		else {
			std::vector<vertex_id_t> neigh_vec(neighbors2->begin(),
					neighbors2->end());
			request_vertices(neigh_vec.data(), neigh_vec.size());
		}
		delete neighbors;
		neighbors = NULL;
	}
}

void local_scan2_vertex::run_on_neighs2_neigh_list(vertex_program &prog,
		const page_vertex &vertex, edge_type type)
{
	auto it = vertex.get_neigh_seq_it(type);
	PAGE_FOREACH(vertex_id_t, id, it) {
		auto it1 = neighbors2->find(id);
		// If its neighbor is also two hops away from the current vertex,
		// we count each edge only once.
		if (it1 != neighbors2->end()) {
			if (*it1 > vertex.get_id())
				local_scan2++;
		}
	} PAGE_FOREACH_END
}

void local_scan2_vertex::run_on_neighbor2(vertex_program &prog,
		const page_vertex &vertex)
{
	num_fetched++;
	run_on_neighs2_neigh_list(prog, vertex, IN_EDGE);
	run_on_neighs2_neigh_list(prog, vertex, OUT_EDGE);

	if (num_fetched == neighbors2->size()) {
		delete neighbors2;
		neighbors2 = NULL;
	}
}

}

namespace fg
{

FG_vector<size_t>::ptr compute_local_scan2(FG_graph::ptr fg)
{
	graph_index::ptr index = NUMA_graph_index<local_scan2_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);

	BOOST_LOG_TRIVIAL(info) << "local scan starts";
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
		<< boost::format("It takes %1% seconds to compute all local scan")
		% time_diff(start, end);

	FG_vector<size_t>::ptr vec = FG_vector<size_t>::create(graph);
	graph->query_on_all(vertex_query::ptr(
				new save_query<size_t, local_scan2_vertex>(vec)));
	return vec;
}

}
