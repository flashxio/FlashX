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

#include <signal.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include <vector>

#include "graph_engine.h"
#include "graph_config.h"
#include "FGlib.h"

using namespace safs;
using namespace fg;

namespace
{

edge_type traverse_edge = edge_type::OUT_EDGE;

/*
 * Vertex program for BFS on a directed graph.
 */
class bfs_dvertex: public compute_directed_vertex
{
	bool visited;
public:
	bfs_dvertex(vertex_id_t id): compute_directed_vertex(id) {
		visited = false;
	}

	bool has_visited() const {
		return visited;
	}

	void set_visited(bool visited) {
		this->visited = visited;
	}

	void run(vertex_program &prog) {
		if (!has_visited()) {
			directed_vertex_request req(prog.get_vertex_id(*this),
					traverse_edge);
			request_partial_vertices(&req, 1);
		}
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}
};

void bfs_dvertex::run(vertex_program &prog, const page_vertex &vertex)
{
	assert(!has_visited());
	set_visited(true);

	int num_dests = vertex.get_num_edges(traverse_edge);
	if (num_dests == 0)
		return;

	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	if (traverse_edge == BOTH_EDGES) {
		edge_seq_iterator it = vertex.get_neigh_seq_it(IN_EDGE, 0,
				num_dests);
		prog.activate_vertices(it);
		it = vertex.get_neigh_seq_it(OUT_EDGE, 0, num_dests);
		prog.activate_vertices(it);
	}
	else {
		edge_seq_iterator it = vertex.get_neigh_seq_it(traverse_edge, 0,
				num_dests);
		prog.activate_vertices(it);
	}
}

/*
 * Vertex program for BFS on an undirected graph.
 */
class bfs_uvertex: public compute_vertex
{
	bool visited;
public:
	bfs_uvertex(vertex_id_t id): compute_vertex(id) {
		visited = false;
	}

	bool has_visited() const {
		return visited;
	}

	void run(vertex_program &prog) {
		if (!has_visited()) {
			vertex_id_t id = prog.get_vertex_id(*this);
			request_vertices(&id, 1);
		}
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}
};

void bfs_uvertex::run(vertex_program &prog, const page_vertex &vertex)
{
	assert(!has_visited());
	visited = true;

	int num_dests = vertex.get_num_edges(edge_type::BOTH_EDGES);
	if (num_dests == 0)
		return;

	// We need to add the neighbors of the vertex to the queue of
	// the next level.
#ifdef USE_ARRAY
	stack_array<vertex_id_t, 1024> neighs(num_dests);
	vertex.read_edges(edge_type::BOTH_EDGES, neighs.data(), num_dests);
	prog.activate_vertices(neighs.data(), num_dests);
#else
	edge_seq_iterator it = vertex.get_neigh_seq_it(edge_type::BOTH_EDGES, 0, num_dests);
	prog.activate_vertices(it);
#endif
}

template<class vertex_type>
class count_vertex_query: public vertex_query
{
	size_t num_visited;
public:
	count_vertex_query() {
		num_visited = 0;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		vertex_type &bfs_v = (vertex_type &) v;
		if (bfs_v.has_visited())
			num_visited++;
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		count_vertex_query *cvq = (count_vertex_query *) q.get();
		num_visited += cvq->num_visited;
	}

	virtual ptr clone() {
		return vertex_query::ptr(new count_vertex_query());
	}

	size_t get_num_visited() const {
		return num_visited;
	}
};

}

size_t bfs(FG_graph::ptr fg, vertex_id_t start_vertex, edge_type traverse_e)
{
	bool directed = fg->get_graph_header().is_directed_graph();
	graph_index::ptr index;
	if (directed)
		index = NUMA_graph_index<bfs_dvertex>::create(fg->get_graph_header());
	else
		index = NUMA_graph_index<bfs_uvertex>::create(fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);

	traverse_edge = traverse_e;
	printf("BFS starts\n");
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	graph->start(&start_vertex, 1);
	graph->wait4complete();

	size_t num_visited;
	vertex_query::ptr cvq;
	if (directed) {
		cvq = vertex_query::ptr(new count_vertex_query<bfs_dvertex>());
		graph->query_on_all(cvq);
		num_visited = ((count_vertex_query<bfs_dvertex> *) cvq.get())->get_num_visited();
	}
	else {
		cvq = vertex_query::ptr(new count_vertex_query<bfs_uvertex>());
		graph->query_on_all(cvq);
		num_visited = ((count_vertex_query<bfs_uvertex> *) cvq.get())->get_num_visited();
	}

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	return num_visited;
}
