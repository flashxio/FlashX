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

#include <vector>
#include <unordered_map>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"
#include "FG_vector.h"
#include "FGlib.h"

namespace {

enum wcc_stage_t
{
	FIND_COMPONENTS,
	REMOVE_EMPTY,
} wcc_stage;

class component_message: public vertex_message
{
	int id;
public:
	component_message(vertex_id_t id): vertex_message(
			sizeof(component_message), true) {
		this->id = id;
	}

	vertex_id_t get_id() const {
		return id;
	}
};

class wcc_vertex: public compute_vertex
{
	bool updated;
	bool empty;
	vertex_id_t component_id;
public:
	wcc_vertex(vertex_id_t id): compute_vertex(id) {
		component_id = id;
		updated = true;
		empty = false;
	}

	bool is_empty(graph_engine &graph) const {
		return empty;
	}

	bool belong2component() const {
		return component_id != UINT_MAX;
	}

	vertex_id_t get_component_id() const {
		return component_id;
	}

	void run(vertex_program &prog) {
		vertex_id_t id = get_id();
		if (wcc_stage == wcc_stage_t::FIND_COMPONENTS) {
			if (updated) {
				request_vertices(&id, 1);
				updated = false;
			}
		}
		else {
			request_vertex_headers(&id, 1);
		}
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &, const vertex_message &msg1) {
		component_message &msg = (component_message &) msg1;
		if (msg.get_id() < component_id) {
			updated = true;
			component_id = msg.get_id();
		}
	}

	void run_on_vertex_header(vertex_program &, const vertex_header &header) {
		assert(get_id() == header.get_id());
		empty = (header.get_num_edges() == 0);
	}

	vertex_id_t get_result() const {
		if (!empty)
			return get_component_id();
		else
			return INVALID_VERTEX_ID;
	}
};

void wcc_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	int num_dests = vertex.get_num_edges(BOTH_EDGES);
	edge_seq_iterator it = vertex.get_neigh_seq_it(BOTH_EDGES, 0, num_dests);
	component_message msg(component_id);
	prog.multicast_msg(it, msg);
}

}

#include "save_result.h"
FG_vector<vertex_id_t>::ptr compute_wcc(FG_graph::ptr fg)
{
	graph_index::ptr index = NUMA_graph_index<wcc_vertex>::create(
			fg->get_index_file());
	graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
			index, fg->get_configs());
	printf("weakly connected components starts\n");
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	struct timeval start, end;
	gettimeofday(&start, NULL);
	wcc_stage = wcc_stage_t::FIND_COMPONENTS;
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("WCC takes %f seconds\n", time_diff(start, end));

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif

	wcc_stage = wcc_stage_t::REMOVE_EMPTY;
	graph->start_all();
	graph->wait4complete();
	FG_vector<vertex_id_t>::ptr vec = FG_vector<vertex_id_t>::create(graph);
	graph->query_on_all(vertex_query::ptr(
				new save_query<vertex_id_t, wcc_vertex>(vec)));
	return vec;
}
