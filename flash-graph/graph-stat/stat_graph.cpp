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

#include <vector>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"
#include "stat.h"

class stat_vertex: public compute_vertex
{
	int num_in_edges;
	int num_out_edges;
public:
	stat_vertex(vertex_id_t id): compute_vertex(id) {
		num_in_edges = 0;
		num_out_edges = 0;
	}

	virtual void run(vertex_program &prog) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex) {
		num_in_edges = vertex.get_num_edges(edge_type::IN_EDGE);
		num_out_edges = vertex.get_num_edges(edge_type::OUT_EDGE);
	}

	virtual void run_on_message(vertex_program &, const vertex_message &msg) {
	}

	int get_num_edges(edge_type type) const {
		switch (type) {
			case edge_type::IN_EDGE:
				return num_in_edges;
			case edge_type::OUT_EDGE:
				return num_out_edges;
			case edge_type::BOTH_EDGES:
				return num_in_edges + num_out_edges;
			default:
				assert(0);
		}
	}
};

const int POWER_CONST = 10;

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "stat conf_file graph_file index_file\n");
		graph_conf.print_help();
		params.print_help();
		exit(-1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];

	config_map configs(conf_file);
	graph_conf.init(configs);

	graph_index::ptr index = NUMA_graph_index<stat_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	const graph_header &header = graph->get_graph_header();
	printf("start the graph algorithm\n");
	graph->start_all();
	graph->wait4complete();

	std::string graph_type_str;
	switch(header.get_graph_type()) {
		case graph_type::DIRECTED:
			graph_type_str = "directed";
			break;
		case graph_type::UNDIRECTED:
			graph_type_str = "undirected";
			break;
		case graph_type::TS_DIRECTED:
			graph_type_str = "time-series directed";
			break;
		case graph_type::TS_UNDIRECTED:
			graph_type_str = "time-series undirected";
			break;
		default:
			assert(0);
	}
	printf("The graph type: %s\n", graph_type_str.c_str());
	printf("There are %ld vertices and %ld edges\n",
			header.get_num_vertices(), header.get_num_edges());

	bool directed = header.is_directed_graph();
	graph_index::const_iterator it = index->begin();
	graph_index::const_iterator end_it = index->end();
	vertex_id_t min_id = INT_MAX;
	vertex_id_t max_id = 0;
	size_t max_num_edges = 0;
	size_t max_num_in_edges = 0;
	size_t max_num_out_edges = 0;
	size_t tot_edges = 0;
	size_t tot_in_edges = 0;
	size_t tot_out_edges = 0;
	size_t num_non_empty_vertices = 0;
	for (; it != end_it; ++it) {
		const stat_vertex &v = (const stat_vertex &) *it;
		if (v.get_id() > max_id)
			max_id = v.get_id();
		if (v.get_id() < min_id)
			min_id = v.get_id();
		if (v.get_num_edges(edge_type::BOTH_EDGES) > 0)
			num_non_empty_vertices++;
		if (directed) {
			tot_in_edges += v.get_num_edges(edge_type::IN_EDGE);
			tot_out_edges += v.get_num_edges(edge_type::OUT_EDGE);
			if (max_num_edges < (size_t) v.get_num_edges(edge_type::BOTH_EDGES))
				max_num_edges = v.get_num_edges(edge_type::BOTH_EDGES);
			if (max_num_in_edges < (size_t) v.get_num_edges(edge_type::IN_EDGE))
				max_num_in_edges = v.get_num_edges(edge_type::IN_EDGE);
			if (max_num_out_edges < (size_t) v.get_num_edges(edge_type::OUT_EDGE))
				max_num_out_edges = v.get_num_edges(edge_type::OUT_EDGE);
		}
		else {
			tot_edges += v.get_num_edges(edge_type::IN_EDGE);
			assert(v.get_num_edges(edge_type::IN_EDGE)
					== v.get_num_edges(edge_type::OUT_EDGE));
			if (max_num_edges < (size_t) v.get_num_edges(edge_type::IN_EDGE))
				max_num_edges = v.get_num_edges(edge_type::IN_EDGE);
		}
	}

	printf("min id: %ld, max id: %ld\n", (long) min_id, (long) max_id);
	printf("There are %ld non-empty vertices\n", num_non_empty_vertices);
	if (directed) {
		assert(tot_in_edges == tot_out_edges);
		printf("There are %ld edges\n", tot_in_edges);
		printf("max edges of a vertex: %ld, max in-edges: %ld, max out-edges: %ld\n",
				max_num_edges, max_num_in_edges, max_num_out_edges);
	}
	else {
		printf("There are %ld edges\n", tot_edges);
		printf("max edges of a vertex: %ld\n", max_num_edges);
	}

	it = index->begin();
	end_it = index->end();
	int num_buckets = ceil(log(max_num_edges) / log(POWER_CONST));
	log_histogram hist_edges(num_buckets);
	num_buckets = ceil(log(max_num_in_edges) / log(POWER_CONST));
	log_histogram hist_in_edges(num_buckets);
	num_buckets = ceil(log(max_num_out_edges) / log(POWER_CONST));
	log_histogram hist_out_edges(num_buckets);

	for (; it != end_it; ++it) {
		const stat_vertex &v = (const stat_vertex &) *it;
		if (directed) {
			hist_edges.add_value(v.get_num_edges(edge_type::BOTH_EDGES));
			hist_in_edges.add_value(v.get_num_edges(edge_type::IN_EDGE));
			hist_out_edges.add_value(v.get_num_edges(edge_type::OUT_EDGE));
		}
		else
			hist_edges.add_value(v.get_num_edges(edge_type::IN_EDGE));
	}
	printf("edge histogram\n");
	hist_edges.print(stdout);
	if (directed) {
		printf("in-edges histogram: \n");
		hist_in_edges.print(stdout);
		printf("out-edges histogram: \n");
		hist_out_edges.print(stdout);
	}
}
