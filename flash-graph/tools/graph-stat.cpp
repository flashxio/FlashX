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

#include "graph_config.h"
#include "stat.h"
#include "FGlib.h"

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

	FG_graph::ptr fg = FG_graph::create(graph_file, index_file, conf_file);
	const graph_header &header = fg->get_graph_header();

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
	size_t num_non_empty_vertices = 0;

	printf("There are %ld non-empty vertices\n", num_non_empty_vertices);
	if (directed) {
		FG_vector<vsize_t>::ptr in_degrees = get_degree(fg, IN_EDGE);
		FG_vector<vsize_t>::ptr out_degrees = get_degree(fg, OUT_EDGE);
		size_t tot_in_edges = in_degrees->sum<size_t>();
		size_t tot_out_edges = out_degrees->sum<size_t>();
		assert(tot_in_edges == tot_out_edges);
		vsize_t max_num_in_edges = in_degrees->max();
		vsize_t max_num_out_edges = out_degrees->max();
		log_histogram hist_in_edges = in_degrees->log_hist(POWER_CONST);
		log_histogram hist_out_edges = out_degrees->log_hist(POWER_CONST);

		in_degrees->add_in_place(out_degrees);
		vsize_t max_num_edges = in_degrees->max();
		log_histogram hist_edges = in_degrees->log_hist(POWER_CONST);
		printf("There are %ld edges\n", tot_in_edges);
		printf("max edges of a vertex: %ld, max in-edges: %ld, max out-edges: %ld\n",
				(size_t) max_num_edges, (size_t) max_num_in_edges,
				(size_t) max_num_out_edges);

		printf("edge histogram\n");
		hist_edges.print(stdout);
		printf("in-edges histogram: \n");
		hist_in_edges.print(stdout);
		printf("out-edges histogram: \n");
		hist_out_edges.print(stdout);
	}
	else {
		FG_vector<vsize_t>::ptr degrees = get_degree(fg, BOTH_EDGES);
		size_t tot_edges = degrees->sum<size_t>();
		vsize_t max_num_edges = degrees->max();
		printf("There are %ld edges\n", tot_edges);
		printf("max edges of a vertex: %ld\n", (size_t) max_num_edges);

		log_histogram hist_edges = degrees->log_hist(POWER_CONST);
		printf("edge histogram\n");
		hist_edges.print(stdout);
	}
}
