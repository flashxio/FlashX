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
#include "ts_graph.h"

const int POWER_CONST = 10;

void print_usage()
{
	fprintf(stderr, "stat [options] conf_file graph_file index_file\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	int num_opts = 0;

	std::string start_time_str;
	std::string time_unit_str;
	time_t time_interval = 1;
	time_t start_time = 0;
	while ((opt = getopt(argc, argv, "t:l:u:")) != -1) {
		num_opts++;
		switch (opt) {
			case 't':
				start_time_str = optarg;
				num_opts++;
				break;
			case 'l':
				time_interval = atol(optarg);
				num_opts++;
				break;
			case 'u':
				time_unit_str = optarg;
				num_opts++;
				break;
			default:
				print_usage();
		}
	}
	argv += 1 + num_opts;
	argc -= 1 + num_opts;

	if (argc < 3) {
		print_usage();
		exit(-1);
	}

	std::string conf_file = argv[0];
	std::string graph_file = argv[1];
	std::string index_file = argv[2];

	config_map configs(conf_file);
	graph_conf.init(configs);

	bool on_ts = false;
	if (!start_time_str.empty()) {
		on_ts = true;
		if (!time_unit_str.empty()) {
			if (time_unit_str == "hour")
				time_interval *= HOUR_SECS;
			else if (time_unit_str == "day")
				time_interval *= DAY_SECS;
			else if (time_unit_str == "month")
				time_interval *= MONTH_SECS;
			else
				fprintf(stderr, "a wrong time unit: %s\n", time_unit_str.c_str());
		}

		start_time = conv_str_to_time(start_time_str);
		printf("start time: %ld, interval: %ld\n", start_time, time_interval);
	}

	FG_graph::ptr fg = FG_graph::create(graph_file, index_file, conf_file);
	graph_header header = get_graph_header(fg);
	header = get_graph_header(fg);

	std::string graph_type_str;
	switch(header.get_graph_type()) {
		case graph_type::DIRECTED:
			graph_type_str = "directed";
			break;
		case graph_type::UNDIRECTED:
			graph_type_str = "undirected";
			break;
		default:
			assert(0);
	}
	printf("The graph type: %s\n", graph_type_str.c_str());
	printf("There are %ld vertices and %ld edges\n",
			header.get_num_vertices(), header.get_num_edges());

	struct non_empty_func {
		size_t operator()(vsize_t v) {
			return !!v;
		}
	};

	bool directed = header.is_directed_graph();
	if (directed) {
		FG_vector<vsize_t>::ptr in_degrees;
		FG_vector<vsize_t>::ptr out_degrees;
		if (on_ts) {
			in_degrees = get_ts_degree(fg, IN_EDGE, start_time, time_interval);
			out_degrees = get_ts_degree(fg, OUT_EDGE, start_time, time_interval);
		}
		else {
			in_degrees = get_degree(fg, IN_EDGE);
			out_degrees = get_degree(fg, OUT_EDGE);
		}
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
		printf("There are %ld non-empty vertices\n",
				in_degrees->aggregate<non_empty_func, size_t>(non_empty_func()));

		printf("edge histogram\n");
		hist_edges.print(stdout);
		printf("in-edges histogram: \n");
		hist_in_edges.print(stdout);
		printf("out-edges histogram: \n");
		hist_out_edges.print(stdout);
	}
	else {
		FG_vector<vsize_t>::ptr degrees;
		if (on_ts)
			degrees = get_ts_degree(fg, BOTH_EDGES, start_time, time_interval);
		else
			degrees = get_degree(fg, BOTH_EDGES);
		size_t tot_edges = degrees->sum<size_t>();
		vsize_t max_num_edges = degrees->max();
		printf("There are %ld edges\n", tot_edges);
		printf("max edges of a vertex: %ld\n", (size_t) max_num_edges);
		printf("There are %ld non-empty vertices\n",
				degrees->aggregate<non_empty_func, size_t>(non_empty_func()));

		log_histogram hist_edges = degrees->log_hist(POWER_CONST);
		printf("edge histogram\n");
		hist_edges.print(stdout);
	}
}
