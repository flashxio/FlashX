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
#include <google/profiler.h>

#include <vector>
#include <unordered_map>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"
#include "stat.h"

atomic_number<long> num_visits;

enum wcc_stage_t
{
	REMOVE_EMPTY,
	FIND_COMPONENTS,
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
	vertex_id_t component_id;
public:
	wcc_vertex() {
		component_id = UINT_MAX;
		updated = true;
	}

	wcc_vertex(vertex_id_t id, const vertex_index &index1): compute_vertex(
			id, index1) {
		component_id = id;
		updated = true;
	}

	bool is_empty(graph_engine &graph) const {
		return graph.get_vertex_edges(get_id()) == 0;
	}

	bool belong2component() const {
		return component_id != UINT_MAX;
	}

	vertex_id_t get_component_id() const {
		return component_id;
	}

	void run(vertex_program &prog) {
		if (updated) {
			vertex_id_t id = get_id();
			request_vertices(&id, 1);
			updated = false;
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

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

void print_usage()
{
	fprintf(stderr,
			"wcc [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-s size: the output min component size\n");
	fprintf(stderr, "-o file: output the component size to the file\n");
	fprintf(stderr, "-p: preload the graph\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	std::string output_file;
	size_t min_comp_size = 0;
	int num_opts = 0;
	bool preload = false;
	while ((opt = getopt(argc, argv, "c:s:o:p")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 's':
				min_comp_size = atoi(optarg);
				num_opts++;
				break;
			case 'o':
				output_file = optarg;
				num_opts++;
				break;
			case 'p':
				preload = true;
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
	configs.add_options(confs);

	signal(SIGINT, int_handler);

	graph_index::ptr index = NUMA_graph_index<wcc_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	if (preload)
		graph->preload_graph();
	printf("weakly connected components starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("It takes %f seconds\n", time_diff(start, end));

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();

	typedef std::unordered_map<vertex_id_t, size_t> comp_map_t;
	comp_map_t comp_counts;
	graph_index::const_iterator it = index->begin();
	graph_index::const_iterator end_it = index->end();
	for (; it != end_it; ++it) {
		const wcc_vertex &v = (const wcc_vertex &) *it;
		if (v.is_empty(*graph))
			continue;

		comp_map_t::iterator map_it = comp_counts.find(v.get_component_id());
		if (map_it == comp_counts.end()) {
			comp_counts.insert(std::pair<vertex_id_t, size_t>(
						v.get_component_id(), 1));
		}
		else {
			map_it->second++;
		}
	}
	size_t max_comp_size = 0;
	BOOST_FOREACH(comp_map_t::value_type v, comp_counts) {
		if (max_comp_size < v.second)
			max_comp_size = v.second;
	}
	printf("There are %ld components, and largest comp has %ld vertices\n",
			comp_counts.size(), max_comp_size);

	if (!output_file.empty()) {
		FILE *f = fopen(output_file.c_str(), "w");
		assert(f);
		std::unordered_map<vertex_id_t, log_histogram> comp_hist_map;
		BOOST_FOREACH(comp_map_t::value_type &p, comp_counts) {
			if (p.second >= min_comp_size) {
				comp_hist_map.insert(std::pair<vertex_id_t, log_histogram>(p.first,
							log_histogram(10)));
				fprintf(f, "component %u: %ld\n", p.first, p.second);
			}
		}

		it = index->begin();
		end_it = index->end();
		for (; it != end_it; ++it) {
			const wcc_vertex &v = (const wcc_vertex &) *it;
			if (v.is_empty(*graph))
				continue;

			std::unordered_map<vertex_id_t, log_histogram>::iterator map_it
				= comp_hist_map.find(v.get_component_id());
			if (map_it != comp_hist_map.end())
				map_it->second.add_value(graph->get_vertex_edges(v.get_id()));
		}
		for (std::unordered_map<vertex_id_t, log_histogram>::iterator it
				= comp_hist_map.begin(); it != comp_hist_map.end(); it++) {
			fprintf(f, "comp %u:\n", it->first);
			it->second.print(f);
		}

		fclose(f);
	}
}
