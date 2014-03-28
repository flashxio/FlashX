/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SA-GraphLib.
 *
 * SA-GraphLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SA-GraphLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SA-GraphLib.  If not, see <http://www.gnu.org/licenses/>.
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
	vertex_id_t component_id;
	vertex_id_t neigh_min;
	int num_edges;
public:
	wcc_vertex() {
		component_id = UINT_MAX;
		neigh_min = UINT_MAX;
		num_edges = 0;
	}

	wcc_vertex(vertex_id_t id, const vertex_index *index1): compute_vertex(
			id, index1) {
		component_id = UINT_MAX;
		neigh_min = id;
		const directed_vertex_index *index = (const directed_vertex_index *) index1;
		num_edges = (index->get_num_in_edges(id) + index->get_num_out_edges(id));
		if (num_edges == 0) {
			component_id = id;
			neigh_min = id;
		}
	}

	bool get_num_edges() const {
		return num_edges;
	}

	bool belong2component() const {
		return component_id != UINT_MAX;
	}

	vertex_id_t get_component_id() const {
		return component_id;
	}

	void run(graph_engine &graph) {
		if (neigh_min < component_id) {
			component_id = neigh_min;
			vertex_id_t id = get_id();
			request_vertices(&id, 1);
		}
	}

	void run(graph_engine &graph, const page_vertex &vertex);

	void run_on_messages(graph_engine &,
			const vertex_message *msgs[], int num) {
		for (int i = 0; i < num; i++) {
			component_message *msg = (component_message *) msgs[i];
			neigh_min = min(neigh_min, msg->get_id());
		}
	}
};

void wcc_vertex::run(graph_engine &graph, const page_vertex &vertex)
{
	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	page_byte_array::const_iterator<vertex_id_t> end_it
		= vertex.get_neigh_end(BOTH_EDGES);
	stack_array<vertex_id_t, 1024> dest_buf(vertex.get_num_edges(BOTH_EDGES));
	int num_dests = 0;
	for (page_byte_array::const_iterator<vertex_id_t> it
			= vertex.get_neigh_begin(BOTH_EDGES); it != end_it; ++it) {
		vertex_id_t id = *it;
		dest_buf[num_dests++] = id;
	}
	component_message msg(component_id);
	graph.multicast_msg(dest_buf.data(), num_dests, msg);
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
	while ((opt = getopt(argc, argv, "c:s:o:")) != -1) {
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
	graph_conf.init(configs);
	graph_conf.print();

	signal(SIGINT, int_handler);
	init_io_system(configs);

	graph_index *index = NUMA_graph_index<wcc_vertex>::create(index_file,
			graph_conf.get_num_threads(), params.get_num_nodes());
	graph_engine *graph = graph_engine::create(graph_conf.get_num_threads(),
			params.get_num_nodes(), graph_file, index);
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
	graph_engine::destroy(graph);
	destroy_io_system();

	typedef std::unordered_map<vertex_id_t, size_t> comp_map_t;
	comp_map_t comp_counts;
	graph_index::const_iterator it = index->begin();
	graph_index::const_iterator end_it = index->end();
	for (; it != end_it; ++it) {
		const wcc_vertex &v = (const wcc_vertex &) *it;
		if (v.get_num_edges() == 0)
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
	printf("There are %ld components\n", comp_counts.size());

	FILE *f = stdout;
	if (!output_file.empty()) {
		f = fopen(output_file.c_str(), "w");
		assert(f);
	}
	BOOST_FOREACH(comp_map_t::value_type &p, comp_counts) {
		if (p.second >= min_comp_size)
			fprintf(f, "component %u: %ld\n", p.first, p.second);
	}
	if (!output_file.empty()) {
		fclose(f);
	}
}
