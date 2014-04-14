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

class bitmap
{
	uint32_t v;
public:
	bitmap() {
		v = 0;
	}

	void set(int bit) {
		v |= 0x1U << bit;
	}

	bitmap &operator|=(const bitmap &map) {
		v |= map.v;
		return *this;
	}

	bool operator!=(const bitmap &map) const {
		return v != map.v;
	}

	uint32_t get_value() const {
		return v;
	}
};

class diameter_message: public vertex_message
{
	bitmap bfs_map;
public:
	diameter_message(bitmap bfs_map): vertex_message(
			sizeof(diameter_message), true) {
		this->bfs_map = bfs_map;
	}

	const bitmap &get_bfs_map() const {
		return bfs_map;
	}
};

class diameter_vertex: public compute_vertex
{
	bitmap bfs_map;
	bool updated;
public:
	diameter_vertex() {
		updated = false;
	}

	diameter_vertex(vertex_id_t id,
			const vertex_index *index): compute_vertex(id, index) {
		updated = false;
	}

	void init(int bfs_id) {
		bfs_map.set(bfs_id);
		updated = true;
		printf("v%u gets %x\n", get_id(), bfs_map.get_value());
	}

	void run(graph_engine &graph) {
		if (updated) {
			updated = false;
			vertex_id_t id = get_id();
			request_vertices(&id, 1);
		}
	}

	void run(graph_engine &graph, const page_vertex &vertex);

	void run_on_message(graph_engine &, const vertex_message &msg) {
		const diameter_message &dmsg = (const diameter_message &) msg;
		bitmap old_map = bfs_map;
		bfs_map |= dmsg.get_bfs_map();
		if (bfs_map != old_map)
			updated = true;
	}
};

void diameter_vertex::run(graph_engine &graph, const page_vertex &vertex)
{
	int num_dests = vertex.get_num_edges(BOTH_EDGES);
	if (num_dests == 0)
		return;

	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	edge_seq_iterator it = vertex.get_neigh_seq_it(BOTH_EDGES, 0, num_dests);
	diameter_message msg(bfs_map);
	graph.multicast_msg(it, msg);
}

class diameter_initiator: public vertex_initiator
{
	std::unordered_map<vertex_id_t, int> start_vertices;
public:
	diameter_initiator(std::vector<vertex_id_t> &vertices) {
		for (size_t i = 0; i < vertices.size(); i++) {
			start_vertices.insert(std::pair<vertex_id_t, int>(vertices[i], i));
		}
	}

	void init(compute_vertex &v) {
		diameter_vertex &dv = (diameter_vertex &) v;
		std::unordered_map<vertex_id_t, int>::const_iterator it
			= start_vertices.find(v.get_id());
		assert(it != start_vertices.end());
		dv.init(it->second);
	}
};

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

void print_usage()
{
	fprintf(stderr,
			"diameter_estimator [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-p: preload the graph\n");
	fprintf(stderr, "-n: the number of paprallel BFS\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_opts = 0;
	bool preload = false;
	size_t num_bfs = 1;
	while ((opt = getopt(argc, argv, "c:pn:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 'p':
				preload = true;
				break;
			case 'n':
				num_bfs = atoi(optarg);
				assert(num_bfs < sizeof(bitmap) * 8);
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

	graph_index *index = NUMA_graph_index<diameter_vertex>::create(index_file,
			graph_conf.get_num_threads(), params.get_num_nodes());

	graph_engine *graph = graph_engine::create(graph_conf.get_num_threads(),
			params.get_num_nodes(), graph_file, index);
	if (preload)
		graph->preload_graph();
	printf("BFS starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	struct timeval start, end;
	gettimeofday(&start, NULL);
	std::vector<vertex_id_t> start_vertices;
	while (start_vertices.size() < num_bfs) {
		vertex_id_t id = random() % graph->get_max_vertex_id();
		// We should skip the empty vertices.
		if (graph->get_vertex_edges(id) == 0)
			continue;

		start_vertices.push_back(id);
	}
	graph->start(start_vertices.data(), start_vertices.size(),
			vertex_initiator::ptr(new diameter_initiator(start_vertices)));
	graph->wait4complete();
	gettimeofday(&end, NULL);

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph_engine::destroy(graph);
	destroy_io_system();
	printf("It takes %f seconds\n", time_diff(start, end));
}
