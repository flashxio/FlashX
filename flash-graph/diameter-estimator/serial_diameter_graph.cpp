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

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"

edge_type traverse_edge = edge_type::OUT_EDGE;
int start_level;

class bfs_vertex: public compute_directed_vertex
{
	bool visited;
	short dist;
public:
	bfs_vertex() {
		visited = false;
		dist = 0;
	}

	bfs_vertex(vertex_id_t id,
			const vertex_index &index): compute_directed_vertex(id, index) {
		visited = false;
		dist = 0;
	}

	void reset() {
		visited = false;
		dist = 0;
	}

	bool has_visited() const {
		return visited;
	}

	short get_dist() const {
		return dist;
	}

	void run(vertex_program &prog) {
		if (!has_visited()) {
			visited = true;
			dist = prog.get_graph().get_curr_level() - start_level;
			directed_vertex_request req(get_id(), traverse_edge);
			request_partial_vertices(&req, 1);
		}
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}
};

void bfs_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	int num_dests = vertex.get_num_edges(traverse_edge);
	if (num_dests == 0)
		return;

	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	edge_seq_iterator it = vertex.get_neigh_seq_it(traverse_edge, 0, num_dests);
	prog.activate_vertices(it);
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
			"bfs [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-p: preload the graph\n");
	fprintf(stderr, "-b: traverse with both in-edges and out-edges\n");
	fprintf(stderr, "-n: the number of BFS\n");
	fprintf(stderr, "-o: the output file\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_bfs = 1;
	int num_opts = 0;
	bool preload = false;
	std::string output_file;
	while ((opt = getopt(argc, argv, "c:pbo:n:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 'p':
				preload = true;
				break;
			case 'b':
				traverse_edge = edge_type::BOTH_EDGES;
				break;
			case 'n':
				num_bfs = atoi(optarg);
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

	signal(SIGINT, int_handler);

	graph_index::ptr index = NUMA_graph_index<bfs_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	if (preload)
		graph->preload_graph();
	printf("BFS starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	std::vector<vertex_id_t> start_vertices;
	while (start_vertices.size() < (size_t) num_bfs) {
		vertex_id_t id = random() % graph->get_max_vertex_id();
		// We should skip the empty vertices.
		if (graph->get_vertex_edges(id) == 0)
			continue;

		start_vertices.push_back(id);
	}

	std::vector<short> dist_vec(graph->get_num_vertices());
	for (size_t i = 0; i < start_vertices.size(); i++) {
		struct timeval start, end;
		gettimeofday(&start, NULL);
		start_level = graph->get_curr_level();
		printf("The starting level is %d\n", start_level);
		vertex_id_t start_vertex = start_vertices[i];
		graph->start(&start_vertex, 1);
		graph->wait4complete();
		gettimeofday(&end, NULL);

		// Clean up all vertices.
		graph_index::const_iterator it = index->begin();
		graph_index::const_iterator end_it = index->end();
		size_t i = 0;
		for (; it != end_it; ++it) {
			bfs_vertex &v = (bfs_vertex &) *it;
			dist_vec[i] = max(dist_vec[i], v.get_dist());
			i++;
#if 0
			if (v.get_dist(0) != USHRT_MAX)
				fprintf(stderr, "v%u: %d\n", v.get_id(), v.get_dist(0));
#endif
			v.reset();
		}
		printf("BFS on v%u takes %f seconds\n", start_vertex,
				time_diff(start, end));
	}

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();

	if (!output_file.empty()) {
		FILE *f = fopen(output_file.c_str(), "w");
		assert(f);
		BOOST_FOREACH(short dist, dist_vec) {
			fprintf(f, "%d\n", dist);
		}
		fclose(f);
	}
}
