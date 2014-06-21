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
#include <google/profiler.h>
#endif

#include <vector>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"

class bfs_vertex: public compute_vertex
{
	bool visited;
public:
	bfs_vertex(vertex_id_t id): compute_vertex(id) {
		visited = false;
	}

	bool has_visited() const {
		return visited;
	}

	void run(vertex_program &prog) {
		if (!has_visited()) {
			vertex_id_t id = get_id();
			request_vertices(&id, 1);
		}
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}
};

void bfs_vertex::run(vertex_program &prog, const page_vertex &vertex)
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

class count_vertex_query: public vertex_query
{
	size_t num_visited;
public:
	count_vertex_query() {
		num_visited = 0;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		bfs_vertex &bfs_v = (bfs_vertex &) v;
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

void int_handler(int sig_num)
{
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	exit(0);
}

void print_usage()
{
	fprintf(stderr,
			"bfs [options] conf_file graph_file index_file start_vertex\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_opts = 0;
	while ((opt = getopt(argc, argv, "c:b")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			default:
				print_usage();
		}
	}
	argv += 1 + num_opts;
	argc -= 1 + num_opts;

	if (argc < 4) {
		print_usage();
		exit(-1);
	}

	std::string conf_file = argv[0];
	std::string graph_file = argv[1];
	std::string index_file = argv[2];
	vertex_id_t start_vertex = atoi(argv[3]);

	config_map configs(conf_file);
	configs.add_options(confs);

	signal(SIGINT, int_handler);

	graph_index::ptr index = NUMA_graph_index<bfs_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	printf("BFS starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start(&start_vertex, 1);
	graph->wait4complete();
	gettimeofday(&end, NULL);

	vertex_query::ptr cvq(new count_vertex_query());
	graph->query_on_all(cvq);
	size_t num_visited = ((count_vertex_query *) cvq.get())->get_num_visited();

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	printf("BFS from vertex %ld visits %ld vertices. It takes %f seconds\n",
			(unsigned long) start_vertex, num_visited, time_diff(start, end));
}
