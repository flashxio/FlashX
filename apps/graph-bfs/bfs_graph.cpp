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

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"

edge_type traverse_edge = edge_type::OUT_EDGE;

class bfs_vertex: public compute_directed_vertex
{
	enum {
		VISITED,
	};

	atomic_flags<int> flags;
public:
	bfs_vertex() {
	}

	bfs_vertex(vertex_id_t id,
			const vertex_index &index): compute_directed_vertex(id, index) {
	}

	bool has_visited() const {
		return flags.test_flag(VISITED);
	}

	bool set_visited(bool visited) {
		if (visited)
			return flags.set_flag(VISITED);
		else
			return flags.clear_flag(VISITED);
	}

	void run(vertex_program &prog) {
		if (!has_visited()) {
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
	assert(!has_visited());
	set_visited(true);

	int num_dests = vertex.get_num_edges(traverse_edge);
	if (num_dests == 0)
		return;

	// We need to add the neighbors of the vertex to the queue of
	// the next level.
#ifdef USE_ARRAY
	stack_array<vertex_id_t, 1024> neighs(num_dests);
	vertex.read_edges(traverse_edge, neighs.data(), num_dests);
	prog.activate_vertices(neighs.data(), num_dests);
#else
	edge_seq_iterator it = vertex.get_neigh_seq_it(traverse_edge, 0, num_dests);
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
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

void print_usage()
{
	fprintf(stderr,
			"bfs [options] conf_file graph_file index_file start_vertex\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-p: preload the graph\n");
	fprintf(stderr, "-b: traverse with both in-edges and out-edges\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_opts = 0;
	bool preload = false;
	while ((opt = getopt(argc, argv, "c:pb")) != -1) {
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
	if (preload)
		graph->preload_graph();
	printf("BFS starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start(&start_vertex, 1);
	graph->wait4complete();
	gettimeofday(&end, NULL);

	vertex_query::ptr cvq(new count_vertex_query());
	graph->query_on_all(cvq);
	size_t num_visited = ((count_vertex_query *) cvq.get())->get_num_visited();

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	printf("BFS from vertex %ld visits %ld vertices. It takes %f seconds\n",
			(unsigned long) start_vertex, num_visited, time_diff(start, end));
}
