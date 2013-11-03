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

#include "thread.h"
#include "io_interface.h"

#include "bfs_graph.h"
#include "graph_config.h"

atomic_integer num_visited_vertices;

void bfs_vertex::run(graph_engine &graph, const ext_mem_vertex vertices[],
			int num)
{
	vertex_id_t max_id = graph.get_max_vertex_id();
	vertex_id_t min_id = graph.get_min_vertex_id();

	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	std::vector<vertex_id_t> activated_vertices;
	for (int j = 0; j < this->get_num_edges(OUT_EDGE); j++) {
		vertex_id_t id = this->get_edge(OUT_EDGE, j).get_to();
		assert(id >= min_id && id <= max_id);
		bfs_vertex &info = (bfs_vertex &) graph.get_vertex(id);
		// If the vertex has been visited, we can skip it.
		if (info.has_visited())
			continue;
		if (info.set_visited(true))
			continue;
		activated_vertices.push_back(id);
	}
	if (activated_vertices.size() > 0)
		num_visited_vertices.inc(activated_vertices.size());

	graph.activate_vertices(activated_vertices.data(),
			activated_vertices.size());
}

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

int main(int argc, char *argv[])
{
	if (argc < 6) {
		fprintf(stderr, "bfs conf_file graph_file index_file start_vertex directed\n");
		graph_conf.print_help();
		params.print_help();
		exit(-1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];
	vertex_id_t start_vertex = atoi(argv[4]);
	bool directed = atoi(argv[5]);

	config_map configs(conf_file);
	configs.add_options(argv + 6, argc - 6);
	graph_conf.init(configs);
	graph_conf.print();

	signal(SIGINT, int_handler);
	init_io_system(configs);

	graph_index *index = graph_index_impl<bfs_vertex>::create(index_file);
	bfs_graph *graph = bfs_graph::create(graph_conf.get_num_threads(),
			params.get_num_nodes(), graph_file, index, directed);
	printf("BFS starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start(&start_vertex, 1);
	graph->wait4complete();
	gettimeofday(&end, NULL);

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph->cleanup();
	printf("BFS from vertex %ld visits %d vertices. It takes %f seconds\n",
			start_vertex, num_visited_vertices.get(),
			time_diff(start, end));
}
