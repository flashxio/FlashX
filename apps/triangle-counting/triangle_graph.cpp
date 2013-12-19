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

#include "graph_engine.h"
#include "graph_config.h"

atomic_number<long> num_triangles;

class triangle_vertex: public compute_vertex
{
	// The number of required vertices to join with this vertex.
	int num_required;
	// The number of vertices that have joined with the vertex.
	int num_joined;
	// The number of vertices that have been asked to fetch.
	int num_fetched;
	std::vector<vertex_id_t> in_edges;
	std::vector<vertex_id_t> out_edges;
public:
	triangle_vertex(): compute_vertex(-1, -1, 0) {
		num_required = 0;
		num_joined = 0;
		num_fetched = 0;
	}

	triangle_vertex(vertex_id_t id, off_t off, int size): compute_vertex(
			id, off, size) {
		num_required = 0;
		num_joined = 0;
		num_fetched = 0;
	}

	int count_triangle(const ext_mem_vertex &neighbor) const;

	virtual bool has_required_vertices() const {
		return num_fetched < num_required;
	}

	virtual vertex_id_t get_next_required_vertex() {
		return out_edges[num_fetched++];
	}

	void run(graph_engine &graph, const page_vertex *vertices[], int num);
};

void triangle_vertex::run(graph_engine &graph, const page_vertex *vertices[],
			int num)
{
	// If there aren't neighbors passed to the vertex's user code,
	// it's because the vertex doesn't have neighbors.
	if (num == 0) {
		assert(in_edges.size() == 0);
		assert(out_edges.size() == 0);
		// A vertex has to have in-edges and out-edges in order to form
		// a triangle. so we can simply skip the vertices that don't have
		// either of them.
		if (get_num_edges(edge_type::OUT_EDGE) == 0
				|| get_num_edges(edge_type::IN_EDGE) == 0)
			return;

		this->get_all_edges(edge_type::IN_EDGE, in_edges);
		this->get_all_edges(edge_type::OUT_EDGE, out_edges);
		num_required = out_edges.size();
		num_joined = 0;
		num_fetched = 0;
		return;
	}

	int num_local_triangles = 0;
	num_joined++;
	for (int i = 0; i < num; i++) {
		const page_vertex *v = vertices[i];
		std::vector<vertex_id_t>::const_iterator this_it = in_edges.begin();
		std::vector<vertex_id_t>::const_iterator this_end = in_edges.end();
		page_byte_array::const_iterator<vertex_id_t> other_it
			= v->get_neigh_begin(edge_type::OUT_EDGE);
		page_byte_array::const_iterator<vertex_id_t> other_end
			= v->get_neigh_end(edge_type::OUT_EDGE);
		while (this_it != this_end && other_it != other_end) {
			vertex_id_t this_neighbor = *this_it;
			vertex_id_t neigh_neighbor = *other_it;
			if (this_neighbor == neigh_neighbor) {
				num_local_triangles++;
				++this_it;
				++other_it;
			}
			else if (this_neighbor < neigh_neighbor)
				++this_it;
			else
				++other_it;
		}
	}
	if (num_local_triangles > 0)
		num_triangles.inc(num_local_triangles);

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (num_joined == num_required) {
		in_edges = std::vector<vertex_id_t>();
		out_edges = std::vector<vertex_id_t>();
	}
}

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

int main(int argc, char *argv[])
{
	if (argc < 5) {
		fprintf(stderr, "triangle-counting conf_file graph_file index_file directed\n");
		graph_conf.print_help();
		params.print_help();
		exit(-1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];
	bool directed = atoi(argv[4]);

	config_map configs(conf_file);
	configs.add_options(argv + 5, argc - 5);
	graph_conf.init(configs);
	graph_conf.print();

	signal(SIGINT, int_handler);
	init_io_system(configs);

	graph_index *index = graph_index_impl<triangle_vertex>::create(index_file);
	graph_engine *graph = graph_engine::create(graph_conf.get_num_threads(),
			params.get_num_nodes(), graph_file, index, directed);
	graph->set_required_neighbor_type(edge_type::OUT_EDGE);
	printf("triangle counting starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph->cleanup();
	printf("there are %ld triangles. It takes %f seconds\n",
			num_triangles.get(), time_diff(start, end));
}
