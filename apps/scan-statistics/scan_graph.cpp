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

const double BIN_SEARCH_RATIO = 100;

atomic_number<long> num_working_vertices;
atomic_number<long> num_completed_vertices;

class count_msg: public vertex_message
{
	int num;
public:
	count_msg(int num) {
		this->num = num;
	}

	int get_num() const {
		return num;
	}
};

class scan_vertex: public compute_vertex
{
	// The number of required vertices to join with this vertex.
	int num_required;
	// The number of vertices that have joined with the vertex.
	int num_joined;
	// The number of vertices that have been asked to fetch.
	int num_fetched;
	// The number of edges in its neighborhood.
	atomic_integer num_edges;
	// All neighbors (in both in-edges and out-edges)
	std::vector<vertex_id_t> *neighbors;
public:
	scan_vertex(): compute_vertex(-1, -1, 0) {
		num_required = 0;
		num_joined = 0;
		num_fetched = 0;
		neighbors = NULL;
	}

	scan_vertex(vertex_id_t id, off_t off, int size): compute_vertex(
			id, off, size) {
		num_required = 0;
		num_joined = 0;
		num_fetched = 0;
		neighbors = NULL;
	}

	int count_edges(const page_vertex *v) const;

	virtual bool has_required_vertices() const {
		return num_fetched < num_required;
	}

	virtual vertex_id_t get_next_required_vertex() {
		return neighbors->at(num_fetched++);
	}

	// We only use the neighbors whose ID is smaller than this vertex.
	int get_required_edges(edge_type type, std::vector<vertex_id_t> &edges) {
		page_byte_array::const_iterator<vertex_id_t> it = get_neigh_begin(type);
		page_byte_array::const_iterator<vertex_id_t> end = get_neigh_end(type);
		int num = 0;
		for (; it != end; ++it) {
			vertex_id_t id = *it;
			if (id != get_id()) {
				edges.push_back(id);
				num++;
			}
		}
		return num;
	}

	int get_num_edges_in_neigh() const {
		return num_edges.get();
	}

	void run(graph_engine &graph);

	void run_on_neighbors(graph_engine &graph,
			const page_vertex *vertices[], int num);

	void run_on_messages(graph_engine &graph,
			const vertex_message *msgs[], int num) {
	}
};

int scan_vertex::count_edges(const page_vertex *v) const
{
	int num_local_edges = 0;
	assert(v->get_id() != this->get_id());

	if (v->get_num_edges(edge_type::BOTH_EDGES) == 0)
		return 0;

	/*
	 * We search for triangles with two different ways:
	 * binary search if two adjacency lists have very different sizes,
	 * scan otherwise.
	 *
	 * when binary search for multiple neighbors, we can reduce binary search
	 * overhead by using the new end in the search range. We can further reduce
	 * overhead by searching in a reverse order (start from the largest neighbor).
	 * Since vertices of smaller ID has more neighbors, it's more likely
	 * that a neighbor is in the beginning of the adjacency list, and
	 * the search range will be narrowed faster.
	 */

	assert(neighbors);
	// If the neighbor vertex has way more edges than this vertex.
	if (v->get_num_edges(edge_type::BOTH_EDGES) / neighbors->size() > BIN_SEARCH_RATIO) {
		page_byte_array::const_iterator<vertex_id_t> other_it
			= v->get_neigh_begin(edge_type::BOTH_EDGES);
		page_byte_array::const_iterator<vertex_id_t> other_end
			= v->get_neigh_end(edge_type::BOTH_EDGES);
		for (int i = neighbors->size() - 1; i >= 0; i--) {
			vertex_id_t this_neighbor = neighbors->at(i);
			// We need to skip loops.
			if (this_neighbor != v->get_id()
					&& this_neighbor != this->get_id()) {
				page_byte_array::const_iterator<vertex_id_t> first
					= std::lower_bound(other_it, other_end, this_neighbor);
				if (first != other_end && this_neighbor == *first) {
					num_local_edges++;
				}
				other_end = first;
			}
		}
	}
	// If this vertex has way more edges than the neighbor vertex.
	else if (neighbors->size() / v->get_num_edges(edge_type::BOTH_EDGES)
			> BIN_SEARCH_RATIO) {
		// TODO the same as above.
		page_byte_array::const_iterator<vertex_id_t> other_it
			= v->get_neigh_begin(edge_type::BOTH_EDGES);
		page_byte_array::const_iterator<vertex_id_t> other_end
			= v->get_neigh_end(edge_type::BOTH_EDGES);
		while (other_it != other_end) {
			vertex_id_t neigh_neighbor = *other_it;
			if (neigh_neighbor != v->get_id()
					&& neigh_neighbor != this->get_id()) {
				std::vector<vertex_id_t>::const_iterator first = std::lower_bound(
						neighbors->begin(), neighbors->end(), neigh_neighbor);
				if (first != neighbors->end() && neigh_neighbor == *first) {
					num_local_edges++;
				}
			}
			++other_it;
		}
	}
	else {
		std::vector<vertex_id_t>::const_iterator this_it = neighbors->begin();
		std::vector<vertex_id_t>::const_iterator this_end = neighbors->end();
		page_byte_array::const_iterator<vertex_id_t> other_it
			= v->get_neigh_begin(edge_type::BOTH_EDGES);
		page_byte_array::const_iterator<vertex_id_t> other_end
			= v->get_neigh_end(edge_type::BOTH_EDGES);
		while (this_it != this_end && other_it != other_end) {
			vertex_id_t this_neighbor = *this_it;
			vertex_id_t neigh_neighbor = *other_it;
			if (this_neighbor == neigh_neighbor) {
				// skip loop
				if (neigh_neighbor != v->get_id()
						&& neigh_neighbor != this->get_id()) {
					num_local_edges++;
				}
				++this_it;
				++other_it;
			}
			else if (this_neighbor < neigh_neighbor) {
				++this_it;
			}
			else
				++other_it;
		}
	}
	return num_local_edges;
}

void scan_vertex::run(graph_engine &graph)
{
	assert(neighbors == NULL);
	assert(num_joined == 0);
	assert(num_fetched == 0);

	long ret = num_working_vertices.inc(1);
	if (ret % 100000 == 0)
		printf("%ld working vertices\n", ret);
	if (get_num_edges(edge_type::BOTH_EDGES) == 0) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		return;
	}

	neighbors = new std::vector<vertex_id_t>();
	this->get_required_edges(edge_type::BOTH_EDGES, *neighbors);
	num_required = neighbors->size();
	num_edges.inc(neighbors->size());
}

void scan_vertex::run_on_neighbors(graph_engine &graph,
		const page_vertex *vertices[], int num)
{
	num_joined++;
	for (int i = 0; i < num; i++) {
		int ret = count_edges(vertices[i]);
		// If we find triangles with the neighbor, notify the neighbor
		// as well.
		if (ret > 0) {
			num_edges.inc(ret);
		}
	}

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (num_joined == num_required) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);

//		printf("v%ld (%ld neighbors): %d\n", get_id(), neighbors->size(),
//				num_edges.get());
		delete neighbors;
		neighbors = NULL;
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
		fprintf(stderr, "scan-statistics conf_file graph_file index_file directed [output_file]\n");
		graph_conf.print_help();
		params.print_help();
		exit(-1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];
	bool directed = atoi(argv[4]);
	assert(directed);
	std::string output_file;
	if (argc == 6) {
		output_file = argv[5];
		argc--;
	}

	config_map configs(conf_file);
	configs.add_options(argv + 5, argc - 5);
	graph_conf.init(configs);
	graph_conf.print();

	signal(SIGINT, int_handler);
	init_io_system(configs);

	graph_index *index = graph_index_impl<scan_vertex>::create(
			index_file, directed);
	graph_engine *graph = graph_engine::create(
			graph_conf.get_num_threads(), params.get_num_nodes(),
			graph_file, index, directed);
	// TODO I need to redefine this interface.
	graph->set_required_neighbor_type(edge_type::BOTH_EDGES);
	printf("scan statistics starts\n");
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
	printf("It takes %f seconds\n", time_diff(start, end));
	printf("There are %ld vertices\n", index->get_num_vertices());
	printf("process %ld vertices and complete %ld vertices\n",
			num_working_vertices.get(), num_completed_vertices.get());

	if (!output_file.empty()) {
		FILE *f = fopen(output_file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			return -1;
		}
		std::vector<vertex_id_t> vertices;
		index->get_all_vertices(vertices);
		for (size_t i = 0; i < index->get_num_vertices(); i++) {
			scan_vertex &v = (scan_vertex &) index->get_vertex(vertices[i]);
			fprintf(f, "v%ld: %d\n", v.get_id(), v.get_num_edges_in_neigh());
		}
		fclose(f);
	}
}
