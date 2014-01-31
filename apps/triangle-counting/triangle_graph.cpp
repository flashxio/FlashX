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
	count_msg(int num): vertex_message(sizeof(count_msg), false) {
		this->num = num;
	}

	int get_num() const {
		return num;
	}
};

class triangle_vertex: public compute_vertex
{
	// The number of required vertices to join with this vertex.
	int num_required;
	// The number of vertices that have joined with the vertex.
	int num_joined;
	// The number of vertices that have been asked to fetch.
	int num_fetched;
	// The number of triangles per vertex.
	atomic_integer num_pv_triangles;
	// It contains part of in-edges.
	// We only use the neighbors whose ID is smaller than this vertex.
	std::vector<vertex_id_t> *in_edges;
	// The vector contains the number of the neighbors' triangles shared
	// with this vertex. It only keeps the triangles of neighbors in the
	// in-edges.
	std::vector<int> *triangles;
	// It contains part of out-edges.
	// We only read neighbors whose ID is smaller than this vertex.
	std::vector<vertex_id_t> *out_edges;
public:
	triangle_vertex(): compute_vertex(-1, -1, 0) {
		num_required = 0;
		num_joined = 0;
		num_fetched = 0;
		in_edges = NULL;
		triangles = NULL;
		out_edges = NULL;
	}

	triangle_vertex(vertex_id_t id, off_t off, int size): compute_vertex(
			id, off, size) {
		num_required = 0;
		num_joined = 0;
		num_fetched = 0;
		in_edges = NULL;
		triangles = NULL;
		out_edges = NULL;
	}

	int count_triangles(const page_vertex *v) const;

	virtual bool has_required_vertices() const {
		return num_fetched < num_required;
	}

	virtual vertex_id_t get_next_required_vertex() {
		return out_edges->at(num_fetched++);
	}

	int get_num_triangles() const {
		return num_pv_triangles.get();
	}

	bool run(graph_engine &graph, const page_vertex *vertex);

	bool run_on_neighbors(graph_engine &graph,
			const page_vertex *vertices[], int num);

	void run_on_messages(graph_engine &graph,
			const vertex_message *msgs[], int num) {
		int sum = 0;
		for (int i = 0; i < num; i++)
			sum += ((count_msg *) msgs[i])->get_num();
		num_pv_triangles.inc(sum);
	}
};

int triangle_vertex::count_triangles(const page_vertex *v) const
{
	int num_local_triangles = 0;
	assert(v->get_id() != this->get_id());

	if (v->get_num_edges(edge_type::OUT_EDGE) == 0)
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

	assert(in_edges);
	// If the neighbor vertex has way more edges than this vertex.
	if (v->get_num_edges(edge_type::OUT_EDGE) / in_edges->size() > BIN_SEARCH_RATIO) {
		page_byte_array::const_iterator<vertex_id_t> other_it
			= v->get_neigh_begin(edge_type::OUT_EDGE);
		page_byte_array::const_iterator<vertex_id_t> other_end
			= v->get_neigh_end(edge_type::OUT_EDGE);
		for (int i = in_edges->size() - 1; i >= 0; i--) {
			vertex_id_t this_neighbor = in_edges->at(i);
			// We need to skip loops.
			if (this_neighbor != v->get_id()
					&& this_neighbor != this->get_id()) {
				page_byte_array::const_iterator<vertex_id_t> first
					= std::lower_bound(other_it, other_end, this_neighbor);
				if (first != other_end && this_neighbor == *first) {
					num_local_triangles++;
					(triangles->at(i))++;
				}
				other_end = first;
			}
		}
	}
	// If this vertex has way more edges than the neighbor vertex.
	else if (in_edges->size() / v->get_num_edges(edge_type::OUT_EDGE)
			> BIN_SEARCH_RATIO) {
		// TODO the same as above.
		page_byte_array::const_iterator<vertex_id_t> other_it
			= v->get_neigh_begin(edge_type::OUT_EDGE);
		page_byte_array::const_iterator<vertex_id_t> other_end
			= v->get_neigh_end(edge_type::OUT_EDGE);
		while (other_it != other_end) {
			vertex_id_t neigh_neighbor = *other_it;
			if (neigh_neighbor != v->get_id()
					&& neigh_neighbor != this->get_id()) {
				std::vector<vertex_id_t>::const_iterator first = std::lower_bound(
						in_edges->begin(), in_edges->end(), neigh_neighbor);
				if (first != in_edges->end() && neigh_neighbor == *first) {
					num_local_triangles++;
					int distance = first - in_edges->begin();
					(triangles->at(distance))++;
				}
			}
			++other_it;
		}
	}
	else {
		std::vector<vertex_id_t>::const_iterator this_it = in_edges->begin();
		std::vector<int>::iterator count_it = triangles->begin();
		std::vector<vertex_id_t>::const_iterator this_end = in_edges->end();
		page_byte_array::const_iterator<vertex_id_t> other_it
			= v->get_neigh_begin(edge_type::OUT_EDGE);
		page_byte_array::const_iterator<vertex_id_t> other_end
			= v->get_neigh_end(edge_type::OUT_EDGE);
		while (this_it != this_end && other_it != other_end) {
			vertex_id_t this_neighbor = *this_it;
			vertex_id_t neigh_neighbor = *other_it;
			if (this_neighbor == neigh_neighbor) {
				// skip loop
				if (neigh_neighbor != v->get_id()
						&& neigh_neighbor != this->get_id()) {
					num_local_triangles++;
					(*count_it)++;
				}
				++this_it;
				++other_it;
				++count_it;
			}
			else if (this_neighbor < neigh_neighbor) {
				++this_it;
				++count_it;
			}
			else
				++other_it;
		}
	}
	return num_local_triangles;
}

// We only use the neighbors whose ID is smaller than this vertex.
static int get_required_edges(const page_vertex *vertex, edge_type type,
		std::vector<vertex_id_t> &edges)
{
	page_byte_array::const_iterator<vertex_id_t> it
		= vertex->get_neigh_begin(type);
	page_byte_array::const_iterator<vertex_id_t> end
		= vertex->get_neigh_end(type);
	int num = 0;
	for (; it != end; ++it) {
		vertex_id_t id = *it;
		if (id < vertex->get_id()) {
			edges.push_back(id);
			num++;
		}
	}
	return num;
}

bool triangle_vertex::run(graph_engine &graph, const page_vertex *vertex)
{
	assert(in_edges == NULL);
	assert(out_edges == NULL);
	assert(num_joined == 0);
	assert(num_fetched == 0);

	long ret = num_working_vertices.inc(1);
	if (ret % 100000 == 0)
		printf("%ld working vertices\n", ret);
	// A vertex has to have in-edges and out-edges in order to form
	// a triangle. so we can simply skip the vertices that don't have
	// either of them.
	if (vertex->get_num_edges(edge_type::OUT_EDGE) == 0
			|| vertex->get_num_edges(edge_type::IN_EDGE) == 0) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		return true;
	}

	in_edges = new std::vector<vertex_id_t>();
	out_edges = new std::vector<vertex_id_t>();
	get_required_edges(vertex, edge_type::IN_EDGE, *in_edges);
	get_required_edges(vertex, edge_type::OUT_EDGE, *out_edges);
	num_required = out_edges->size();

	if (in_edges->empty() || out_edges->empty()) {
		num_required = 0;
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		delete in_edges;
		delete out_edges;
		in_edges = NULL;
		out_edges = NULL;
		return true;
	}
	else {
		triangles = new std::vector<int>(in_edges->size());
	}
	return false;
}

bool triangle_vertex::run_on_neighbors(graph_engine &graph,
		const page_vertex *vertices[], int num)
{
	num_joined++;
	for (int i = 0; i < num; i++) {
		int ret = count_triangles(vertices[i]);
		// If we find triangles with the neighbor, notify the neighbor
		// as well.
		if (ret > 0) {
			num_pv_triangles.inc(ret);
			count_msg msg(ret);
			graph.send_msg(vertices[i]->get_id(), msg);
		}
	}

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (num_joined == num_required) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);

		// Inform all neighbors in the in-edges.
		for (size_t i = 0; i < triangles->size(); i++) {
			// Inform the neighbor if they share triangles.
			if (triangles->at(i) > 0) {
				count_msg msg(triangles->at(i));
				graph.send_msg(in_edges->at(i), msg);
			}
		}

		delete triangles;
		delete in_edges;
		delete out_edges;
		in_edges = NULL;
		out_edges = NULL;
		return true;
	}
	return false;
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

	int min_vertex_size;
	if (directed)
		min_vertex_size = sizeof(ext_mem_directed_vertex);
	else
		min_vertex_size = sizeof(ext_mem_undirected_vertex);

	graph_index *index = graph_index_impl<triangle_vertex>::create(
			index_file, min_vertex_size);
	ext_mem_vertex_interpreter *interpreter;
	if (directed)
		interpreter = new ext_mem_directed_vertex_interpreter();
	else
		interpreter = new ext_mem_undirected_vertex_interpreter();
	graph_engine *graph = graph_engine::create(
			graph_conf.get_num_threads(), params.get_num_nodes(),
			graph_file, index, interpreter, directed);
	graph->set_required_neighbor_type(edge_type::OUT_EDGE);
	printf("triangle counting starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all();
	gettimeofday(&end, NULL);
	float activate_time = time_diff(start, end);

	start = end;
	graph->wait4complete();
	gettimeofday(&end, NULL);

	// Count the total number of triangles in the graph.
	graph_index::const_iterator it = index->begin();
	graph_index::const_iterator end_it = index->end();
	long num_triangles = 0;
	for (; it != end_it; ++it) {
		const triangle_vertex &v = (const triangle_vertex &) *it;
		num_triangles += v.get_num_triangles();
	}

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph_engine::destroy(graph);
	printf("There are %ld vertices\n", index->get_num_vertices());
	printf("process %ld vertices and complete %ld vertices\n",
			num_working_vertices.get(), num_completed_vertices.get());
	printf("there are %ld triangles.\n", num_triangles);
	printf("It takes %f seconds to activate all and %f seconds to finish\n",
			activate_time, time_diff(start, end));
}
