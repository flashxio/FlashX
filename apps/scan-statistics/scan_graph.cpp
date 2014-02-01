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

#include <set>
#include <vector>

#include "thread.h"
#include "io_interface.h"

#include "graph_engine.h"
#include "graph_config.h"

const double BIN_SEARCH_RATIO = 10;

atomic_number<long> num_working_vertices;
atomic_number<long> num_completed_vertices;

class scan_vertex: public compute_vertex
{
	// The number of vertices that have joined with the vertex.
	int num_joined;
	std::vector<vertex_id_t>::const_iterator fetch_it;
	atomic_integer num_edges;
	// All neighbors (in both in-edges and out-edges)
	std::vector<vertex_id_t> *neighbors;
public:
	scan_vertex(): compute_vertex(-1, -1, 0) {
		num_joined = 0;
		neighbors = NULL;
	}

	scan_vertex(vertex_id_t id, off_t off, int size): compute_vertex(
			id, off, size) {
		num_joined = 0;
		neighbors = NULL;
	}

	int get_result() const {
		return num_edges.get();
	}

	virtual bool has_required_vertices() const {
		if (neighbors == NULL)
			return false;
		return fetch_it != neighbors->end();
	}

	virtual vertex_id_t get_next_required_vertex() {
		vertex_id_t id = *fetch_it;
		fetch_it++;
		return id;
	}

	int count_edges(const page_vertex *v,
			const std::vector<vertex_id_t> *neighbors);
	int count_edges(const page_vertex *v,
			const std::vector<vertex_id_t> *neighbors, edge_type type);

	bool run(graph_engine &graph, const page_vertex *vertex);

	bool run_on_neighbors(graph_engine &graph,
			const page_vertex *vertices[], int num);

	void run_on_messages(graph_engine &graph,
			const vertex_message *msgs[], int num) {
	}
};

int scan_vertex::count_edges(const page_vertex *v,
		const std::vector<vertex_id_t> *neighbors, edge_type type)
{
	int num_local_edges = 0;
	int num_v_edges = v->get_num_edges(type);
	if (num_v_edges == 0)
		return 0;

	page_byte_array::const_iterator<vertex_id_t> other_it
		= v->get_neigh_begin(type);
#if 0
	page_byte_array::const_iterator<edge_count> other_data_it
		= v->get_edge_data_begin<edge_count>(type);
#endif
	page_byte_array::const_iterator<vertex_id_t> other_end
		= std::lower_bound(other_it, v->get_neigh_end(type),
				v->get_id());

	std::vector<vertex_id_t>::const_iterator this_it = neighbors->begin();
	std::vector<vertex_id_t>::const_iterator this_end
		= std::lower_bound(this_it, neighbors->end(), v->get_id());

	if (num_v_edges / neighbors->size() > BIN_SEARCH_RATIO) {
		for (; this_it != this_end; this_it++) {
			vertex_id_t this_neighbor = *this_it;
			// We need to skip loops.
			if (this_neighbor == v->get_id()
					|| this_neighbor == this->get_id()) {
				continue;
			}

			page_byte_array::const_iterator<vertex_id_t> first
				= std::lower_bound(other_it, other_end, this_neighbor);
			// found it.
			if (first != other_end && !(this_neighbor < *first)) {
				do {
#if 0
					page_byte_array::const_iterator<edge_count> data_it
						= other_data_it;
					data_it += first - other_it;
					// Edges in the v's neighbor lists may duplicated.
					// The duplicated neighbors need to be counted
					// multiple times.
					num_local_edges += (*data_it).get_count();
					++data_it;
#endif
					num_local_edges++;
					++first;
				} while (first != other_end && this_neighbor == *first);
			}
		}
	}
	else if (neighbors->size() / num_v_edges > BIN_SEARCH_RATIO) {
		while (other_it != other_end) {
			vertex_id_t neigh_neighbor = *other_it;
			if (neigh_neighbor != v->get_id()
					&& neigh_neighbor != this->get_id()) {
				if (std::binary_search(this_it, this_end, neigh_neighbor)) {
#if 0
					num_local_edges += (*other_data_it).get_count();
#endif
					num_local_edges++;
				}
			}
			++other_it;
#if 0
			++other_data_it;
#endif
		}
	}
	else {
		while (other_it != other_end && this_it != this_end) {
			vertex_id_t this_neighbor = *this_it;
			vertex_id_t neigh_neighbor = *other_it;
			if (neigh_neighbor == v->get_id()
					|| neigh_neighbor == this->get_id()) {
				++other_it;
#if 0
				++other_data_it;
#endif
				continue;
			}
			if (this_neighbor == neigh_neighbor) {
				do {
					// Edges in the v's neighbor lists may duplicated.
					// The duplicated neighbors need to be counted
					// multiple times.
#if 0
					num_local_edges += (*other_data_it).get_count();
					++other_data_it;
#endif
					num_local_edges++;
					++other_it;
				} while (other_it != other_end && this_neighbor == *other_it);
				++this_it;
			}
			else if (this_neighbor < neigh_neighbor) {
				++this_it;
			}
			else {
				++other_it;
#if 0
				++other_data_it;
#endif
			}
		}
	}
	return num_local_edges;
}

int scan_vertex::count_edges(const page_vertex *v,
		const std::vector<vertex_id_t> *neighbors)
{
	if (v->get_num_edges(edge_type::BOTH_EDGES) == 0
			|| neighbors->empty())
		return 0;

	return count_edges(v, neighbors, edge_type::IN_EDGE)
		+ count_edges(v, neighbors, edge_type::OUT_EDGE);
}

template<class InputIterator1, class InputIterator2, class Skipper,
	class OutputIterator>
int unique_merge(InputIterator1 it1, InputIterator1 last1,
		InputIterator2 it2, InputIterator2 last2, Skipper skip,
		OutputIterator result)
{
	OutputIterator result_begin = result;
	while (it1 != last1 && it2 != last2) {
		if (*it1 > *it2) {
			typename std::iterator_traits<InputIterator2>::value_type v = *it2;
			if (!skip(v))
				*(result++) = v;
			while (*it2 == v && it2 != last2)
				++it2;
		}
		else if (*it1 < *it2) {
			typename std::iterator_traits<InputIterator1>::value_type v = *it1;
			if (!skip(v))
				*(result++) = v;
			while (*it1 == v && it1 != last1)
				++it1;
		}
		else {
			typename std::iterator_traits<InputIterator1>::value_type v = *it1;
			if (!skip(v))
				*(result++) = v;
			while (*it1 == v && it1 != last1)
				++it1;
			while (*it2 == v && it2 != last2)
				++it2;
		}
	}

	while (it1 != last1) {
		typename std::iterator_traits<InputIterator1>::value_type v = *it1;
		if (!skip(v))
			*(result++) = v;
		while (*it1 == v && it1 != last1)
			++it1;
	}

	while (it2 != last2) {
		typename std::iterator_traits<InputIterator2>::value_type v = *it2;
		if (!skip(v))
			*(result++) = v;
		while (*it2 == v && it2 != last2)
			++it2;
	}
	return result - result_begin;
}

bool scan_vertex::run(graph_engine &graph, const page_vertex *vertex)
{
	assert(neighbors == NULL);
	assert(num_joined == 0);
	assert(num_edges.get() == 0);

	long ret = num_working_vertices.inc(1);
	if (ret % 100000 == 0)
		printf("%ld working vertices\n", ret);
	if (vertex->get_num_edges(edge_type::BOTH_EDGES) == 0) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		return true;
	}

	neighbors = new std::vector<vertex_id_t>(vertex->get_num_edges(
				edge_type::BOTH_EDGES));

	class skip_self {
		vertex_id_t id;
	public:
		skip_self(vertex_id_t id) {
			this->id = id;
		}

		bool operator()(vertex_id_t id) {
			return this->id == id;
		}
	};

	int num_neighbors = unique_merge(
			vertex->get_neigh_begin(edge_type::IN_EDGE),
			vertex->get_neigh_end(edge_type::IN_EDGE),
			vertex->get_neigh_begin(edge_type::OUT_EDGE),
			vertex->get_neigh_end(edge_type::OUT_EDGE),
			skip_self(vertex->get_id()),
			neighbors->begin());
	neighbors->resize(num_neighbors);

	page_byte_array::const_iterator<vertex_id_t> it = vertex->get_neigh_begin(
			edge_type::BOTH_EDGES);
#if 0
	page_byte_array::const_iterator<edge_count> data_it
		= vertex->get_edge_data_begin<edge_count>(
				edge_type::BOTH_EDGES);
#endif
	page_byte_array::const_iterator<vertex_id_t> end = vertex->get_neigh_end(
			edge_type::BOTH_EDGES);
#if 0
	for (; it != end; ++it, ++data_it) {
#endif
	for (; it != end; ++it) {
		vertex_id_t id = *it;
		// Ignore loops
		if (id != vertex->get_id()) {
#if 0
			num_edges.inc((*data_it).get_count());
#endif
			num_edges.inc(1);
		}
	}

	if (neighbors->size() == 0) {
		delete neighbors;
		neighbors = NULL;
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		return true;
	}

	fetch_it = neighbors->begin();
	return false;
}

bool scan_vertex::run_on_neighbors(graph_engine &graph,
		const page_vertex *vertices[], int num)
{
	num_joined += num;
	assert(neighbors);
	for (int i = 0; i < num; i++) {
		int ret = count_edges(vertices[i], neighbors);
		if (ret > 0)
			num_edges.inc(ret);
	}

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (num_joined == (int) neighbors->size()) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);

		delete neighbors;
		neighbors = NULL;
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
		fprintf(stderr,
				"scan-statistics conf_file graph_file index_file directed [output_file]\n");
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
			index_file, sizeof(ext_mem_directed_vertex));
	graph_engine *graph = graph_engine::create(
			graph_conf.get_num_threads(), params.get_num_nodes(), graph_file,
			index, new ext_mem_directed_vertex_interpreter(), directed);
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
	graph_engine::destroy(graph);
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
		graph_index::const_iterator it = index->begin();
		graph_index::const_iterator end_it = index->end();
		for (; it != end_it; ++it) {
			const scan_vertex &v = (const scan_vertex &) *it;
			fprintf(f, "\"%ld\" %d\n", (unsigned long) v.get_id(), v.get_result());
		}
		fclose(f);
	}
}
