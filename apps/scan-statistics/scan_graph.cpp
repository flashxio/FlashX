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

class weighted_edge
{
	vertex_id_t id;
	int num_dups;
public:
	weighted_edge() {
		id = 0;
		num_dups = 0;
	}

	weighted_edge(vertex_id_t id) {
		this->id = id;
		num_dups = 1;
	}

	weighted_edge(vertex_id_t id, int num_dups) {
		this->id = id;
		this->num_dups = num_dups;
	}

	vertex_id_t get_id() const {
		return id;
	}

	int get_num_dups() const {
		return num_dups;
	}

	bool operator<(const weighted_edge &e) const {
		return this->id < e.id;
	}

	bool operator==(vertex_id_t id) const {
		return this->id == id;
	}

	bool operator==(const weighted_edge &e) const {
		return this->id == e.id;
	}
};

class comp_edge
{
public:
	bool operator()(const weighted_edge &e1, const weighted_edge &e2) {
		return e1.get_id() < e2.get_id();
	}
};

class scan_vertex: public compute_vertex
{
	// The number of vertices that have joined with the vertex.
	int num_joined;
	// The number of vertices required to join with the vertex.
	int num_required;
	std::vector<weighted_edge>::const_iterator fetch_it;
	std::vector<weighted_edge>::const_iterator fetch_end;
	atomic_integer num_edges;
	// All neighbors (in both in-edges and out-edges)
	std::vector<weighted_edge> *neighbors;
	std::vector<int> *edge_counts;
public:
	scan_vertex(): compute_vertex(-1, -1, 0) {
		num_joined = 0;
		num_required = 0;
		neighbors = NULL;
		edge_counts = NULL;
	}

	scan_vertex(vertex_id_t id, off_t off, int size): compute_vertex(
			id, off, size) {
		num_joined = 0;
		num_required = 0;
		neighbors = NULL;
		edge_counts = NULL;
	}

	int get_result() const {
		return num_edges.get();
	}

	virtual bool has_required_vertices() const {
		if (neighbors == NULL)
			return false;
		return fetch_it != fetch_end;
	}

	virtual vertex_id_t get_next_required_vertex() {
		vertex_id_t id = fetch_it->get_id();
		fetch_it++;
		return id;
	}

	int count_edges(graph_engine &graph, const page_vertex *v);
	int count_edges(graph_engine &graph, const page_vertex *v,
			edge_type type, std::vector<vertex_id_t> &common_neighs);

	bool run(graph_engine &graph, const page_vertex *vertex);

	bool run_on_neighbors(graph_engine &graph,
			const page_vertex *vertices[], int num);

	void run_on_messages(graph_engine &graph,
			const vertex_message *msgs[], int num) {
		for (int i = 0; i < num; i++) {
			const count_msg *msg = (const count_msg *) msgs[i];
			num_edges.inc(msg->get_num());
		}
	}
};

int scan_vertex::count_edges(graph_engine &graph, const page_vertex *v,
		edge_type type, std::vector<vertex_id_t> &common_neighs)
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

	std::vector<weighted_edge>::const_iterator this_it = neighbors->begin();
	std::vector<weighted_edge>::const_iterator this_end = neighbors->end();
	this_end = std::lower_bound(this_it, this_end,
			weighted_edge(v->get_id()), comp_edge());

	if (num_v_edges / neighbors->size() > BIN_SEARCH_RATIO) {
		for (; this_it != this_end; this_it++) {
			vertex_id_t this_neighbor = this_it->get_id();
			// We need to skip loops.
			if (this_neighbor == v->get_id()
					|| this_neighbor == this->get_id()) {
				continue;
			}

			page_byte_array::const_iterator<vertex_id_t> first
				= std::lower_bound(other_it, other_end, this_neighbor);
			// found it.
			if (first != other_end && !(this_neighbor < *first)) {
				int num_dups = 0;
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
					num_dups++;
					num_local_edges++;
					++first;
				} while (first != other_end && this_neighbor == *first);
				common_neighs.push_back(this_neighbor);
			}
		}
	}
	else if (neighbors->size() / num_v_edges > BIN_SEARCH_RATIO) {
		while (other_it != other_end) {
			vertex_id_t neigh_neighbor = *other_it;
			if (neigh_neighbor != v->get_id()
					&& neigh_neighbor != this->get_id()) {
				std::vector<weighted_edge>::const_iterator first
					= std::lower_bound(this_it, this_end,
							weighted_edge(neigh_neighbor), comp_edge());
				if (first != this_end && neigh_neighbor == first->get_id()) {
#if 0
					num_local_edges += (*other_data_it).get_count();
#endif
					num_local_edges++;
					common_neighs.push_back(first->get_id());
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
			vertex_id_t this_neighbor = this_it->get_id();
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
				common_neighs.push_back(this_it->get_id());
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

template<class InputIterator1, class InputIterator2, class Skipper,
	class Merger, class OutputIterator>
int unique_merge(InputIterator1 it1, InputIterator1 last1,
		InputIterator2 it2, InputIterator2 last2, Skipper skip,
		Merger merge, OutputIterator result)
{
	OutputIterator result_begin = result;
	while (it1 != last1 && it2 != last2) {
		if (*it2 < *it1) {
			typename std::iterator_traits<OutputIterator>::value_type v = *it2;
			++it2;
			while (it2 != last2 && v == *it2) {
				v = merge(v, *it2);
				++it2;
			}
			if (!skip(v))
				*(result++) = v;
		}
		else if (*it1 < *it2) {
			typename std::iterator_traits<OutputIterator>::value_type v = *it1;
			++it1;
			while (it1 != last1 && v == *it1) {
				v = merge(v, *it1);
				++it1;
			}
			if (!skip(v))
				*(result++) = v;
		}
		else {
			typename std::iterator_traits<OutputIterator>::value_type v = *it1;
			v = merge(v, *it2);
			++it2;
			while (it2 != last2 && v == *it2) {
				v = merge(v, *it2);
				++it2;
			}
			++it1;
			while (it1 != last1 && v == *it1) {
				v = merge(v, *it1);
				++it1;
			}
			if (!skip(v))
				*(result++) = v;
		}
	}

	while (it1 != last1) {
		typename std::iterator_traits<OutputIterator>::value_type v = *it1;
		++it1;
		while (it1 != last1 && v == *it1) {
			v = merge(v, *it1);
			++it1;
		}
		if (!skip(v))
			*(result++) = v;
	}

	while (it2 != last2) {
		typename std::iterator_traits<OutputIterator>::value_type v = *it2;
		++it2;
		while (it2 != last2 && v == *it2) {
			v = merge(v, *it2);
			++it2;
		}
		if (!skip(v))
			*(result++) = v;
	}
	return result - result_begin;
}

int scan_vertex::count_edges(graph_engine &graph, const page_vertex *v)
{
	assert(!neighbors->empty());
	if (v->get_num_edges(edge_type::BOTH_EDGES) == 0)
		return 0;

	std::vector<vertex_id_t> common_neighs1;
	std::vector<vertex_id_t> common_neighs2;
	int ret = count_edges(graph, v, edge_type::IN_EDGE, common_neighs1)
		+ count_edges(graph, v, edge_type::OUT_EDGE, common_neighs2);

	class skip_self {
	public:
		bool operator()(vertex_id_t id) {
			return false;
		}
	};

	class merge_edge {
	public:
		vertex_id_t operator()(vertex_id_t id1, vertex_id_t id2) {
			assert(id1 == id2);
			return id1;
		}
	};

	std::vector<vertex_id_t> common_neighs(common_neighs1.size()
			+ common_neighs2.size());
	int num_neighbors = unique_merge(
			common_neighs1.begin(), common_neighs1.end(),
			common_neighs2.begin(), common_neighs2.end(),
			skip_self(), merge_edge(), common_neighs.begin());
	common_neighs.resize(num_neighbors);

	// The number of duplicated edges between v and this vertex.
	std::vector<weighted_edge>::const_iterator v_it = std::lower_bound(
			neighbors->begin(), neighbors->end(),
			weighted_edge(v->get_id()), comp_edge());
	assert(v_it != neighbors->end());
	int v_idx = v_it - neighbors->begin();
	int num_v_dups = v_it->get_num_dups();
	assert(num_v_dups > 0);
	int num_edges = 0;
	for (std::vector<vertex_id_t>::const_iterator it = common_neighs.begin();
			it != common_neighs.end(); it++) {
		std::vector<weighted_edge>::const_iterator v_it = std::lower_bound(
				neighbors->begin(), neighbors->end(),
				weighted_edge(*it), comp_edge());
		assert(v_it != neighbors->end());
		int v_idx = v_it - neighbors->begin();
		num_edges += v_it->get_num_dups();

		edge_counts->at(v_idx) += num_v_dups;
	}
	if (num_edges > 0)
		edge_counts->at(v_idx) += num_edges;
	return ret;
}

bool scan_vertex::run(graph_engine &graph, const page_vertex *vertex)
{
	assert(neighbors == NULL);
	assert(num_joined == 0);

	long ret = num_working_vertices.inc(1);
	if (ret % 100000 == 0)
		printf("%ld working vertices\n", ret);
	if (vertex->get_num_edges(edge_type::BOTH_EDGES) == 0) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		return true;
	}

	neighbors = new std::vector<weighted_edge>(vertex->get_num_edges(
				edge_type::BOTH_EDGES));

	class skip_self {
		vertex_id_t id;
	public:
		skip_self(vertex_id_t id) {
			this->id = id;
		}

		bool operator()(weighted_edge &e) {
			return this->id == e.get_id();
		}

		bool operator()(vertex_id_t id) {
			return this->id == id;
		}
	};

	class merge_edge {
	public:
		weighted_edge operator()(const weighted_edge &e1,
				const weighted_edge &e2) {
			assert(e1.get_id() == e2.get_id());
			return weighted_edge(e1.get_id(),
					e1.get_num_dups() + e2.get_num_dups());
		}
	};

	int num_neighbors = unique_merge(
			vertex->get_neigh_begin(edge_type::IN_EDGE),
			vertex->get_neigh_end(edge_type::IN_EDGE),
			vertex->get_neigh_begin(edge_type::OUT_EDGE),
			vertex->get_neigh_end(edge_type::OUT_EDGE),
			skip_self(vertex->get_id()), merge_edge(),
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
	int tmp = 0;
	for (; it != end; ++it) {
		vertex_id_t id = *it;
		// Ignore loops
		if (id != vertex->get_id()) {
#if 0
			num_edges.inc((*data_it).get_count());
#endif
			tmp++;
		}
	}
	num_edges.inc(tmp);

	if (neighbors->size() == 0) {
		delete neighbors;
		neighbors = NULL;
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		return true;
	}

	fetch_it = neighbors->begin();
	// We use the same optimization as triangle counting.
	// We only count edges in the neighborhood of a vertex with the largest ID.
	fetch_end = std::lower_bound(neighbors->begin(), neighbors->end(),
			weighted_edge(get_id()), comp_edge());
	num_required = fetch_end - fetch_it;
	if (num_required == 0) {
		delete neighbors;
		neighbors = NULL;
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		return true;
	}
	neighbors->resize(num_required);
	fetch_end = neighbors->end();
	edge_counts = new std::vector<int>(num_required);
	return false;
}

bool scan_vertex::run_on_neighbors(graph_engine &graph,
		const page_vertex *vertices[], int num)
{
	num_joined += num;
	assert(neighbors);
	for (int i = 0; i < num; i++) {
		int ret = count_edges(graph, vertices[i]);
		if (ret > 0)
			num_edges.inc(ret);
	}

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (num_joined == num_required) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);

		// Inform all neighbors in the in-edges.
		for (size_t i = 0; i < edge_counts->size(); i++) {
			if (edge_counts->at(i) > 0) {
				count_msg msg(edge_counts->at(i));
				graph.send_msg(neighbors->at(i).get_id(), msg);
			}
		}

		delete neighbors;
		delete edge_counts;
		neighbors = NULL;
		edge_counts = NULL;
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
	if (argc < 4) {
		fprintf(stderr,
				"scan-statistics conf_file graph_file index_file [output_file]\n");
		graph_conf.print_help();
		params.print_help();
		exit(-1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];
	std::string output_file;
	if (argc == 5) {
		output_file = argv[4];
		argc--;
	}

	config_map configs(conf_file);
	configs.add_options(argv + 5, argc - 5);
	graph_conf.init(configs);
	graph_conf.print();

	signal(SIGINT, int_handler);
	init_io_system(configs);

	graph_index *index = graph_index_impl<scan_vertex>::create(
			index_file);
	graph_engine *graph = graph_engine::create(
			graph_conf.get_num_threads(), params.get_num_nodes(), graph_file,
			index);
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
