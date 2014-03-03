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

/**
 * This program computes the overlap of each pair of vertices
 * given by the user.
 */

#include <signal.h>
#include <google/profiler.h>

#include <tr1/unordered_set>

#include "graph_engine.h"

enum overlap_stage_t
{
	CONSTRUCT_NEIGHBORS,
	COMP_OVERLAP,
} overlap_stage;

/**
 * This contains all vertices that we want to compute pair-wise overlap.
 */
std::vector<vertex_id_t> overlap_vertices;

template<class InputIterator1, class InputIterator2, class Skipper,
	class Merger, class OutputIterator>
size_t unique_merge(InputIterator1 it1, InputIterator1 last1,
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

class skip_self
{
	vertex_id_t id;
public:
	skip_self(vertex_id_t id) {
		this->id = id;
	}

	bool operator()(vertex_id_t id) {
		return this->id == id;
	}
};

class merge_edge
{
public:
	vertex_id_t operator()(vertex_id_t e1, vertex_id_t e2) {
		assert(e1 == e2);
		return e1;
	}
};

size_t get_unique_neighbors(const page_vertex &vertex,
		std::vector<vertex_id_t> &neighbors)
{
	neighbors.resize(vertex.get_num_edges(edge_type::BOTH_EDGES));
	size_t num_neighbors = unique_merge(
			vertex.get_neigh_begin(edge_type::IN_EDGE),
			vertex.get_neigh_end(edge_type::IN_EDGE),
			vertex.get_neigh_begin(edge_type::OUT_EDGE),
			vertex.get_neigh_end(edge_type::OUT_EDGE),
			skip_self(vertex.get_id()), merge_edge(),
			neighbors.begin());
	neighbors.resize(num_neighbors);
	return num_neighbors;
}

size_t get_common_vertices(const std::vector<vertex_id_t> &vertices1,
		const std::vector<vertex_id_t> &vertices2)
{
	size_t common = 0;
	size_t i, j;
	// We assume both vectors are sorted in the ascending order.
	// We only need to scan them together to extract the common elements.
	for (i = 0, j = 0; i < vertices1.size() && j < vertices2.size();) {
		if (vertices1[i] == vertices2[j]) {
			common++;
			i++;
			j++;
		}
		else if (vertices1[i] > vertices2[j])
			j++;
		else
			i++;

	}
	return common;
}

size_t get_union_vertices(const std::vector<vertex_id_t> &vertices1,
		const std::vector<vertex_id_t> &vertices2)
{
	class skip_none {
	public:
		bool operator()(vertex_id_t id) {
			return false;
		}
	};

	class null_iterator: public std::iterator<std::random_access_iterator_tag, vertex_id_t> {
		vertex_id_t id;
		int idx;
	public:
		typedef typename std::iterator<std::random_access_iterator_tag,
				vertex_id_t>::difference_type difference_type;

		null_iterator() {
			idx = 0;
			id = 0;
		}

		vertex_id_t &operator*() {
			return id;
		}

		null_iterator operator++(int) {
			null_iterator it = *this;
			idx++;
			return it;
		}

		difference_type operator-(const null_iterator &it) {
			return idx - it.idx;
		}
	};

	size_t num_eles = unique_merge(vertices1.begin(), vertices1.end(),
			vertices2.begin(), vertices2.end(), skip_none(), merge_edge(),
			null_iterator());
	return num_eles;
}

class union_set
{
	std::tr1::unordered_set<vertex_id_t> set;
	pthread_mutex_t lock;
public:
	union_set() {
		pthread_mutex_init(&lock, NULL);
	}

	void add(const std::vector<vertex_id_t> &vec) {
		pthread_mutex_lock(&lock);
		set.insert(vec.begin(), vec.end());
		pthread_mutex_unlock(&lock);
	}

	size_t get_size() const {
		return set.size();
	}
} vertex_union;

class intersection_set
{
	int num_added;
	std::vector<vertex_id_t> set;
	pthread_mutex_t lock;
public:
	intersection_set() {
		num_added = 0;
		pthread_mutex_init(&lock, NULL);
	}

	void add(const std::vector<vertex_id_t> &vec);

	size_t get_size() const {
		return set.size();
	}

	vertex_id_t get(int idx) const {
		return set[idx];
	}
} vertex_intersection;

void intersection_set::add(const std::vector<vertex_id_t> &vec)
{
	std::vector<vertex_id_t> copy(vec);
	std::sort(copy.begin(), copy.end());

	std::vector<vertex_id_t> new_set;
	pthread_mutex_lock(&lock);
	num_added++;
	if (num_added == 1) {
		set = copy;
	}
	else {
		size_t i, j;
		// We assume both vectors are sorted in the ascending order.
		// We only need to scan them together to extract the common elements.
		for (i = 0, j = 0; i < set.size() && j < copy.size();) {
			if (set[i] == copy[j]) {
				new_set.push_back(set[i]);
				i++;
				j++;
			}
			else if (set[i] > copy[j])
				j++;
			else
				i++;

		}
		set = new_set;
	}
	pthread_mutex_unlock(&lock);
}

class overlap_vertex: public compute_vertex
{
	std::vector<vertex_id_t> *neighborhood;
public:
	overlap_vertex() {
		neighborhood = NULL;
	}

	overlap_vertex(vertex_id_t id, const vertex_index *index1): compute_vertex(
			id, index1) {
		neighborhood = NULL;
	}

	void run(graph_engine &graph) {
		switch(overlap_stage) {
			case overlap_stage_t::CONSTRUCT_NEIGHBORS:
				run_stage1(graph);
				break;
			case overlap_stage_t::COMP_OVERLAP:
				run_stage2(graph);
				break;
			default:
				assert(0);
		}
	}

	void run_stage1(graph_engine &graph) {
		vertex_id_t id = get_id();
		graph.request_vertices(*this, &id, 1);
	}

	void run_stage2(graph_engine &graph) {
		BOOST_FOREACH(vertex_id_t id, overlap_vertices) {
			if (id == get_id())
				continue;
			overlap_vertex &neigh = (overlap_vertex &) graph.get_vertex(id);
			size_t common = get_common_vertices(*neighborhood, *neigh.neighborhood);
			size_t vunion = get_union_vertices(*neighborhood, *neigh.neighborhood);
			printf("v%u:v%u, common: %ld, union: %ld, overlap: %f\n",
					get_id(), id, common, vunion,
					((double) common) / vunion);
		}
	}

	void run(graph_engine &graph, const page_vertex &vertex) {
		assert(vertex.get_id() == get_id());
		run_on_itself(graph, vertex);
	}

	void run_on_itself(graph_engine &graph, const page_vertex &vertex) {
		neighborhood = new std::vector<vertex_id_t>();
		get_unique_neighbors(vertex, *neighborhood);
		assert(std::is_sorted(neighborhood->begin(), neighborhood->end()));

		std::vector<vertex_id_t>::iterator it = std::lower_bound(
				neighborhood->begin(), neighborhood->end(), get_id());
		if (it != neighborhood->end())
			assert(*it != get_id());
		neighborhood->insert(it, get_id());
		assert(std::is_sorted(neighborhood->begin(), neighborhood->end()));
	}

	void run_on_messages(graph_engine &,
			const vertex_message *msgs[], int num) {
	}
};

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

int read_vertices(const std::string &file, std::vector<vertex_id_t> &vertices)
{
	FILE *f = fopen(file.c_str(), "r");
	assert(f);
	ssize_t ret;
	char *line = NULL;
	size_t line_size = 0;
	while ((ret = getline(&line, &line_size, f)) > 0) {
		if (line[ret - 1] == '\n')
			line[ret - 1] = 0;
		vertex_id_t id = atol(line);
		printf("%u\n", id);
		vertices.push_back(id);
	}
	fclose(f);
	return vertices.size();
}

int main(int argc, char *argv[])
{
	std::string output_file;
	std::string confs;

	if (argc < 3) {
		fprintf(stderr,
				"overlap conf_file graph_file index_file vertex_file\n");
		exit(-1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];
	std::string vertex_file = argv[4];

	read_vertices(vertex_file, overlap_vertices);

	config_map configs(conf_file);
	configs.add_options(confs);
	graph_conf.init(configs);
	graph_conf.print();

	signal(SIGINT, int_handler);
	init_io_system(configs);

	graph_index *index = NUMA_graph_index<overlap_vertex>::create(
			index_file, graph_conf.get_num_threads(), params.get_num_nodes());
	graph_engine *graph = graph_engine::create(
			graph_conf.get_num_threads(), params.get_num_nodes(), graph_file,
			index);
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	struct timeval start, end;

	gettimeofday(&start, NULL);
	overlap_stage = overlap_stage_t::CONSTRUCT_NEIGHBORS;
	graph->start(overlap_vertices.data(), overlap_vertices.size());
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to construct neighborhoods\n",
			time_diff(start, end));

	gettimeofday(&start, NULL);
	overlap_stage = overlap_stage_t::COMP_OVERLAP;
	graph->start(overlap_vertices.data(), overlap_vertices.size());
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to compute overlaps\n", time_diff(start, end));

	printf("All vertices have %ld common neighbors and cover %ld vertices\n",
			vertex_intersection.get_size(), vertex_union.get_size());
	for (size_t i = 0; i < vertex_intersection.get_size(); i++)
		printf("%u\n", vertex_intersection.get(i));

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph_engine::destroy(graph);
}
