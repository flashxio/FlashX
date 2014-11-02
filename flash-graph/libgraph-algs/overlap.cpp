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

/**
 * This program computes the overlap of each pair of vertices
 * given by the user.
 */

#include <signal.h>

#include <tr1/unordered_set>

#include "FGlib.h"
#include "graph_engine.h"

namespace {
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

class overlap_vertex: public compute_vertex
{
	std::vector<vertex_id_t> *neighborhood;
public:
	overlap_vertex(vertex_id_t id): compute_vertex(id) {
		neighborhood = NULL;
	}

	void run(vertex_program &prog) {
		switch(overlap_stage) {
			case overlap_stage_t::CONSTRUCT_NEIGHBORS:
				run_stage1(prog);
				break;
			case overlap_stage_t::COMP_OVERLAP:
				run_stage2(prog);
				break;
			default:
				ABORT_MSG("wrong overlap stage");
		}
	}

	void run_stage1(vertex_program &prog) {
		vertex_id_t id = prog.get_vertex_id(*this);
		request_vertices(&id, 1);
	}
	void run_stage2(vertex_program &prog);

	void run(vertex_program &prog, const page_vertex &vertex) {
		assert(vertex.get_id() == prog.get_vertex_id(*this));
		run_on_itself(prog, vertex);
	}

	void run_on_itself(vertex_program &prog, const page_vertex &vertex) {
		neighborhood = new std::vector<vertex_id_t>();
		get_unique_neighbors(vertex, *neighborhood);
		assert(std::is_sorted(neighborhood->begin(), neighborhood->end()));

		vertex_id_t this_id = prog.get_vertex_id(*this);
		std::vector<vertex_id_t>::iterator it = std::lower_bound(
				neighborhood->begin(), neighborhood->end(), this_id);
		if (it != neighborhood->end())
			assert(*it != this_id);
		neighborhood->insert(it, this_id);
		assert(std::is_sorted(neighborhood->begin(), neighborhood->end()));
	}

	void run_on_message(vertex_program &, const vertex_message &msg) {
	}
};

class overlap_vertex_program: public vertex_program_impl<overlap_vertex>
{
	std::vector<std::vector<double> > &overlap_matrix;
public:
	overlap_vertex_program(
			std::vector<std::vector<double> > &_overlaps): overlap_matrix(_overlaps) {
	}

	void set_overlap(vertex_id_t id, const std::vector<double> &overlaps) {
		std::vector<vertex_id_t>::const_iterator it = std::lower_bound(
				overlap_vertices.begin(), overlap_vertices.end(), id);
		assert(it != overlap_vertices.end());
		assert(*it == id);
		off_t idx = it - overlap_vertices.begin();
		overlap_matrix[idx] = overlaps;
	}
};

class overlap_vertex_program_creater: public vertex_program_creater
{
	std::vector<std::vector<double> > &overlaps;
public:
	overlap_vertex_program_creater(
			std::vector<std::vector<double> > &_overlaps): overlaps(_overlaps) {
	}

	vertex_program::ptr create() const {
		return vertex_program::ptr(new overlap_vertex_program(overlaps));
	}
};

void overlap_vertex::run_stage2(vertex_program &prog)
{
	std::vector<double> overlaps(overlap_vertices.size());

	vertex_id_t this_id = prog.get_vertex_id(*this);
	for (size_t i = 0; i < overlap_vertices.size(); i++) {
		vertex_id_t id = overlap_vertices[i];
		if (id == this_id)
			continue;

		overlap_vertex &neigh = (overlap_vertex &) prog.get_graph().get_vertex(id);
		size_t common = get_common_vertices(*neighborhood, *neigh.neighborhood);
		size_t vunion = get_union_vertices(*neighborhood, *neigh.neighborhood);
		overlaps[i] = ((double) common) / vunion;
	}
	overlap_vertex_program &overlap_prog = (overlap_vertex_program &) prog;
	overlap_prog.set_overlap(this_id, overlaps);
}

}

void compute_overlap(FG_graph::ptr fg, const std::vector<vertex_id_t> &vids,
		std::vector<std::vector<double> > &overlap_matrix)
{
	assert(std::is_sorted(vids.begin(), vids.end()));
	graph_index::ptr index = NUMA_graph_index<overlap_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);
	overlap_vertices = vids;
	overlap_matrix.resize(vids.size());
	for (size_t i = 0; i < overlap_matrix.size(); i++)
		overlap_matrix[i].resize(vids.size());

	struct timeval start, end;
	gettimeofday(&start, NULL);

	overlap_stage = overlap_stage_t::CONSTRUCT_NEIGHBORS;
	graph->start(overlap_vertices.data(), overlap_vertices.size());
	graph->wait4complete();

	overlap_stage = overlap_stage_t::COMP_OVERLAP;
	graph->start(overlap_vertices.data(), overlap_vertices.size(),
			vertex_initializer::ptr(), vertex_program_creater::ptr(
				new overlap_vertex_program_creater(overlap_matrix)));
	graph->wait4complete();

	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("It takes %1% seconds to compute overlaps")
		% time_diff(start, end);
}
