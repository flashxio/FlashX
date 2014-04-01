#ifndef __SCAN_GRAPH_H__
#define __SCAN_GRAPH_H__

/**
 * Copyright 2014 Da Zheng
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

#include <memory>

#include "graphlab/cuckoo_set_pow2.hpp"
#include "graph_engine.h"

/**
 * The edge has two attributes:
 * The number of duplicated edges with the neighbor;
 * The number of edges that contributes to the neighbor's neighborhood.
 */
class attributed_neighbor
{
	vertex_id_t id;
	int num_dups;
public:
	attributed_neighbor() {
		id = -1;
		num_dups = 0;
	}

	attributed_neighbor(vertex_id_t id) {
		this->id = id;
		num_dups = 1;
	}

	attributed_neighbor(vertex_id_t id, int num_dups) {
		this->id = id;
		this->num_dups = num_dups;
	}

	vertex_id_t get_id() const {
		return id;
	}

	int get_num_dups() const {
		return num_dups;
	}

	bool operator<(const attributed_neighbor &e) const {
		return this->id < e.id;
	}

	bool operator==(vertex_id_t id) const {
		return this->id == id;
	}

	bool operator==(const attributed_neighbor &e) const {
		return this->id == e.id;
	}
};

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

class neighbor_list
{
public:
	class index_entry
	{
		vertex_id_t id;
		uint32_t idx;
	public:
		index_entry() {
			id = -1;
			idx = -1;
		}

		index_entry(vertex_id_t id) {
			this->id = id;
			this->idx = -1;
		}

		index_entry(vertex_id_t id, uint32_t idx) {
			this->id = id;
			this->idx = idx;
		}

		vertex_id_t get_id() const {
			return id;
		}

		uint32_t get_idx() const {
			return idx;
		}

		bool operator==(const index_entry &e) const {
			return id == e.get_id();
		}
	};

	class index_hash
	{
		boost::hash<vertex_id_t> id_hash;
	public:
		size_t operator()(const index_entry &e) const {
			return id_hash(e.get_id());
		}
	};

	typedef graphlab::cuckoo_set_pow2<index_entry, 3, size_t,
			index_hash> edge_set_t;

protected:
	// The vertex Id that the neighbor list belongs to.
	vertex_id_t id;
	std::vector<vertex_id_t> id_list;
	std::vector<int> num_dup_list;
	edge_set_t *neighbor_set;
public:
	class id_iterator: public std::iterator<std::random_access_iterator_tag, vertex_id_t>
	{
		std::vector<vertex_id_t>::const_iterator it;
	public:
		typedef typename std::iterator<std::random_access_iterator_tag,
				vertex_id_t>::difference_type difference_type;

		id_iterator() {
		}

		id_iterator(const std::vector<vertex_id_t> &v) {
			it = v.begin();
		}

		difference_type operator-(const id_iterator &it) const {
			return this->it - it.it;
		}

		vertex_id_t operator*() const {
			return *it;
		}

		id_iterator &operator++() {
			it++;
			return *this;
		}

		id_iterator operator++(int) {
			id_iterator ret = *this;
			it++;
			return ret;
		}

		bool operator==(const id_iterator &it) const {
			return it.it == this->it;
		}
		
		bool operator!=(const id_iterator &it) const {
			return it.it != this->it;
		}

		id_iterator &operator+=(size_t num) {
			it += num;
			return *this;
		}
	};

	neighbor_list(const page_vertex &vertex,
			const std::vector<attributed_neighbor> &neighbors) {
		this->id = vertex.get_id();
		size_t num_neighbors = neighbors.size();
		id_list.resize(num_neighbors);
		num_dup_list.resize(num_neighbors);
		for (size_t i = 0; i < num_neighbors; i++) {
			id_list[i] = neighbors[i].get_id();
			num_dup_list[i] = neighbors[i].get_num_dups();
		}
		neighbor_set = NULL;
		if (num_neighbors > 0) {
			neighbor_set = new edge_set_t(index_entry(),
					0, 2 * neighbors.size());
			for (size_t i = 0; i < neighbors.size(); i++)
				neighbor_set->insert(index_entry(neighbors[i].get_id(), i));
		}
	}

	virtual ~neighbor_list() {
		if (neighbor_set)
			delete neighbor_set;
	}

	vertex_id_t get_neighbor_id(size_t idx) const {
		return id_list[idx];
	}

	bool contains(vertex_id_t id) {
		edge_set_t::const_iterator it = neighbor_set->find(id);
		return it != neighbor_set->end();
	}

	id_iterator get_id_begin() const {
		return id_iterator(id_list);
	}

	id_iterator get_id_end() const {
		id_iterator ret(id_list);
		ret += id_list.size();
		return ret;
	}

	size_t size() const {
		return id_list.size();
	}

	bool empty() const {
		return id_list.empty();
	}

	size_t get_neighbors(std::vector<vertex_id_t> &neighbors) {
		neighbors = id_list;
		return neighbors.size();
	}

	/**
	 * This return the vertex Id that the neighbor list belongs to.
	 */
	vertex_id_t get_id() const {
		return id;
	}

	virtual size_t count_edges(const page_vertex *v);
	virtual size_t count_edges(const page_vertex *v, edge_type type,
			std::vector<vertex_id_t> *common_neighs);
	virtual size_t count_edges_hash(const page_vertex *v,
			page_byte_array::const_iterator<vertex_id_t> other_it,
			page_byte_array::const_iterator<vertex_id_t> other_end,
			std::vector<vertex_id_t> *common_neighs);
#if 0
	virtual size_t count_edges_bin_search_this(const page_vertex *v,
			neighbor_list::id_iterator this_it,
			neighbor_list::id_iterator this_end,
			page_byte_array::const_iterator<vertex_id_t> other_it,
			page_byte_array::const_iterator<vertex_id_t> other_end);
#endif
	virtual size_t count_edges_bin_search_other(const page_vertex *v,
			neighbor_list::id_iterator this_it,
			neighbor_list::id_iterator this_end,
			page_byte_array::const_iterator<vertex_id_t> other_it,
			page_byte_array::const_iterator<vertex_id_t> other_end,
			std::vector<vertex_id_t> *common_neighs);
	virtual size_t count_edges_scan(const page_vertex *v,
			neighbor_list::id_iterator this_it,
			neighbor_list::id_iterator this_end,
			page_byte_array::seq_const_iterator<vertex_id_t> other_it,
			std::vector<vertex_id_t> *common_neighs);
};

/**
 * This data structure contains all data required when the vertex is
 * computing local scan.
 * It doesn't need to exist before or after the vertex computes local scan.
 */
struct runtime_data_t
{
	// All neighbors (in both in-edges and out-edges)
	std::unique_ptr<neighbor_list> neighbors;
	// The number of vertices that have joined with the vertex.
	unsigned num_joined;
	size_t local_scan;

	runtime_data_t(std::unique_ptr<neighbor_list> neighbors) {
		this->neighbors = std::move(neighbors);
		num_joined = 0;
		local_scan = 0;
	}
};

enum multi_func_flags
{
	EST_LOCAL,
	REAL_LOCAL,
	POINTER,
	NUM_FLAGS,
};

class multi_func_value
{
	static const int VALUE_BITS = sizeof(size_t) * 8 - NUM_FLAGS;
	static const size_t FLAGS_MASK = ((1UL << VALUE_BITS) - 1);
	size_t value;

	void set_flag(int flag) {
		value |= 1UL << (VALUE_BITS + flag);
	}

	bool has_flag(int flag) const {
		return value & (1UL << (VALUE_BITS + flag));
	}
public:
	multi_func_value() {
		value = 0;
	}

	void set_est_local(size_t num) {
		value = num;
		set_flag(EST_LOCAL);
	}

	bool has_est_local() const {
		return has_flag(EST_LOCAL);
	}

	size_t get_est_local() const {
		assert(has_flag(EST_LOCAL));
		return value & FLAGS_MASK;
	}

	void set_real_local(size_t num) {
		value = num;
		set_flag(REAL_LOCAL);
	}

	void inc_real_local(size_t num) {
		assert(REAL_LOCAL);
		value += num;
	}

	bool has_real_local() const {
		return has_flag(REAL_LOCAL);
	}

	size_t get_real_local() const {
		assert(has_flag(REAL_LOCAL));
		return value & FLAGS_MASK;
	}
};

class scan_vertex: public compute_vertex
{
protected:
	multi_func_value local_value;

#ifdef PV_STAT
	// For testing
	size_t num_all_edges;
	size_t scan_bytes;
	size_t rand_jumps;
	size_t min_comps;
	long time_us;
	struct timeval vertex_start;
#endif
public:
	runtime_data_t *data;

	scan_vertex() {
		data = NULL;

#ifdef PV_STAT
		num_all_edges = 0;
		scan_bytes = 0;
		rand_jumps = 0;
		min_comps = 0;
		time_us = 0;
#endif
	}

	scan_vertex(vertex_id_t id, const vertex_index *index): compute_vertex(
			id, index) {
		data = NULL;

#ifdef PV_STAT
		num_all_edges = 0;
		scan_bytes = 0;
		rand_jumps = 0;
		min_comps = 0;
		time_us = 0;
#endif
	}

#ifdef PV_STAT
	size_t get_scan_bytes() const {
		return scan_bytes;
	}

	size_t get_rand_jumps() const {
		return rand_jumps;
	}
#endif

	bool has_local_scan() const {
		return local_value.has_real_local();
	}

	size_t get_local_scan() const {
		return local_value.get_real_local();
	}

	void run(graph_engine &graph, const page_vertex &vertex) {
		if (vertex.get_id() == get_id())
			run_on_itself(graph, vertex);
		else
			run_on_neighbor(graph, vertex);
	}

	virtual void run_on_itself(graph_engine &graph, const page_vertex &vertex);
	void run_on_neighbor(graph_engine &graph, const page_vertex &vertex);
	virtual void finding_triangles_end(graph_engine &graph) = 0;

	virtual runtime_data_t *create_runtime(graph_engine &graph,
			const page_vertex *vertex);
	virtual void destroy_runtime(runtime_data_t *data);

	void run_on_messages(graph_engine &graph,
			const vertex_message *msgs[], int num) {
	}
};

extern atomic_number<long> num_working_vertices;
extern atomic_number<long> num_completed_vertices;

#endif
