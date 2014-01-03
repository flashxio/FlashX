#ifndef __EXT_VERTEX_H__
#define __EXT_VERTEX_H__

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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <vector>
#include <algorithm>

#include "container.h"
#include "cache.h"

typedef unsigned long vertex_id_t;

enum edge_type {
	NONE,
	IN_EDGE,
	OUT_EDGE,
	BOTH_EDGES,
};

/**
 * This class contains the basic information of a vertex in the memory.
 */
class in_mem_vertex_info
{
	vertex_id_t id;
	off_t off;
	int size;
public:
	in_mem_vertex_info() {
		off = 0;
		size = 0;
		id = -1;
	}

	in_mem_vertex_info(vertex_id_t id, off_t off, int size) {
		this->id = id;
		this->off = off;
		this->size = size;
	}

	off_t get_ext_mem_off() const {
		return off;
	}

	int get_ext_mem_size() const {
		return size;
	}

	vertex_id_t get_id() const {
		return id;
	}
};

class empty_data
{
};

template<class data_type = empty_data>
class edge
{
	bool has_data;
	vertex_id_t from;
	vertex_id_t to;
	data_type data;
public:
	edge() {
		this->from = -1;
		this->to = -1;
		has_data = false;
	}

	edge(vertex_id_t from, vertex_id_t to) {
		this->from = from;
		this->to = to;
		has_data = false;
	}

	edge(vertex_id_t from, vertex_id_t to, const data_type &data) {
		this->from = from;
		this->to = to;
		this->data = data;
		has_data = true;
	}

	vertex_id_t get_from() const {
		return from;
	}

	vertex_id_t get_to() const {
		return to;
	}

	bool has_edge_data() const {
		return has_data;
	}

	const data_type &get_data() const {
		assert(has_data);
		return data;
	}
};

template<class edge_data_type>
class in_mem_directed_vertex;

template<class edge_data_type>
class in_mem_undirected_vertex;

/**
 * This vertex represents a directed vertex stored in the external memory.
 */
class ext_mem_directed_vertex
{
	vertex_id_t id;
	int num_in_edges;
	int num_out_edges;
	vertex_id_t neighbors[0];

	void set_edge_data(bool has_data) {
		if (has_data) {
			// the max long is 0x80000000 in a 64-bit machine.
			id = id | LONG_MIN;
		}
		else {
			// the max long is 0x7FFFFFFF in a 64-bit machine.
			// This is more portable.
			id = id & LONG_MAX;
		}
	}

	bool has_edge_data() const {
		return id >> (sizeof(id) * 8 - 1);
	}

	void set_id(vertex_id_t id) {
		bool has_data = has_edge_data();
		this->id = id;
		set_edge_data(has_data);
	}

	template<class edge_data_type = empty_data>
	edge_data_type *get_edge_data_begin(edge_type type) {
		switch (type) {
			case edge_type::IN_EDGE:
				return (edge_data_type *) (neighbors + num_in_edges + num_out_edges);
			case edge_type::OUT_EDGE:
				return get_edge_data_begin<edge_data_type>(edge_type::IN_EDGE)
					+ num_in_edges;
			default:
				assert(0);
		}
	}
public:
	static ext_mem_directed_vertex *deserialize(char *buf, int size) {
		assert((unsigned) size >= sizeof(ext_mem_directed_vertex));
		ext_mem_directed_vertex *v = (ext_mem_directed_vertex *) buf;
		assert((unsigned) size >= sizeof(ext_mem_directed_vertex)
				+ (v->num_in_edges + v->num_out_edges) * sizeof(v->neighbors[0]));
		return v;
	}

	template<class edge_data_type = empty_data>
	static int serialize(const in_mem_directed_vertex<edge_data_type> &in_v,
			char *buf, int size) {
		int mem_size = in_v.get_serialize_size();
		assert(size >= mem_size);
		ext_mem_directed_vertex *ext_v = (ext_mem_directed_vertex *) buf;
		ext_v->set_id(in_v.get_id());
		ext_v->set_edge_data(in_v.has_edge_data());
		ext_v->num_in_edges = in_v.get_num_in_edges();
		ext_v->num_out_edges = in_v.get_num_out_edges();

		vertex_id_t *neighbors = ext_v->neighbors;
		for (int i = 0; i < in_v.get_num_in_edges(); i++) {
			neighbors[i] = in_v.get_in_edge(i).get_from();
		}
		for (int i = 0; i < in_v.get_num_out_edges(); i++) {
			neighbors[i + in_v.get_num_in_edges()] = in_v.get_out_edge(i).get_to();
		}

		// serialize edge data
		if (in_v.has_edge_data()) {
			for (int i = 0; i < in_v.get_num_in_edges(); i++) {
				edge<edge_data_type> e = in_v.get_in_edge(i);
				ext_v->get_edge_data_begin<edge_data_type>(edge_type::IN_EDGE)[i]
					= e.get_data();
			}
			for (int i = 0; i < in_v.get_num_out_edges(); i++) {
				edge<edge_data_type> e = in_v.get_out_edge(i);
				ext_v->get_edge_data_begin<edge_data_type>(edge_type::OUT_EDGE)[i]
					= e.get_data();
			}
		}
		return mem_size;
	}

	int get_num_edges(edge_type type) const {
		if (type == IN_EDGE)
			return get_num_in_edges();
		else if (type == OUT_EDGE)
			return get_num_out_edges();
		else
			return get_num_in_edges() + get_num_out_edges();
	}

	vertex_id_t get_neighbor(edge_type type, int idx) const {
		if (type == IN_EDGE)
			return neighbors[idx];
		else if (type == OUT_EDGE)
			return neighbors[num_in_edges + idx];
		else
			assert(0);
	}

	int get_num_in_edges() const {
		return num_in_edges;
	}

	int get_num_out_edges() const {
		return num_out_edges;
	}

	const vertex_id_t get_id() const {
		return id & LONG_MAX;
	}

	bool is_edge_list_sorted(edge_type type) const;

	int get_neighbors(edge_type type,
			fifo_queue<vertex_id_t> &neighbors) const {
		if (type == edge_type::IN_EDGE || type == edge_type::BOTH_EDGES) {
			int ret = neighbors.add((vertex_id_t *) this->neighbors, num_in_edges);
			assert(ret == num_in_edges);
		}
		if (type == edge_type::OUT_EDGE || type == edge_type::BOTH_EDGES) {
			int ret = neighbors.add((vertex_id_t *) this->neighbors + num_in_edges,
					num_out_edges);
			assert(ret == num_out_edges);
		}
		return neighbors.get_num_entries();
	}
};

/**
 * This vertex represents an undirected vertex in the external memory.
 */
class ext_mem_undirected_vertex
{
	vertex_id_t id;
	int num_edges;
	vertex_id_t neighbors[0];

	void set_edge_data(bool has_data) {
		if (has_data) {
			// the max long is 0x80000000 in a 64-bit machine.
			id = id | LONG_MIN;
		}
		else {
			// the max long is 0x7FFFFFFF in a 64-bit machine.
			// This is more portable.
			id = id & LONG_MAX;
		}
	}

	bool has_edge_data() const {
		return id >> (sizeof(id) * 8 - 1);
	}

	void set_id(vertex_id_t id) {
		bool has_data = has_edge_data();
		this->id = id;
		set_edge_data(has_data);
	}

	template<class edge_data_type = empty_data>
	edge_data_type *get_edge_data_begin() {
		return (edge_data_type *) (neighbors + num_edges);
	}
public:
	static ext_mem_undirected_vertex *deserialize(char *buf, int size) {
		assert((unsigned) size >= sizeof(ext_mem_undirected_vertex));
		ext_mem_undirected_vertex *v = (ext_mem_undirected_vertex *) buf;
		assert((unsigned) size >= sizeof(ext_mem_undirected_vertex)
				+ sizeof(v->neighbors[0]) * v->num_edges);
		return v;
	}

	template<class edge_data_type = empty_data>
	static int serialize(const in_mem_undirected_vertex<edge_data_type> &v,
			char *buf, int size) {
		int mem_size = v.get_serialize_size();
		assert(size >= mem_size);
		ext_mem_undirected_vertex *ext_v = (ext_mem_undirected_vertex *) buf;
		ext_v->set_id(v.get_id());
		ext_v->set_edge_data(v.has_edge_data());
		ext_v->num_edges = v.get_num_edges();

		vertex_id_t *neighbors = ext_v->neighbors;
		for (int i = 0; i < v.get_num_edges(); i++) {
			neighbors[i] = v.get_edge(i).get_to();
		}

		// serialize edge data
		if (v.has_edge_data()) {
			for (int i = 0; i < v.get_num_edges(); i++) {
				edge<edge_data_type> e = v.get_edge(i);
				ext_v->get_edge_data_begin<edge_data_type>()[i] = e.get_data();
			}
		}

		return mem_size;
	}

	int get_num_edges(edge_type type) const {
		return num_edges;
	}

	vertex_id_t get_neighbor(edge_type type, int idx) const {
		return neighbors[idx];
	}

	vertex_id_t get_id() const {
		return id & LONG_MAX;
	}

	bool is_edge_list_sorted(edge_type type) const;

	int get_neighbors(edge_type type,
			fifo_queue<vertex_id_t> &neighbors) const {
		int ret = neighbors.add((vertex_id_t *) this->neighbors, num_edges);
		assert(ret == num_edges);
		return neighbors.get_num_entries();
	}
};

class page_vertex
{
public:
	virtual int get_num_edges(edge_type type) const = 0;
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			edge_type type) const = 0;
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			edge_type type) const = 0;
	virtual vertex_id_t get_id() const = 0;
	bool contain_edge(edge_type type, vertex_id_t id) const {
		return std::binary_search(get_neigh_begin(type), get_neigh_end(type), id);
	}
};

/**
 * Time-series page vertex
 */
class TS_page_vertex: public page_vertex
{
public:
	virtual int get_num_edges(int timestamp, edge_type type) const = 0;
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			int timestamp, edge_type type) const = 0;
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			int timestamp, edge_type type) const = 0;
	bool contain_edge(int timestamp, edge_type type, vertex_id_t id) const {
		return std::binary_search(get_neigh_begin(timestamp, type),
				get_neigh_end(timestamp, type), id);
	}
};

/**
 * This vertex represents a directed vertex stored in the page cache.
 */
class page_directed_vertex: public page_vertex
{
	vertex_id_t id;
	int num_in_edges;
	int num_out_edges;
	const page_byte_array &array;
public:
	page_directed_vertex(const page_byte_array &arr): array(arr) {
		unsigned size = arr.get_size();
		// We only want to know the header of the vertex, so we don't need to
		// know what data type an edge has.
		assert((unsigned) size >= sizeof(ext_mem_directed_vertex));
		ext_mem_directed_vertex v = arr.get<ext_mem_directed_vertex>(0);
		assert((unsigned) size >= sizeof(ext_mem_directed_vertex)
				+ (v.get_num_in_edges() + v.get_num_out_edges()) * sizeof(vertex_id_t));

		id = v.get_id();
		num_in_edges = v.get_num_in_edges();
		num_out_edges = v.get_num_out_edges();
	}

	int get_num_edges(edge_type type) const {
		if (type == IN_EDGE)
			return num_in_edges;
		else if (type == OUT_EDGE)
			return num_out_edges;
		else if (type == BOTH_EDGES)
			return num_in_edges + num_out_edges;
		else
			assert(0);
	}

	page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			edge_type type) const {
		if (type == IN_EDGE || type == BOTH_EDGES)
			return array.begin<vertex_id_t>(sizeof(ext_mem_directed_vertex));
		else if (type == OUT_EDGE)
			return array.begin<vertex_id_t>(sizeof(ext_mem_directed_vertex)
				+ num_in_edges * sizeof(vertex_id_t));
		else
			assert(0);
	}

	page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			edge_type type) const {
		page_byte_array::const_iterator<vertex_id_t> it = get_neigh_begin(type);
		it += get_num_edges(type);
		return it;
	}

	vertex_id_t get_id() const {
		return id;
	}
};

/**
 * This vertex represents an undirected vertex in the page cache.
 */
class page_undirected_vertex: public page_vertex
{
	vertex_id_t id;
	int num_edges;
	const page_byte_array &array;
public:
	page_undirected_vertex(const page_byte_array &arr): array(arr) {
		unsigned size = arr.get_size();
		assert((unsigned) size >= sizeof(ext_mem_undirected_vertex));
		// We only want to know the header of the vertex, so we don't need to
		// know what data type an edge has.
		ext_mem_undirected_vertex v = arr.get<ext_mem_undirected_vertex>(0);
		assert((unsigned) size >= sizeof(ext_mem_undirected_vertex)
				+ sizeof(vertex_id_t) * v.get_num_edges(BOTH_EDGES));

		id = v.get_id();
		num_edges = v.get_num_edges(BOTH_EDGES);
	}

	int get_num_edges(edge_type type) const {
		return num_edges;
	}

	page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			edge_type type) const {
		return array.begin<vertex_id_t>(sizeof(ext_mem_undirected_vertex));
	}

	page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			edge_type type) const {
		page_byte_array::const_iterator<vertex_id_t> it = get_neigh_begin(type);
		it += num_edges;
		return it;
	}

	vertex_id_t get_id() const {
		return id;
	}
};

/**
 * This is a data structure to interpret a time-series directed vertex
 * in the external memory.
 * The memory layout of a time-series directed vertex is as follows:
 *	id
 *	the total number of edges
 *	edge offsets of each timestamp (in-edge offset and out-edge offset)
 *	edges
 *
 * The size of this object is defined at the runtime, so it can only be
 * allocated from the heap.
 */
class TS_page_directed_vertex: public TS_page_vertex
{
	struct edge_off
	{
		int in_off;
		int out_off;
	};

	struct TS_directed_vertex_header
	{
		vertex_id_t id;
		int num_edges;
	};

	vertex_id_t id;
	int num_timestamps;
	// The total number of edges in the vertex
	int num_edges;
	const page_byte_array &array;
	// The edge offsets of each timestamp in the edge list
	// They are offsets in edges (not in bytes).
	edge_off ts_edge_offs[0];

	TS_page_directed_vertex(int num_timestamps,
			const page_byte_array &arr): array(arr) {
		this->num_timestamps = num_timestamps;

		unsigned size = arr.get_size();
		assert((unsigned) size >= sizeof(TS_directed_vertex_header));
		TS_directed_vertex_header v = arr.get<TS_directed_vertex_header>(0);
		assert((unsigned) size >= sizeof(TS_directed_vertex_header)
				+ num_timestamps * sizeof(edge_off));
		assert((unsigned) size >= sizeof(TS_directed_vertex_header)
				+ num_timestamps * sizeof(edge_off)
				+ v.num_edges * sizeof(vertex_id_t));

		id = v.id;
		this->num_edges = v.num_edges;
		arr.memcpy(sizeof(TS_directed_vertex_header), (char *) ts_edge_offs,
				sizeof(edge_off) * num_timestamps);
	}

	// This object is not allowed to be copied.
	// Disable the copy constructor and the assign operator.
	TS_page_directed_vertex(const TS_page_directed_vertex &);
	TS_page_directed_vertex &operator=(const TS_page_directed_vertex &);
public:
	// We don't allow this object to be allocated in the stack.
	TS_page_directed_vertex *create(int num_timestamps,
			const page_byte_array &arr) {
		return NULL;
	}

	virtual int get_num_edges(edge_type type) const {
		assert(0);
	}

	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			edge_type type) const {
		assert(0);
	}

	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			edge_type type) const {
		assert(0);
	}

	virtual vertex_id_t get_id() const {
		return id;
	}

	virtual int get_num_edges(int timestamp, edge_type type) const {
		switch (type) {
			case edge_type::IN_EDGE:
				return ts_edge_offs[timestamp].out_off
					- ts_edge_offs[timestamp].in_off;
			case edge_type::OUT_EDGE:
				if (timestamp == num_timestamps - 1)
					return num_edges - ts_edge_offs[timestamp].out_off;
				else
					return ts_edge_offs[timestamp + 1].in_off
						- ts_edge_offs[timestamp].out_off;
			case edge_type::BOTH_EDGES:
				if (timestamp == num_timestamps - 1)
					return num_edges - ts_edge_offs[timestamp].in_off;
				else
					return ts_edge_offs[timestamp + 1].in_off
						- ts_edge_offs[timestamp].in_off;
			default:
				assert(0);
		}
	}

	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			int timestamp, edge_type type) const {
		// The start location of the edge list.
		page_byte_array::const_iterator<vertex_id_t> it
			= array.begin<vertex_id_t>(sizeof(TS_directed_vertex_header)
				+ num_timestamps * sizeof(edge_off));
		if (type == edge_type::IN_EDGE)
			it += ts_edge_offs[timestamp].in_off;
		else if (type == edge_type::OUT_EDGE)
			it += ts_edge_offs[timestamp].out_off;
		else
			assert(0);
		return it;
	}

	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			int timestamp, edge_type type) const {
		// The start location of the edge list.
		page_byte_array::const_iterator<vertex_id_t> it
			= array.begin<vertex_id_t>(sizeof(TS_directed_vertex_header)
				+ num_timestamps * sizeof(edge_off));
		if (type == edge_type::IN_EDGE)
			it += ts_edge_offs[timestamp].out_off;
		else if (type == edge_type::OUT_EDGE) {
			if (timestamp == num_timestamps - 1)
				it += num_edges;
			else
				it += ts_edge_offs[timestamp + 1].in_off;
		}
		else
			assert(0);
		return it;
	}
};

/**
 * This is the size of a page vertex (either directed or undirected).
 * It's mainly used for allocating a buffer from the stack for a page vertex.
 */
const int STACK_PAGE_VERTEX_SIZE = sizeof(page_directed_vertex);

template<class edge_data_type = empty_data>
class in_mem_directed_vertex
{
	vertex_id_t id;
	bool has_data;
	std::vector<vertex_id_t> out_edges;
	std::vector<vertex_id_t> in_edges;
	std::vector<edge_data_type> out_data;
	std::vector<edge_data_type> in_data;
public:
	in_mem_directed_vertex(vertex_id_t id) {
		this->id = id;
		has_data = false;
	}

	vertex_id_t get_id() const {
		return id;
	}

	bool has_edge_data() const {
		return has_data;
	}

	void add_in_edge(vertex_id_t id) {
		assert(!has_edge_data());
		in_edges.push_back(id);
	}

	void add_out_edge(vertex_id_t id) {
		assert(!has_edge_data());
		out_edges.push_back(id);
	}

	void add_in_edge(vertex_id_t id, const edge_data_type &data) {
		assert(has_edge_data());
		in_edges.push_back(id);
		in_data.push_back(data);
	}

	void add_out_edge(vertex_id_t id, const edge_data_type &data) {
		assert(has_edge_data());
		out_edges.push_back(id);
		out_data.push_back(data);
	}

	int get_num_in_edges() const {
		return in_edges.size();
	}

	int get_num_out_edges() const {
		return out_edges.size();
	}

	const edge<edge_data_type> get_in_edge(int idx) const {
		if (has_edge_data())
			return edge<edge_data_type>(in_edges[idx], id, in_data[idx]);
		else
			return edge<edge_data_type>(in_edges[idx], id);
	}

	const edge<edge_data_type> get_out_edge(int idx) const {
		if (has_edge_data())
			return edge<edge_data_type>(id, out_edges[idx], out_data[idx]);
		else
			return edge<edge_data_type>(id, out_edges[idx]);
	}

	int get_serialize_size() const {
		int size = sizeof(ext_mem_directed_vertex)
			+ sizeof(vertex_id_t) * (get_num_in_edges() + get_num_out_edges());
		if (has_edge_data())
			size += sizeof(edge_data_type) * (get_num_in_edges()
					+ get_num_out_edges());
		return size;
	}
};

template<class edge_data_type = empty_data>
class in_mem_undirected_vertex
{
	bool has_data;
	vertex_id_t id;
	std::vector<vertex_id_t> edges;
	std::vector<edge_data_type> data_arr;
public:
	in_mem_undirected_vertex(vertex_id_t id) {
		this->id = id;
		has_data = false;
	}

	vertex_id_t get_id() const {
		return id;
	}

	bool has_edge_data() const {
		return has_data;
	}

	void add_edge(vertex_id_t id) {
		assert(!has_edge_data());
		edges.push_back(id);
	}

	void add_edge(vertex_id_t id, const edge_data_type &data) {
		assert(has_edge_data());
		edges.push_back(id);
		data_arr.push_back(data);
	}

	int get_num_edges() const {
		return edges.size();
	}

	bool has_edge(vertex_id_t id) const {
		for (size_t i = 0; i < edges.size(); i++)
			if (edges[i] == id)
				return true;
		return false;
	}

	const edge<edge_data_type> get_edge(int idx) const {
		if (has_edge_data())
			return edge<edge_data_type>(id, edges[idx], data_arr[idx]);
		else
			return edge<edge_data_type>(id, edges[idx]);
	}

	int get_serialize_size() const {
		int size = sizeof(ext_mem_undirected_vertex)
			+ sizeof(vertex_id_t) * get_num_edges();
		if (has_edge_data())
			size += sizeof(edge_data_type) * get_num_edges();
		return size;
	}
};

#endif
