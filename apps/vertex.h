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

typedef unsigned int vertex_id_t;
const vertex_id_t MAX_VERTEX_ID = UINT_MAX;

enum edge_type {
	NONE,
	IN_EDGE,
	OUT_EDGE,
	BOTH_EDGES,
};

class vertex_index;

/**
 * This class contains the basic information of a vertex in the memory.
 */
class in_mem_vertex_info
{
	vertex_id_t id;
	int size;
	off_t off;
public:
	in_mem_vertex_info() {
		off = 0;
		size = 0;
		id = -1;
	}

	in_mem_vertex_info(vertex_id_t id, const vertex_index *index);

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

template<class edge_data_type>
class edge_const_iterator
{
	bool is_in_edge;
	vertex_id_t id;
	const vertex_id_t *ptr;
	const edge_data_type *data_ptr;
public:
	edge_const_iterator(vertex_id_t id,
			const vertex_id_t *edge_ptr,
			const edge_data_type *data_ptr, bool is_in_edge) {
		this->is_in_edge = is_in_edge;
		this->id = id;
		this->ptr = edge_ptr;
		this->data_ptr = data_ptr;
	}

	edge<edge_data_type> operator*() const {
		if (data_ptr) {
			if (is_in_edge)
				return edge<edge_data_type>(*ptr, id, *data_ptr);
			else
				return edge<edge_data_type>(id, *ptr, *data_ptr);
		}
		else {
			if (is_in_edge)
				return edge<edge_data_type>(*ptr, id);
			else
				return edge<edge_data_type>(id, *ptr);
		}
	}

	edge_const_iterator &operator++() {
		ptr++;
		if (data_ptr)
			data_ptr++;
		return *this;
	}

	bool operator==(const edge_const_iterator &it) const {
		return this->ptr == it.ptr;
	}

	bool operator!=(const edge_const_iterator &it) const {
		return this->ptr != it.ptr;
	}

	edge_const_iterator &operator+=(int num) {
		ptr += num;
		if (data_ptr)
			data_ptr += num;
		return *this;
	}
};

/**
 * This vertex represents a directed vertex stored in the external memory.
 */
class ext_mem_directed_vertex
{
	vertex_id_t id;
	int edge_data_size;
	int num_in_edges;
	int num_out_edges;
	vertex_id_t neighbors[0];

	bool has_edge_data() const {
		return edge_data_size > 0;
	}

	void set_id(vertex_id_t id) {
		this->id = id;
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
	static int get_header_size() {
		return offsetof(ext_mem_directed_vertex, neighbors);
	}

	static ext_mem_directed_vertex *deserialize(char *buf, int size) {
		assert((unsigned) size >= ext_mem_directed_vertex::get_header_size());
		ext_mem_directed_vertex *v = (ext_mem_directed_vertex *) buf;
		assert((unsigned) size >= ext_mem_directed_vertex::get_header_size()
				+ (v->num_in_edges + v->num_out_edges) * sizeof(v->neighbors[0]));
		return v;
	}

	static ext_mem_directed_vertex *merge(
			const std::vector<const ext_mem_directed_vertex *> &vertices,
			char *vertex_buf, size_t buf_size);

	template<class edge_data_type = empty_data>
	static int serialize(const in_mem_directed_vertex<edge_data_type> &in_v,
			char *buf, int size) {
		int mem_size = in_v.get_serialize_size();
		assert(size >= mem_size);
		ext_mem_directed_vertex *ext_v = (ext_mem_directed_vertex *) buf;
		ext_v->set_id(in_v.get_id());
		ext_v->num_in_edges = in_v.get_num_in_edges();
		ext_v->num_out_edges = in_v.get_num_out_edges();
		if (in_v.has_edge_data())
			ext_v->edge_data_size = sizeof(edge_data_type);
		else
			ext_v->edge_data_size = 0;

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

	int get_num_in_edges() const {
		return num_in_edges;
	}

	int get_num_out_edges() const {
		return num_out_edges;
	}

	const vertex_id_t get_id() const {
		return id;
	}

	template<class edge_data_type = empty_data>
	edge_const_iterator<edge_data_type> get_out_edge_begin() const {
		return edge_const_iterator<edge_data_type>(get_id(),
				neighbors + num_in_edges, NULL, false);
	}

	template<class edge_data_type = empty_data>
	edge_const_iterator<edge_data_type> get_out_edge_end() const {
		edge_const_iterator<edge_data_type> it
			= get_out_edge_begin<edge_data_type>();
		it += get_num_out_edges();
		return it;
	}

	template<class edge_data_type = empty_data>
	edge_const_iterator<edge_data_type> get_in_edge_begin() const {
		return edge_const_iterator<edge_data_type>(get_id(),
				neighbors, NULL, true);
	}

	template<class edge_data_type = empty_data>
	edge_const_iterator<edge_data_type> get_in_edge_end() const {
		edge_const_iterator<edge_data_type> it
			= get_in_edge_begin<edge_data_type>();
		it += get_num_in_edges();
		return it;
	}

	size_t get_size() const {
		size_t size = ext_mem_directed_vertex::get_header_size()
			+ sizeof(vertex_id_t) * (num_in_edges + num_out_edges);
		if (has_edge_data())
			size += edge_data_size * (num_in_edges + num_out_edges);
		return size;
	}
};

/**
 * This vertex represents an undirected vertex in the external memory.
 */
class ext_mem_undirected_vertex
{
	vertex_id_t id;
	int edge_data_size;
	int num_edges;
	vertex_id_t neighbors[0];

	bool has_edge_data() const {
		return edge_data_size > 0;
	}

	void set_id(vertex_id_t id) {
		this->id = id;
	}

	template<class edge_data_type = empty_data>
	edge_data_type *get_edge_data_begin() {
		return (edge_data_type *) (neighbors + num_edges);
	}
public:
	static int get_header_size() {
		return offsetof(ext_mem_undirected_vertex, neighbors);
	}

	static ext_mem_undirected_vertex *deserialize(char *buf, int size) {
		assert(size >= ext_mem_undirected_vertex::get_header_size());
		ext_mem_undirected_vertex *v = (ext_mem_undirected_vertex *) buf;
		assert((unsigned) size >= ext_mem_undirected_vertex::get_header_size()
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
		ext_v->num_edges = v.get_num_edges();
		if (v.has_edge_data())
			ext_v->edge_data_size = sizeof(edge_data_type);
		else
			ext_v->edge_data_size = 0;

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
		return id;
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
	virtual bool is_complete() const {
		return true;
	}
	bool contain_edge(edge_type type, vertex_id_t id) const {
		return std::binary_search(get_neigh_begin(type), get_neigh_end(type), id);
	}
	virtual void print() const {
	}
};

/**
 * These two ranges are defined as [first, second),
 * i.e., inclusive in the beginning and exclusive in the end.
 */
typedef std::pair<int, int> timestamp_pair;
typedef std::pair<off_t, off_t> offset_pair;

/**
 * Time-series page vertex
 */
class TS_page_vertex: public page_vertex
{
public:
	virtual int get_num_edges() const = 0;
	virtual int get_num_edges(int timestamp, edge_type type) const = 0;
	virtual int get_num_timestamps() const = 0;
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			int timestamp, edge_type type) const = 0;
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			int timestamp, edge_type type) const = 0;
	// This method should translate the timestamp range to the absolute
	// location of the adjacency list in the timestamp range.
	virtual offset_pair get_edge_list_offset(
			const timestamp_pair &range) const = 0;
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
		int size = arr.get_size();
		// We only want to know the header of the vertex, so we don't need to
		// know what data type an edge has.
		assert(size >= ext_mem_directed_vertex::get_header_size());
		ext_mem_directed_vertex v = arr.get<ext_mem_directed_vertex>(0);
		assert((unsigned) size >= ext_mem_directed_vertex::get_header_size()
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
			return array.begin<vertex_id_t>(
					ext_mem_directed_vertex::get_header_size());
		else if (type == OUT_EDGE)
			return array.begin<vertex_id_t>(
					ext_mem_directed_vertex::get_header_size()
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
		int size = arr.get_size();
		assert(size >= ext_mem_undirected_vertex::get_header_size());
		// We only want to know the header of the vertex, so we don't need to
		// know what data type an edge has.
		ext_mem_undirected_vertex v = arr.get<ext_mem_undirected_vertex>(0);
		assert((unsigned) size >= ext_mem_undirected_vertex::get_header_size()
				+ sizeof(vertex_id_t) * v.get_num_edges(BOTH_EDGES));

		id = v.get_id();
		num_edges = v.get_num_edges(BOTH_EDGES);
	}

	int get_num_edges(edge_type type) const {
		return num_edges;
	}

	page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			edge_type type) const {
		return array.begin<vertex_id_t>(
				ext_mem_undirected_vertex::get_header_size());
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

struct edge_off
{
	int in_off;
	int out_off;
};

template<class edge_data_type>
class ts_in_mem_directed_vertex;

/**
 * This class represents a time-series directed vertex in the external memory.
 * The memory layout of a time-series directed vertex is as follows:
 *	id
 *	the total number of edges
 *	the number of timestamps
 *	timestamps
 *	edge offsets of existing timestamp (in-edge offset and out-edge offset)
 *	edges
 *	edge data
 */
class ts_ext_mem_directed_vertex
{
	class edge_list {
		vertex_id_t *neighbors;
		int num_in_edges;
		int num_out_edges;
	public:
		edge_list(vertex_id_t *neighbors, int num_in_edges,
				int num_out_edges) {
			this->neighbors = neighbors;
			this->num_in_edges = num_in_edges;
			this->num_out_edges = num_out_edges;
		}

		int get_num_in_edges() const {
			return num_in_edges;
		}

		int get_num_out_edges() const {
			return num_out_edges;
		}

		vertex_id_t *get_in_edge_begin() {
			return neighbors;
		}

		vertex_id_t *get_out_edge_begin() {
			return neighbors + num_in_edges;
		}

		// The edge data is stored right at the end of the edge list.
		template<class edge_data_type>
		edge_data_type *get_in_edge_data() {
			return (edge_data_type *) (neighbors + num_in_edges + num_out_edges);
		}

		template<class edge_data_type>
		edge_data_type *get_out_edge_data() {
			return get_in_edge_data<edge_data_type>() + num_in_edges;
		}
	};

	class const_edge_list {
		const vertex_id_t *neighbors;
		int num_in_edges;
		int num_out_edges;
	public:
		static edge_list to_edge_list(const_edge_list &list) {
			return edge_list((vertex_id_t *) list.neighbors, list.num_in_edges,
					list.num_out_edges);
		}

		const_edge_list(const vertex_id_t *neighbors, int num_in_edges,
				int num_out_edges) {
			this->neighbors = neighbors;
			this->num_in_edges = num_in_edges;
			this->num_out_edges = num_out_edges;
		}

		int get_num_in_edges() const {
			return num_in_edges;
		}

		int get_num_out_edges() const {
			return num_out_edges;
		}

		const vertex_id_t *get_in_edge_begin() const {
			return neighbors;
		}

		const vertex_id_t *get_out_edge_begin() const {
			return neighbors + num_in_edges;
		}

		// The edge data is stored right at the end of the edge list.
		template<class edge_data_type>
		const edge_data_type *get_in_edge_data() const {
			return (edge_data_type *) (neighbors + num_in_edges + num_out_edges);
		}

		template<class edge_data_type>
			const edge_data_type *get_out_edge_data() const {
				return get_in_edge_data<edge_data_type>() + num_in_edges;
			}
	};

	vertex_id_t id;
	int edge_data_size;
	int num_edges;
	int num_timestamps;

	static int get_timestamp_list_size(int num_timestamps) {
		// The edge off list need to align with the word size.
		int num = ROUNDUP(num_timestamps, 8);
		return sizeof(short) * num;
	}

	static int get_timestamp_table_size(int num_timestamps) {
		return get_timestamp_list_size(num_timestamps)
			+ sizeof(edge_off) * num_timestamps;
	}

	void set_edge_data_size(int size) {
		this->edge_data_size = size;
	}

	bool has_edge_data() const {
		return edge_data_size > 0;
	}

	void set_id(vertex_id_t id) {
		this->id = id;
	}

	int find_timestamp(int timestamp) const {
		for (int i = 0; i < num_timestamps; i++)
			if (get_timestamps_begin()[i] == timestamp)
				return i;
		return -1;
	}

	short *get_timestamps_begin() {
		return (short *) (this + 1);
	}

	const short *get_timestamps_begin() const {
		return (short *) (this + 1);
	}

	short get_first_timestamp() const {
		assert(num_timestamps > 0);
		return get_timestamps_begin()[0];
	}

	short get_last_timestamp() const {
		assert(num_timestamps > 0);
		return get_timestamps_begin()[num_timestamps - 1];
	}

	edge_off *get_edge_off_begin() {
		char *p = (char *) get_timestamps_begin();
		return (edge_off *) (p + get_timestamp_list_size(num_timestamps));
	}

	const edge_off *get_edge_off_begin() const {
		const char *p = (const char *) get_timestamps_begin();
		return (const edge_off *) (p + get_timestamp_list_size(num_timestamps));
	}

	vertex_id_t *get_edge_list_begin() {
		char *p = (char *) get_timestamps_begin();
		return (vertex_id_t *) (p + get_timestamp_table_size(num_timestamps));
	}

	const vertex_id_t *get_edge_list_begin() const {
		const char *p = (const char *) get_timestamps_begin();
		return (const vertex_id_t *) (p
				+ get_timestamp_table_size(num_timestamps));
	}

	// The header size of the ts-vertex in the external memory.
	int get_header_size() const {
		return sizeof(ts_ext_mem_directed_vertex)
			+ get_timestamp_table_size(num_timestamps);
	}

	int get_num_in_edges_idx(int idx) const {
		if (idx < 0)
			return 0;
		else
			return get_edge_off_begin()[idx].out_off
				- get_edge_off_begin()[idx].in_off;
	}

	int get_num_out_edges_idx(int idx) const {
		if (idx < 0)
			return 0;

		// The last timestamp
		if (idx == num_timestamps - 1)
			return num_edges - get_edge_off_begin()[idx].out_off;
		else {
			return get_edge_off_begin()[idx + 1].in_off
				- get_edge_off_begin()[idx].out_off;
		}
	}

	/**
	 * These functions define the data layout of the edge list
	 * in the external memory.
	 */

	edge_list get_edge_list(int timestamp) {
		const_edge_list list = get_const_edge_list(timestamp);
		return const_edge_list::to_edge_list(list);
	}

	const_edge_list get_const_edge_list_idx(int idx) const {
		assert(idx < num_timestamps);
		char *edge_list_start = (char *) get_edge_list_begin();
		int num_prev_edges = get_edge_off_begin()[idx].in_off;
		int prev_edge_list_size = num_prev_edges * (sizeof(vertex_id_t)
				+ edge_data_size);
		return const_edge_list((vertex_id_t *) (edge_list_start
					+ prev_edge_list_size), get_num_in_edges_idx(idx),
				get_num_out_edges_idx(idx));
	}

	const_edge_list get_const_edge_list(int timestamp) const {
		int idx = find_timestamp(timestamp);
		if (idx >= 0) {
			return get_const_edge_list_idx(idx);
		}
		else {
			char *edge_list_start = (char *) get_edge_list_begin();
			int prev_edge_list_size = num_edges * (sizeof(vertex_id_t)
					+ edge_data_size);
			return const_edge_list((vertex_id_t *) (edge_list_start
						+ prev_edge_list_size), 0, 0);
		}
	}

	/**
	 * This returns the offset (in bytes) of the edge list of the specified
	 * timestamp in the external memory.
	 * If the timestamp doesn't exist, return the offset of the end
	 * of the vertex.
	 */
	int get_edge_list_offset(int timestamp, edge_type type) const {
		const_edge_list edges = get_const_edge_list(timestamp);
		if (type == edge_type::IN_EDGE)
			return ((char *) edges.get_in_edge_begin()) - ((char *) this);
		else if (type == edge_type::OUT_EDGE)
			return ((char *) edges.get_out_edge_begin()) - ((char *) this);
		else
			assert(0);
	}

	/**
	 * This construct a header of the vertex based on the original vertex
	 * header passed to the method.
	 * @edge_list_off: the relative location (in bytes) of the edge lists
	 * in the vertex.
	 * @edge_list_size: the size (in bytes) of the edge lists.
	 */
	void construct_header(const ts_ext_mem_directed_vertex &header,
			int edge_list_off, int edge_list_size);

	/**
	 * This returns the offset (in bytes) of the edge data list of the specified
	 * timestamp in the external memory.
	 * If the timestamp doesn't exist, return the offset of the end
	 * of the vertex.
	 */
	template<class edge_data_type>
	int get_edge_data_offset(int timestamp, edge_type type) const {
		assert(sizeof(edge_data_type) == edge_data_size);
		const_edge_list edges = get_const_edge_list(timestamp);
		if (type == edge_type::IN_EDGE)
			return ((char *) edges.get_in_edge_data<edge_data_type>())
				- ((char *) this);
		else if (type == edge_type::OUT_EDGE)
			return ((char *) edges.get_out_edge_data<edge_data_type>())
				- ((char *) this);
		else
			assert(0);
	}

	ts_ext_mem_directed_vertex(vertex_id_t id, int num_edges,
			int num_timestamps, int edge_data_size) {
		this->id = id;
		this->num_edges = num_edges;
		this->num_timestamps = num_timestamps;
		this->edge_data_size = edge_data_size;
	}
public:
	template<class edge_data_type = empty_data>
	static int serialize(const ts_in_mem_directed_vertex<edge_data_type> &in_v,
			char *buf, int size) {
		ts_ext_mem_directed_vertex *v = new (buf) ts_ext_mem_directed_vertex(
				in_v.get_id(), in_v.get_num_edges(), in_v.get_num_timestamps(),
				in_v.has_edge_data() ? sizeof(edge_data_type) : 0);
		assert(v->get_size() <= (size_t) size);

		// Generate the timestamp table.
		// It's likely that many vertices don't have edges in most of timestamps.
		std::vector<int> all_timestamps;
		in_v.get_all_timestamps(all_timestamps);
		assert(v->num_timestamps == (int) all_timestamps.size());
		int off = 0;
		for (int i = 0; i < v->num_timestamps; i++) {
			int timestamp = all_timestamps[i];
			v->get_timestamps_begin()[i] = timestamp;
			v->get_edge_off_begin()[i].in_off = off;
			off += in_v.get_num_in_edges(timestamp);
			v->get_edge_off_begin()[i].out_off = off;
			off += in_v.get_num_out_edges(timestamp);
		}
		assert(v->num_edges == off);

		// Generate the edge list.
		for (int i = 0; i < v->num_timestamps; i++) {
			int timestamp = all_timestamps[i];
			int num_in_edges = in_v.get_num_in_edges(timestamp);
			edge_list edges = v->get_edge_list(timestamp);
			edge_const_iterator<edge_data_type> it
				= in_v.get_in_edge_begin(timestamp);
			assert(edges.get_num_in_edges() == num_in_edges);
			for (int j = 0; j < num_in_edges; j++) {
				edge<edge_data_type> e = *it;
				++it;
				edges.get_in_edge_begin()[j] = e.get_from();
				if (v->has_edge_data())
					edges.get_in_edge_data<edge_data_type>()[j] = e.get_data();
			}

			int num_out_edges = in_v.get_num_out_edges(timestamp);
			it = in_v.get_out_edge_begin(timestamp);
			assert(edges.get_num_out_edges() == num_out_edges);
			for (int j = 0; j < num_out_edges; j++) {
				edge<edge_data_type> e = *it;
				++it;
				edges.get_out_edge_begin()[j] = e.get_to();
				if (v->has_edge_data())
					edges.get_out_edge_data<edge_data_type>()[j] = e.get_data();
			}
		}
		return v->get_size();
	}

	static ts_ext_mem_directed_vertex *merge(
			const std::vector<ext_mem_directed_vertex *> &vertices,
			char *buf, size_t size);

	ts_ext_mem_directed_vertex() {
		this->id = 0;
		this->num_edges = 0;
		this->num_timestamps = 0;
		this->edge_data_size = 0;
	}

	vertex_id_t get_id() const {
		return id;
	}

	size_t get_size() const {
		size_t size = sizeof(ts_ext_mem_directed_vertex)
			+ get_timestamp_table_size(num_timestamps)
			+ sizeof(vertex_id_t) * num_edges;
		if (has_edge_data())
			size += edge_data_size * num_edges;
		return size;
	}

	int get_num_edges() const {
		return num_edges;
	}

	int get_num_timestamps() const {
		return num_timestamps;
	}

	int get_edge_data_size() const {
		return edge_data_size;
	}

	int get_num_in_edges(int timestamp) const {
		int idx = find_timestamp(timestamp);
		return get_num_in_edges_idx(idx);
	}

	int get_num_out_edges(int timestamp) const {
		int idx = find_timestamp(timestamp);
		return get_num_out_edges_idx(idx);
	}

	template<class edge_data_type = empty_data>
	edge_const_iterator<edge_data_type> get_out_edge_begin(int timestamp) const {
		const_edge_list edges = get_const_edge_list(timestamp);
		const edge_data_type *data_ptr = NULL;
		const vertex_id_t *ptr = edges.get_out_edge_begin();
		if (has_edge_data())
			data_ptr = edges.get_out_edge_data<edge_data_type>();
		return edge_const_iterator<edge_data_type>(get_id(),
				ptr, data_ptr, false);
	}

	template<class edge_data_type = empty_data>
	edge_const_iterator<edge_data_type> get_out_edge_end(int timestamp) const {
		edge_const_iterator<edge_data_type> it
			= get_out_edge_begin<edge_data_type>(timestamp);
		it += get_num_out_edges(timestamp);
		return it;
	}

	template<class edge_data_type = empty_data>
	edge_const_iterator<edge_data_type> get_in_edge_begin(int timestamp) const {
		const_edge_list edges = get_const_edge_list(timestamp);
		const edge_data_type *data_ptr = NULL;
		const vertex_id_t *ptr = edges.get_in_edge_begin();
		if (has_edge_data())
			data_ptr = edges.get_in_edge_data<edge_data_type>();
		return edge_const_iterator<edge_data_type>(get_id(),
				ptr, data_ptr, true);
	}

	template<class edge_data_type = empty_data>
	edge_const_iterator<edge_data_type> get_in_edge_end(int timestamp) const {
		edge_const_iterator<edge_data_type> it
			= get_in_edge_begin<edge_data_type>(timestamp);
		it += get_num_in_edges(timestamp);
		return it;
	}

	int get_num_in_edges() const {
		int num = 0;
		// TODO I should use a simpler way to compute the result.
		for (int i = 0; i < num_timestamps; i++)
			num += get_num_in_edges(get_timestamps_begin()[i]);
		return num;
	}

	int get_num_out_edges() const {
		int num = 0;
		// TODO I should use a simpler way to compute the result.
		for (int i = 0; i < num_timestamps; i++)
			num += get_num_out_edges(get_timestamps_begin()[i]);
		return num;
	}

	void print() const {
		printf("v%ld has edge data: %d, # timestamps: %d, # edges: %d\n",
				(unsigned long) get_id(), 0, num_timestamps, get_num_edges());
		for (int i = 0; i < num_timestamps; i++) {
			int timestamp = get_timestamps_begin()[i];
			// We need to skip the timestamps without edges.
			if (get_num_in_edges(timestamp) + get_num_out_edges(timestamp) == 0)
				continue;

			printf("timestamp %d\n", timestamp);
			int num_in_edges = get_num_in_edges(timestamp);
			printf("in-edges (%d): ", num_in_edges);
			const vertex_id_t *in_edge_list = get_edge_list_begin()
				+ get_edge_off_begin()[i].in_off;
			for (int j = 0; j < num_in_edges; j++) {
				printf("%ld, ", (unsigned long) in_edge_list[j]);
			}
			printf("\n");
			int num_out_edges = get_num_out_edges(timestamp);
			printf("out-edges (%d): ", num_out_edges);
			const vertex_id_t *out_edge_list = get_edge_list_begin()
				+ get_edge_off_begin()[i].out_off;
			for (int j = 0; j < num_out_edges; j++) {
				printf("%ld, ", (unsigned long) out_edge_list[j]);
			}
			printf("\n");
		}
	}

	template<class edge_data_type>
	friend class ts_in_mem_directed_vertex;
	friend class TS_page_directed_vertex;
};

/**
 * This is a data structure to interpret a time-series directed vertex
 * in the external memory.
 * The memory layout of a time-series directed vertex is as follows:
 *	id
 *	the total number of edges
 *	the number of timestamps
 *	timestamps
 *	edge offsets of existing timestamp (in-edge offset and out-edge offset)
 *	edges
 *	edge data
 *
 * The size of this object is defined at the runtime, so it can only be
 * allocated from the heap.
 */
class TS_page_directed_vertex: public TS_page_vertex
{
	// The entire size of the vertex in the external memory.
	int entire_vertex_size;
	// The location of the vertex in the file.
	off_t vertex_begin_off;
	const page_byte_array *array;
	ts_ext_mem_directed_vertex ext_v;

	TS_page_directed_vertex(const page_byte_array &arr): array(&arr) {
		unsigned size = arr.get_size();
		assert((unsigned) size >= sizeof(ts_ext_mem_directed_vertex));
		ts_ext_mem_directed_vertex v = arr.get<ts_ext_mem_directed_vertex>(0);
		entire_vertex_size = v.get_size();
		vertex_begin_off = arr.get_offset();
		// We have to make sure the array size must be larger than the header
		// of the ts-vertex.
		int header_size = v.get_header_size();
		assert(arr.get_size() >= header_size);
		assert(header_size <= PAGE_SIZE);
		arr.memcpy(0, (char *) &ext_v, header_size);
		// If the vertex isn't complete, we probably only has the header
		// of the vertex. Don't give a uaser a chance to access it.
		if (!is_complete())
			array = NULL;
	}

	/**
	 * This is to create a vertex based on the header of the vertex
	 * read in the last time and the edge list in the array.
	 * The constructed vertex has only part of the vertex.
	 */
	TS_page_directed_vertex(const TS_page_directed_vertex *header,
			const page_byte_array &arr): array(&arr) {
		this->entire_vertex_size = header->entire_vertex_size;
		this->vertex_begin_off = header->vertex_begin_off;
		off_t rel_off = arr.get_offset() - vertex_begin_off;
		assert(rel_off >= header->ext_v.get_header_size());
		ext_v.construct_header(header->ext_v, rel_off, arr.get_size());
	}

	// This object is not allowed to be copied.
	// Disable the copy constructor and the assign operator.
	TS_page_directed_vertex(const TS_page_directed_vertex &);
	TS_page_directed_vertex &operator=(const TS_page_directed_vertex &);

	page_byte_array::const_iterator<vertex_id_t> get_neigh_end() const {
		assert(array);
		return array->begin<vertex_id_t>(array->get_size());
	}

	template<class edge_data_type>
	page_byte_array::const_iterator<edge_data_type> get_edge_data_end() const {
		assert(array);
		return array->begin<edge_data_type>(array->get_size());
	}
public:
	// The size of the vertex object.
	static int get_size(int num_timestamps) {
		return sizeof(TS_page_directed_vertex)
			+ ts_ext_mem_directed_vertex::get_timestamp_table_size(num_timestamps);
	}

	// We create the vertex object in the given buffer.
	static TS_page_directed_vertex *create(const page_byte_array &arr,
			char *buf, int size) {
		return new (buf) TS_page_directed_vertex(arr);
	}

	static TS_page_directed_vertex *create(const TS_page_directed_vertex *header,
			const page_byte_array &arr, char *buf, int size) {
		return new (buf) TS_page_directed_vertex(header, arr);
	}

	/**
	 * This returns the relative offset of the edge lists specified
	 * in the timestamp range in the vertex.
	 */
	virtual offset_pair get_edge_list_offset(const timestamp_pair &range) const;

	virtual bool is_complete() const {
		if (array == NULL)
			return false;
		else
			return array->get_size() >= entire_vertex_size;
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
		return ext_v.get_id();
	}

	virtual int get_num_edges() const {
		return ext_v.get_num_edges();
	}

	virtual int get_num_timestamps() const {
		return ext_v.get_num_timestamps();
	}

	virtual int get_num_edges(int timestamp, edge_type type) const {
		switch (type) {
			case edge_type::IN_EDGE:
				return ext_v.get_num_in_edges(timestamp);
			case edge_type::OUT_EDGE:
				return ext_v.get_num_out_edges(timestamp);
			case edge_type::BOTH_EDGES:
				return ext_v.get_num_in_edges(timestamp)
					+ ext_v.get_num_out_edges(timestamp);
			default:
				assert(0);
		}
	}

	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			int timestamp, edge_type type) const {
		int offset = 0;
		// The start location of the edge list.
		if (type == edge_type::IN_EDGE || type == edge_type::BOTH_EDGES)
			offset = ext_v.get_edge_list_offset(timestamp, edge_type::IN_EDGE);
		else if (type == edge_type::OUT_EDGE)
			offset = ext_v.get_edge_list_offset(timestamp, edge_type::OUT_EDGE);
		else
			assert(0);
		assert(array);
		if (offset < 0)
			return get_neigh_end();
		else {
			if (array->get_size() == entire_vertex_size)
				return array->begin<vertex_id_t>(offset);
			else
				// If the vertex isn't complete, the array only contains
				// the data of edge lists. We need to substract the header
				// size from the offset.
				return array->begin<vertex_id_t>(offset - ext_v.get_header_size());
		}
	}

	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			int timestamp, edge_type type) const {
		// The start location of the edge list.
		page_byte_array::const_iterator<vertex_id_t> it
			= get_neigh_begin(timestamp, type);
		// because the timestamp doesn't exist.
		if (it == get_neigh_end())
			return it;

		it += get_num_edges(timestamp, type);
		return it;
	}

	template<class edge_data_type>
	page_byte_array::const_iterator<edge_data_type> get_edge_data_begin(
			int timestamp, edge_type type) const {
		int offset = 0;
		if (type == edge_type::IN_EDGE || type == edge_type::BOTH_EDGES)
			offset = ext_v.get_edge_data_offset<edge_data_type>(timestamp,
					edge_type::IN_EDGE);
		else if (type == edge_type::OUT_EDGE)
			offset = ext_v.get_edge_data_offset<edge_data_type>(timestamp,
					edge_type::OUT_EDGE);
		else
			assert(0);
		assert(array);
		if (offset < 0)
			return get_edge_data_end<edge_data_type>();
		else {
			if (array->get_size() == entire_vertex_size)
				return array->begin<edge_data_type>(offset);
			else
				// If the vertex isn't complete, the array only contains
				// the data of edge lists. We need to substract the header
				// size from the offset.
				return array->begin<edge_data_type>(offset - ext_v.get_header_size());
		}
	}

	template<class edge_data_type>
	page_byte_array::const_iterator<edge_data_type> get_edge_data_end(
			int timestamp, edge_type type) const {
		// The start location of the edge list.
		page_byte_array::const_iterator<edge_data_type> it
			= get_edge_data_begin<edge_data_type>(timestamp, type);
		// because the timestamp doesn't exist.
		if (it == get_edge_data_end<edge_data_type>())
			return it;

		it += get_num_edges(timestamp, type);
		return it;
	}

	virtual void print() const {
		printf("v%ld has edge data: %d, # timestamps: %d, # edges: %d\n",
				(unsigned long) get_id(), 0, get_num_timestamps(), get_num_edges());
		for (int i = 0; i < get_num_timestamps(); i++) {
			int timestamp = ext_v.get_timestamps_begin()[i];
			// We need to skip the timestamps without edges.
			if (get_num_edges(timestamp, edge_type::IN_EDGE)
					+ get_num_edges(timestamp, edge_type::OUT_EDGE) == 0)
				continue;

			printf("timestamp %d\n", timestamp);
			int num_in_edges = get_num_edges(timestamp, edge_type::IN_EDGE);
			printf("in-edges (%d): ", num_in_edges);
			page_byte_array::const_iterator<vertex_id_t> end_it
				= get_neigh_end(timestamp, edge_type::IN_EDGE);
			for (page_byte_array::const_iterator<vertex_id_t> it
					= get_neigh_begin(timestamp, edge_type::IN_EDGE);
					it != end_it; ++it) {
				printf("%ld, ", (unsigned long) *it);
			}
			printf("\n");
			int num_out_edges = get_num_edges(timestamp, edge_type::OUT_EDGE);
			printf("out-edges (%d): ", num_out_edges);
			end_it = get_neigh_end(timestamp, edge_type::OUT_EDGE);
			for (page_byte_array::const_iterator<vertex_id_t> it
					= get_neigh_begin(timestamp, edge_type::OUT_EDGE);
					it != end_it; ++it) {
				printf("%ld, ", (unsigned long) *it);
			}
			printf("\n");
		}
	}
};

class in_mem_vertex
{
public:
	virtual vertex_id_t get_id() const = 0;
	virtual bool has_edge_data() const = 0;
	virtual int get_serialize_size() const = 0;
	virtual int get_num_edges(edge_type type) const = 0;
};

/**
 * This is the size of a page vertex (either directed or undirected).
 * It's mainly used for allocating a buffer from the stack for a page vertex.
 */
const int STACK_PAGE_VERTEX_SIZE = sizeof(page_directed_vertex);

template<class edge_data_type = empty_data>
class in_mem_directed_vertex: public in_mem_vertex
{
	vertex_id_t id;
	bool has_data;
	std::vector<vertex_id_t> out_edges;
	std::vector<vertex_id_t> in_edges;
	std::vector<edge_data_type> out_data;
	std::vector<edge_data_type> in_data;
public:
	in_mem_directed_vertex(vertex_id_t id, bool has_data) {
		this->id = id;
		this->has_data = has_data;
	}

	vertex_id_t get_id() const {
		return id;
	}

	bool has_edge_data() const {
		return has_data;
	}

	/**
	 * Add an in-edge to the vertex.
	 * We allow to have multiple same edges, but all edges must be added
	 * in a sorted order.
	 */
	void add_in_edge(const edge<edge_data_type> &e) {
		assert(e.get_to() == id);
		if (!in_edges.empty())
			assert(e.get_from() >= in_edges.back());
		in_edges.push_back(e.get_from());
		if (has_edge_data())
			in_data.push_back(e.get_data());
	}

	/**
	 * Add an out-edge to the vertex.
	 * We allow to have multiple same edges, but all edges must be added
	 * in a sorted order.
	 */
	void add_out_edge(const edge<edge_data_type> &e) {
		assert(e.get_from() == id);
		if (!out_edges.empty())
			assert(e.get_to() >= out_edges.back());
		out_edges.push_back(e.get_to());
		if (has_edge_data())
			out_data.push_back(e.get_data());
	}

	int get_num_edges(edge_type type) const {
		switch(type) {
			case edge_type::IN_EDGE:
				return get_num_in_edges();
			case edge_type::OUT_EDGE:
				return get_num_out_edges();
			case edge_type::BOTH_EDGES:
				return get_num_in_edges() + get_num_out_edges();
			default:
				assert(0);
		}
	}

	int get_num_in_edges() const {
		return in_edges.size();
	}

	int get_num_out_edges() const {
		return out_edges.size();
	}

	edge_const_iterator<edge_data_type> get_in_edge_begin() const {
		return edge_const_iterator<edge_data_type>(id,
				in_edges.data(), in_data.data(), true);
	}

	edge_const_iterator<edge_data_type> get_in_edge_end() const {
		edge_const_iterator<edge_data_type> it = get_in_edge_begin();
		it += get_num_in_edges();
		return it;
	}

	edge_const_iterator<edge_data_type> get_out_edge_begin() const {
		return edge_const_iterator<edge_data_type>(id,
				out_edges.data(), out_data.data(), false);
	}

	edge_const_iterator<edge_data_type> get_out_edge_end() const {
		edge_const_iterator<edge_data_type> it = get_out_edge_begin();
		it += get_num_out_edges();
		return it;
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
		int size = ext_mem_directed_vertex::get_header_size()
			+ sizeof(vertex_id_t) * (get_num_in_edges() + get_num_out_edges());
		if (has_edge_data())
			size += sizeof(edge_data_type) * (get_num_in_edges()
					+ get_num_out_edges());
		return size;
	}

	void print() const {
		printf("v%ld has edge data: %d\n", (unsigned long) get_id(), has_edge_data());
		printf("There are %ld in-edges: ", in_edges.size());
		for (size_t i = 0; i < in_edges.size(); i++)
			printf("%ld, ", (unsigned long) in_edges[i]);
		printf("\n");
		printf("There are %ld out-edges: ", out_edges.size());
		for (size_t i = 0; i < out_edges.size(); i++)
			printf("%ld, ", (unsigned long) out_edges[i]);
		printf("\n");
	}

	friend class ts_in_mem_directed_vertex<edge_data_type>;
};

template<class edge_data_type = empty_data>
class in_mem_undirected_vertex: public in_mem_vertex
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

	/**
	 * Add an edge to the vertex.
	 * We allow to have multiple same edges, but all edges must be added
	 * in a sorted order.
	 */
	void add_edge(const edge<edge_data_type> &e) {
		assert(e.get_from() == id);
		if (!edges.empty())
			assert(e.get_to() >= edges.back());
		edges.push_back(e.get_to());
		if (has_edge_data())
			data_arr.push_back(e.get_data());
	}

	int get_num_edges(edge_type type = edge_type::BOTH_EDGES) const {
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
		int size = ext_mem_undirected_vertex::get_header_size()
			+ sizeof(vertex_id_t) * get_num_edges();
		if (has_edge_data())
			size += sizeof(edge_data_type) * get_num_edges();
		return size;
	}
};

template<class edge_data_type = empty_data>
class ts_in_mem_directed_vertex: public in_mem_vertex
{
	struct ts_edge_pair {
		std::vector<vertex_id_t> in_edges;
		std::vector<vertex_id_t> out_edges;
	};

	struct ts_edge_data_pair {
		std::vector<edge_data_type> in_data;
		std::vector<edge_data_type> out_data;
	};

	vertex_id_t id;
	bool has_data;
	std::map<int, ts_edge_pair> ts_edges;
	std::map<int, ts_edge_data_pair> ts_data;
public:
	ts_in_mem_directed_vertex(vertex_id_t id, bool has_data) {
		this->id = id;
		this->has_data = has_data;
	}

	vertex_id_t get_id() const {
		return id;
	}

	bool has_edge_data() const {
		return has_data;
	}

	int get_num_in_edges(int timestamp) const {
		typename std::map<int, ts_edge_pair>::const_iterator it
			= ts_edges.find(timestamp);
		if (it == ts_edges.end())
			return 0;
		else
			return it->second.in_edges.size();
	}

	int get_num_out_edges(int timestamp) const {
		typename std::map<int, ts_edge_pair>::const_iterator it
			= ts_edges.find(timestamp);
		if (it == ts_edges.end())
			return 0;
		else
			return it->second.out_edges.size();
	}

	int get_num_edges(edge_type type) const {
		switch(type) {
			case edge_type::IN_EDGE:
				return get_num_in_edges();
			case edge_type::OUT_EDGE:
				return get_num_out_edges();
			case edge_type::BOTH_EDGES:
				return get_num_in_edges() + get_num_out_edges();
			default:
				assert(0);
		}
	}

	int get_num_in_edges() const {
		int num_edges = 0;
		for (typename std::map<int, ts_edge_pair>::const_iterator it
				= ts_edges.begin(); it != ts_edges.end(); it++)
			num_edges += it->second.in_edges.size();
		return num_edges;
	}

	int get_num_out_edges() const {
		int num_edges = 0;
		for (typename std::map<int, ts_edge_pair>::const_iterator it
				= ts_edges.begin(); it != ts_edges.end(); it++)
			num_edges += it->second.out_edges.size();
		return num_edges;
	}

	int get_num_edges() const {
		return get_num_in_edges() + get_num_out_edges();
	}

	int get_num_timestamps() const {
		return ts_edges.size();
	}

	void get_all_timestamps(std::vector<int> &timestamps) const {
		for (typename std::map<int, ts_edge_pair>::const_iterator it
				= ts_edges.begin(); it != ts_edges.end(); it++) {
			timestamps.push_back(it->first);
		}
	}

	edge_const_iterator<edge_data_type> get_in_edge_begin(
			int timestamp) const {
		typename std::map<int, ts_edge_pair>::const_iterator it = ts_edges.find(
				timestamp);
		typename std::map<int, ts_edge_data_pair>::const_iterator data_it
			= ts_data.find(timestamp);
		return edge_const_iterator<edge_data_type>(id,
				it->second.in_edges.data(), data_it->second.in_data.data(),
				true);
	}

	edge_const_iterator<edge_data_type> get_in_edge_end(
			int timestamp) const {
		edge_const_iterator<edge_data_type> it
			= get_in_edge_begin(timestamp);
		it += get_num_in_edges(timestamp);
		return it;
	}

	edge_const_iterator<edge_data_type> get_out_edge_begin(
			int timestamp) const {
		typename std::map<int, ts_edge_pair>::const_iterator it = ts_edges.find(
				timestamp);
		typename std::map<int, ts_edge_data_pair>::const_iterator data_it
			= ts_data.find(timestamp);
		return edge_const_iterator<edge_data_type>(id,
				it->second.out_edges.data(), data_it->second.out_data.data(),
				false);
	}

	edge_const_iterator<edge_data_type> get_out_edge_end(
			int timestamp) const {
		edge_const_iterator<edge_data_type> it
			= get_out_edge_begin(timestamp);
		it += get_num_out_edges(timestamp);
		return it;
	}

	void add_timestamp(int timestamp,
			const in_mem_directed_vertex<edge_data_type> &v) {
		assert(id == v.get_id());
		assert(has_data == v.has_edge_data());
		assert(ts_edges.find(timestamp) == ts_edges.end());
		assert(v.get_num_in_edges() + v.get_num_out_edges() > 0);
		if (has_data) {
			if (!v.in_edges.empty())
				assert(!v.in_data.empty());
			if (!v.out_edges.empty())
				assert(!v.out_data.empty());
		}

		ts_edge_pair edge_pair;
		edge_pair.in_edges = v.in_edges;
		edge_pair.out_edges = v.out_edges;

		ts_edge_data_pair data_pair;
		data_pair.in_data = v.in_data;
		data_pair.out_data = v.out_data;

		ts_edges.insert(std::pair<int, ts_edge_pair>(timestamp, edge_pair));
		ts_data.insert(std::pair<int, ts_edge_data_pair>(timestamp, data_pair));
	}

	int get_serialize_size() const {
		ts_ext_mem_directed_vertex v(get_id(), get_num_edges(),
				get_num_timestamps(),
				has_edge_data() ? sizeof(edge_data_type) : 0);
		return v.get_size();
	}

	void print() const {
		printf("v%ld has edge data: %d, # timestamps: %d, # edges: %d\n",
				(unsigned long) get_id(), has_edge_data(),
				get_num_timestamps(), get_num_edges());
		for (typename std::map<int, ts_edge_pair>::const_iterator it
				= ts_edges.begin(); it != ts_edges.end(); it++) {
			printf("timestamp %d\n", it->first);
			printf("in-edges (%d): ", get_num_in_edges(it->first));
			for (size_t i = 0; i < it->second.in_edges.size(); i++) {
				printf("%ld, ", (unsigned long) it->second.in_edges[i]);
			}
			printf("\n");
			printf("out-edges (%d): ", get_num_out_edges(it->first));
			for (size_t i = 0; i < it->second.out_edges.size(); i++) {
				printf("%ld, ",  (unsigned long)it->second.out_edges[i]);
			}
			printf("\n");
		}
	}
};

/**
 * The number of duplicated edges.
 * It is used as edge data type.
 */
class edge_count
{
	int num;
public:
	edge_count() {
		num = 1;
	}

	edge_count(int n) {
		num = n;
	}

	int get_count() const {
		return num;
	}
};

#endif
