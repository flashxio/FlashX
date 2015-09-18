#ifndef __EXT_VERTEX_H__
#define __EXT_VERTEX_H__

/*
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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <memory>
#include <vector>
#include <algorithm>

#include "container.h"
#include "cache.h"
#include "FG_basic_types.h"

/**
 \brief Edge type of an edge in the graph.
 */
enum edge_type {
	NONE, /**No edge*/
	IN_EDGE, /**In edges*/
	OUT_EDGE, /**Out edges*/
	BOTH_EDGES, /**Both in and out edges*/
	NUM_TYPES,
};

class vertex_index;

/*
 * \brief This class contains the basic information about a vertex on the disk.
 *
 */
class ext_mem_vertex_info
{
	vertex_id_t id;
	vsize_t size;
	off_t off;
public:
	ext_mem_vertex_info() {
		id = INVALID_VERTEX_ID;
		off = 0;
		size = 0;
	}

	ext_mem_vertex_info(vertex_id_t id, off_t off, size_t size) {
		this->id = id;
		this->off = off;
		this->size = size;
	}

	vertex_id_t get_id() const {
		return id;
	}

	off_t get_off() const {
		return off;
	}

	vsize_t get_size() const {
		return size;
	}

	bool has_edges() const;

	bool is_valid() const {
		return id != INVALID_VERTEX_ID;
	}
};

/*
 * The information of vertex header.
 * It contains the vertex ID and the number of edges.
 */
class vertex_header
{
	vertex_id_t id;
	vsize_t num_edges;
public:
	vertex_header(vertex_id_t id, vsize_t num_edges) {
		this->id = id;
		this->num_edges = num_edges;
	}

	vertex_id_t get_id() const {
		return id;
	}

	vsize_t get_num_edges() const {
		return num_edges;
	}
};

/*
 * The information of directed vertex header.
 * In addition to the vertex header, it has the number of in-edges
 * and out-edges.
 */
class directed_vertex_header: public vertex_header
{
	vsize_t num_in_edges;
	vsize_t num_out_edges;
public:
	directed_vertex_header(vertex_id_t id, vsize_t num_in_edges,
			vsize_t num_out_edges): vertex_header(id,
				num_in_edges + num_out_edges) {
		this->num_in_edges = num_in_edges;
		this->num_out_edges = num_out_edges;
	}

	vsize_t get_num_in_edges() const {
		return num_in_edges;
	}

	vsize_t get_num_out_edges() const {
		return num_out_edges;
	}
};

class empty_data
{
public:
	bool operator==(const empty_data &data) const {
		return true;
	}
};

static inline std::ostream& operator<<(std::ostream& cout, empty_data obj)
{
	return cout;
}

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

	void reverse_dir() {
		vertex_id_t tmp = from;
		from = to;
		to = tmp;
	}
};

template<>
class edge<empty_data>
{
	static empty_data data;
	vertex_id_t from;
	vertex_id_t to;
public:
	edge() {
		this->from = -1;
		this->to = -1;
	}

	edge(vertex_id_t from, vertex_id_t to) {
		this->from = from;
		this->to = to;
	}

	edge(vertex_id_t from, vertex_id_t to, const empty_data &data) {
		this->from = from;
		this->to = to;
	}

	vertex_id_t get_from() const {
		return from;
	}

	vertex_id_t get_to() const {
		return to;
	}

	bool has_edge_data() const {
		return false;
	}

	const empty_data &get_data() const {
		return data;
	}

	void reverse_dir() {
		vertex_id_t tmp = from;
		from = to;
		to = tmp;
	}
};

/**
 * The timestamp of an edge.
 */
class ts_edge_data
{
	time_t timestamp;
public:
	ts_edge_data() {
		timestamp = 0;
	}

	ts_edge_data(time_t timestamp) {
		this->timestamp = timestamp;
	}

	time_t get_timestamp() const {
		return timestamp;
	}

	bool operator==(const ts_edge_data &data) const {
		return this->timestamp == data.timestamp;
	}

	bool operator<(const ts_edge_data &data) const {
		return this->timestamp < data.timestamp;
	}
};

static inline std::ostream& operator<<(std::ostream& cout,
		const ts_edge_data &obj)
{
	return cout << obj.get_timestamp();
}

template<>
class edge<ts_edge_data>
{
	vertex_id_t from;
	vertex_id_t to;
	ts_edge_data data;
public:
	edge() {
		this->from = -1;
		this->to = -1;
	}

	edge(vertex_id_t from, vertex_id_t to) {
		this->from = from;
		this->to = to;
	}

	edge(vertex_id_t from, vertex_id_t to, const ts_edge_data &data) {
		this->from = from;
		this->to = to;
		this->data = data;
	}

	vertex_id_t get_from() const {
		return from;
	}

	vertex_id_t get_to() const {
		return to;
	}

	bool has_edge_data() const {
		return true;
	}

	const ts_edge_data &get_data() const {
		return data;
	}

	void reverse_dir() {
		vertex_id_t tmp = from;
		from = to;
		to = tmp;
	}
};

/**
 * The number of duplicated edges.
 * It is used as edge data type.
 */
class edge_count
{
	uint32_t num;
public:
	edge_count() {
		num = 1;
	}

	edge_count(uint32_t n) {
		num = n;
	}

	uint32_t get_count() const {
		return num;
	}

	bool operator==(const edge_count &count) const {
		return this->num == count.num;
	}
};

static inline std::ostream& operator<<(std::ostream& cout,
		const edge_count &obj)
{
	return cout << obj.get_count();
}

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

template<class T>
struct delete_as_chararr
{
public:
	void operator()(T *obj) const {
		char *char_p = (char *) obj;
		delete [] char_p;
	}
};

class in_mem_vertex;

/*
 * This vertex represents an undirected vertex in the external memory.
 */
class ext_mem_undirected_vertex
{
	vertex_id_t id;
	uint32_t edge_data_size;
	vsize_t num_edges;
	vertex_id_t neighbors[0];

	void set_id(vertex_id_t id) {
		this->id = id;
	}

	/*
	 * The size of the vertex without counting the edge data list.
	 */
	size_t get_size0() const {
		return get_header_size() + sizeof(neighbors[0]) * num_edges;
	}

	char *get_edge_data_addr() const {
		return ((char *) this) + ROUNDUP(get_size0(), edge_data_size);
	}

	template<class edge_data_type = empty_data>
	const edge_data_type *get_edge_data_begin() const {
		assert(sizeof(edge_data_type) == edge_data_size);
		return (edge_data_type *) get_edge_data_addr();
	}

	template<class edge_data_type = empty_data>
	edge_data_type *get_edge_data_begin() {
		assert(sizeof(edge_data_type) == edge_data_size);
		return (edge_data_type *) get_edge_data_addr();
	}
public:
	static size_t get_header_size() {
		return offsetof(ext_mem_undirected_vertex, neighbors);
	}

	static size_t get_edge_data_offset(vsize_t num_edges,
			uint32_t edge_data_size) {
		ext_mem_undirected_vertex v(0, num_edges, edge_data_size);
		return v.get_edge_data_addr() - (char *) &v;
	}

	static ext_mem_undirected_vertex *deserialize(char *buf, size_t size) {
		assert(size >= ext_mem_undirected_vertex::get_header_size());
		ext_mem_undirected_vertex *v = (ext_mem_undirected_vertex *) buf;
		assert(size >= v->get_size());
		return v;
	}

	static size_t serialize(const in_mem_vertex &v, char *buf,
			size_t size, edge_type type);

	static vsize_t vsize2num_edges(size_t vertex_size, size_t edge_data_size) {
		return (vertex_size - get_header_size()) / (sizeof(vertex_id_t)
				+ edge_data_size);
	}

	static size_t num_edges2vsize(vsize_t num_edges, size_t edge_data_size) {
		ext_mem_undirected_vertex v(0, num_edges, edge_data_size);
		return v.get_size();
	}

	ext_mem_undirected_vertex() {
		this->id = 0;
		this->num_edges = 0;
		this->edge_data_size = 0;
	}

	ext_mem_undirected_vertex(vertex_id_t id, vsize_t num_edges,
			uint32_t edge_data_size) {
		this->id = id;
		this->num_edges = num_edges;
		this->edge_data_size = edge_data_size;
	}

	size_t get_size() const {
		if (has_edge_data())
			return ROUNDUP(((size_t) get_edge_data_addr()) - ((size_t) this)
				+ (num_edges) * edge_data_size, sizeof(vertex_id_t));
		else
			return ext_mem_undirected_vertex::get_header_size()
				+ (num_edges) * sizeof(neighbors[0]);
	}

	bool has_edge_data() const {
		return edge_data_size > 0;
	}

	size_t get_edge_data_size() const {
		return edge_data_size;
	}

	size_t get_num_edges() const {
		return num_edges;
	}

	vertex_id_t get_neighbor(size_t idx) const {
		return neighbors[idx];
	}

	template<class edge_data_type = empty_data>
	const edge_data_type &get_edge_data(size_t idx) const {
		return ((edge_data_type *) this->get_edge_data_addr())[idx];
	}

	vertex_id_t get_id() const {
		return id;
	}
};

inline bool ext_mem_vertex_info::has_edges() const
{
	return size > ext_mem_undirected_vertex::get_header_size();
}

/**
 * \brief Vertex representation when in the page cache.
 */
class page_vertex
{
	bool directed;
public:
	page_vertex(bool directed) {
		this->directed = directed;
	}
    /**
     * \brief Get the number of edges connecting the vertex to othe vertices.
     * \return The number of edges conning the vertex.
     * \param type The type of edges a user wishes to iterate 
                over e.g `IN_EDGE`, `OUT_EDGE`.
     */
	virtual size_t get_num_edges(edge_type type) const = 0;
    
    /**
     * \brief Get an STL-style const iterator pointing to the *first* neighbor 
            in a vertex's neighbor list.
     * \return A const iterator pointing to the *first* neighbor in 
            a vertex's neighbor list.
     * \param type The type of edges a user wishes to iterate over
            e.g `IN_EDGE`, `OUT_EDGE`.
     */
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			edge_type type) const = 0;
    
    /**
     * \brief Get an STL-style const iterator pointing to the *end* of
                    a vertex's neighbor list.
     * \return A const iterator pointing to the *end* of a vertex's neighbor list.
     * \param type The type of edges a user wishes to iterate over
                e.g `IN_EDGE`, `OUT_EDGE`.
     */
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			edge_type type) const = 0;
    /**
     * \brief Get a java-style sequential const iterator that iterates
	 *        the neighbors in the specified range.
     * \return A sequential const iterator.
     * \param type The type of edges a user wishes to iterate over e.g `IN_EDGE`, 
     *          `OUT_EDGE`.
	 * \param start The starting offset in the neighbor list iterated by
	 *              the sequential iterator.
	 * \param end The end offset in the neighbor list iterated by the sequential
	 *            iterator.
     */
	virtual page_byte_array::seq_const_iterator<vertex_id_t> get_neigh_seq_it(
			edge_type type, size_t start = 0, size_t end = -1) const = 0;
    
    /**
     * \brief Get the vertex unique ID.
     * \return The vertex unique ID.
     */
	virtual vertex_id_t get_id() const = 0;
    
    /**
     * \brief Read the edges of the specified type.
     * \param type The type of edges a user wishes to read
     *      e.g `IN_EDGE`, `OUT_EDGE`.
	 * \param edges The array of edges returned to a user.
	 * \param num The maximal number of edges read by a user.
     */
	virtual size_t read_edges(edge_type type, vertex_id_t edges[],
			size_t num) const {
		ABORT_MSG("read_edges isn't implemented");
		return 0;
	}

    /**
     * \brief Whether the vertex is directed.
	 * \return true if it's a directed vertex.
     */
	virtual bool is_directed() const {
		return directed;
	}
    
    /**
     * \internal
     */
	virtual void print() const {
	}
};

/*
 * These two ranges are defined as [first, second),
 * i.e., inclusive in the beginning and exclusive in the end.
 */
typedef std::pair<int, int> timestamp_pair;
typedef std::pair<off_t, off_t> offset_pair;

/**
 * Time-series page vertex utilized when doing time series graph analysis
 *
 */
class TS_page_vertex: public page_vertex
{
public:
	TS_page_vertex(bool directed): page_vertex(directed) {
	}

	using page_vertex::get_num_edges;
    
    /**
     * \brief Get the global number of edges associated with a vertex.
     * \return The number of edges associated with a vertex.
     */
	virtual size_t get_num_edges() const = 0;
    
    /**
     * \brief Get the number of edges associated with a vertex at a specific time point.
     * \param timestamp The specific time stamp where you want the vertex metadata evaluated.
     * \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`.
     * \return The number of edges associated with a vertex.
     */
	virtual size_t get_num_edges(int timestamp, edge_type type) const = 0;
    
    /**
     * \brief Get the number of time stamps the vertex has in the graph.
     * \return The number of time stamps the vertex has in the graph.
     */
	virtual int get_num_timestamps() const = 0;
	using page_vertex::get_neigh_begin;
	using page_vertex::get_neigh_end;
    
    /**
     * \brief Get an STL-style const iterator pointing to the *first* element in the
     *         neighbor list of a vertex at a specific time point.
     *  \param timpstamp The time stamp of interest.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`.
     *  \return A const iterator pointing to the *first* element in the
     *         neighbor list of a vertex.
     */
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			int timestamp, edge_type type) const = 0;
    
    /**
     * \brief Get an STL-style const iterator pointing to the *end* of the
     *         neighbor list of a vertex at a specific time point.
     *  \param timpstamp The time stamp of interest.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`.
     *  \return A const iterator pointing to the *end* of the
     *         neighbor list of a vertex.
     */
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			int timestamp, edge_type type) const = 0;

	/** \brief This method should translate the timestamp range to the absolute
     * location of the adjacency list in the timestamp range.
     *
     *  \param range The timestamp range.
	 *  \return The location range of the adjacency list within the time range.
     */
	virtual offset_pair get_edge_list_offset(
			const timestamp_pair &range) const = 0;
};

typedef page_byte_array::const_iterator<vertex_id_t> edge_iterator;
typedef page_byte_array::seq_const_iterator<vertex_id_t>  edge_seq_iterator;

/**
 * This vertex represents a directed vertex stored in the page cache.
 */
class page_directed_vertex: public page_vertex
{
	vertex_id_t id;
	vsize_t num_in_edges;
	vsize_t num_out_edges;
	size_t in_size;
	size_t out_size;
	const page_byte_array *in_array;
	const page_byte_array *out_array;
public:
	static vertex_id_t get_id(const page_byte_array &arr) {
		BOOST_VERIFY(arr.get_size()
				>= ext_mem_undirected_vertex::get_header_size());
		ext_mem_undirected_vertex v = arr.get<ext_mem_undirected_vertex>(0);
		return v.get_id();
	}
    
    /**
	 * \internal
	 * The constructor for a directed vertex in the page cache.
     *  \param arr The byte array containing the directed vertex
	 *             in the page cache.
     */
	page_directed_vertex(const page_byte_array &arr,
			bool in_part): page_vertex(true) {
		size_t size = arr.get_size();
		BOOST_VERIFY(size >= ext_mem_undirected_vertex::get_header_size());
		ext_mem_undirected_vertex v = arr.get<ext_mem_undirected_vertex>(0);

		if (in_part) {
			in_size = v.get_size();
			assert(size >= in_size);
			out_size = 0;
			this->in_array = &arr;
			this->out_array = NULL;
			num_in_edges = v.get_num_edges();
			num_out_edges = 0;
		}
		else {
			out_size = v.get_size();
			in_size = 0;
			assert(size >= out_size);
			this->out_array = &arr;
			this->in_array = NULL;
			num_out_edges = v.get_num_edges();
			num_in_edges = 0;
		}
		id = v.get_id();
	}

	page_directed_vertex(const page_byte_array &in_arr,
			const page_byte_array &out_arr): page_vertex(true) {
		this->in_array = &in_arr;
		this->out_array = &out_arr;

		size_t size = in_arr.get_size();
		BOOST_VERIFY(size >= ext_mem_undirected_vertex::get_header_size());
		ext_mem_undirected_vertex v = in_arr.get<ext_mem_undirected_vertex>(0);
		in_size = v.get_size();
		assert(size >= in_size);
		id = v.get_id();
		num_in_edges = v.get_num_edges();

		size = out_arr.get_size();
		assert(size >= ext_mem_undirected_vertex::get_header_size());
		v = out_arr.get<ext_mem_undirected_vertex>(0);
		out_size = v.get_size();
		assert(size >= out_size);
		assert(id == v.get_id());
		num_out_edges = v.get_num_edges();
	}

	size_t get_in_size() const {
		return in_size;
	}

	size_t get_out_size() const {
		return out_size;
	}
    
    /**
     * \brief Get the number of edges associated with a vertex.
     * \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `OUT_EDGE`.
     * \return The number of edges associated with a vertex.
     */
	size_t get_num_edges(edge_type type) const {
		switch(type) {
			case IN_EDGE:
				return num_in_edges;
			case OUT_EDGE:
				return num_out_edges;
			case BOTH_EDGES:
				return num_in_edges + num_out_edges;
			default:
				abort();
		}
	}
    
    /**
     * \brief Get an STL-style const iterator pointing to the *first* element in the
     *         neighbor list of a vertex.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `OUT_EDGE`.
     *  \return A const iterator pointing to the *first* element in the
     *         neighbor list of a vertex.
     */
	page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			edge_type type) const {
		switch(type) {
			case IN_EDGE:
				assert(in_array);
				return in_array->begin<vertex_id_t>(
						ext_mem_undirected_vertex::get_header_size());
			case OUT_EDGE:
				assert(out_array);
				return out_array->begin<vertex_id_t>(
						ext_mem_undirected_vertex::get_header_size());
			default:
				abort();
		}
	}
    
    /*
     * \brief Get an STL-style const iterator pointing to the *end* of the
     *         neighbor list of a vertex.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `OUT_EDGE`.
     *  \return A const iterator pointing to the *end* of the
     *         neighbor list of a vertex.
     */
	page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			edge_type type) const {
		page_byte_array::const_iterator<vertex_id_t> it = get_neigh_begin(type);
		it += get_num_edges(type);
		return it;
	}
    
    /*
     * \brief Get a java-style sequential const iterator that iterates
	 *        the neighbors in the specified range.
     * \return A sequential const iterator.
     * \param type The type of edges a user wishes to iterate over e.g `IN_EDGE`, 
     *          `OUT_EDGE`.
	 * \param start The starting offset in the neighbor list iterated by
	 *              the sequential iterator.
	 * \param end The end offset in the neighbor list iterated by the sequential
	 *            iterator.
     */
	page_byte_array::seq_const_iterator<vertex_id_t> get_neigh_seq_it(
			edge_type type, size_t start = 0, size_t end = -1) const {
		end = std::min(end, get_num_edges(type));
		assert(start <= end);
		switch(type) {
			case IN_EDGE:
				assert(in_array);
				return in_array->get_seq_iterator<vertex_id_t>(
						ext_mem_undirected_vertex::get_header_size()
						+ start * sizeof(vertex_id_t),
						ext_mem_undirected_vertex::get_header_size()
						+ end * sizeof(vertex_id_t));
			case OUT_EDGE:
				assert(out_array);
				return out_array->get_seq_iterator<vertex_id_t>(
						ext_mem_undirected_vertex::get_header_size()
						+ start * sizeof(vertex_id_t),
						ext_mem_undirected_vertex::get_header_size()
						+ end * sizeof(vertex_id_t));
			default:
				abort();
		}
	}
    
    /**
     * \brief Get an STL-style const iterator pointing to the *first* element in the
     *         edge data list of a vertex.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `OUT_EDGE`.
     *  \return A const iterator pointing to the *first* element in the
     *         edge data list of a vertex.
     */
	template<class edge_data_type>
	page_byte_array::const_iterator<edge_data_type> get_data_begin(
			edge_type type) const {
		switch(type) {
			case IN_EDGE:
				assert(in_array);
				return in_array->begin<edge_data_type>(
						ext_mem_undirected_vertex::get_edge_data_offset(
							num_in_edges, sizeof(edge_data_type)));
			case OUT_EDGE:
				assert(out_array);
				return out_array->begin<edge_data_type>(
						ext_mem_undirected_vertex::get_edge_data_offset(
							num_out_edges, sizeof(edge_data_type)));
			default:
				abort();
		}
	}

    /**
     * \brief Get an STL-style const iterator pointing to the *end* of the
     *         edge data list of a vertex.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `OUT_EDGE`.
     *  \return A const iterator pointing to the *first* element in the
     *         edge data list of a vertex.
     */
	template<class edge_data_type>
	page_byte_array::const_iterator<edge_data_type> get_data_end(
			edge_type type) const {
		page_byte_array::const_iterator<edge_data_type> it
			= get_data_begin<edge_data_type>(type);
		it += get_num_edges(type);
		return it;
	}

	/**
     * \brief Get a java-style sequential const iterator for edge data
	 *        with additional parameters to define the range to iterate.
     * \return A sequential const iterator.
     * \param type The type of edges a user wishes to iterate over e.g `IN_EDGE`, 
     *          `OUT_EDGE`, `OUT_EDGE`.
	 * \param start The starting offset in the edge list.
	 * \param end The end offset in the edge list.
	 */
	template<class edge_data_type>
	page_byte_array::seq_const_iterator<edge_data_type> get_data_seq_it(
			edge_type type, size_t start, size_t end) const {
		off_t edge_end;
		switch(type) {
			case IN_EDGE:
				assert(in_array);
				edge_end = ext_mem_undirected_vertex::get_edge_data_offset(
						num_in_edges, sizeof(edge_data_type));
				return in_array->get_seq_iterator<edge_data_type>(
						edge_end + start * sizeof(edge_data_type),
						edge_end + end * sizeof(edge_data_type));
			case OUT_EDGE:
				assert(out_array);
				edge_end = ext_mem_undirected_vertex::get_edge_data_offset(
						num_out_edges, sizeof(edge_data_type));
				return out_array->get_seq_iterator<edge_data_type>(
						edge_end + start * sizeof(edge_data_type),
						edge_end + end * sizeof(edge_data_type));
			default:
				abort();
		}
	}

    /**
     * \brief Get a java-style sequential const iterator that iterates
	 *        the edge data list.
     * \return A sequential const iterator.
     * \param type The type of edges a user wishes to iterate over e.g `IN_EDGE`, 
     *          `OUT_EDGE`, `OUT_EDGE`.
     */
	template<class edge_data_type>
	page_byte_array::seq_const_iterator<edge_data_type> get_data_seq_it(
			edge_type type) const {
		size_t start = 0;
		size_t end = get_num_edges(type);
		return get_data_seq_it<edge_data_type>(type, start, end);
	}
    
	virtual size_t read_edges(edge_type type, vertex_id_t edges[],
			size_t num) const {
		size_t num_edges;
		switch(type) {
			case IN_EDGE:
				assert(num_in_edges <= num);
				assert(in_array);
				num_edges = num_in_edges;
				in_array->memcpy(ext_mem_undirected_vertex::get_header_size(),
						(char *) edges, sizeof(vertex_id_t) * num_edges);
				break;
			case OUT_EDGE:
				assert(num_out_edges <= num);
				assert(out_array);
				num_edges = num_out_edges;
				out_array->memcpy(ext_mem_undirected_vertex::get_header_size(),
						(char *) edges, sizeof(vertex_id_t) * num_edges);
				break;
			default:
				abort();
		}
		return num_edges;
	}
    
    /** \brief Get the id of the vertex
     *  \return The vertex id
     */
	vertex_id_t get_id() const {
		return id;
	}

    /** \brief Determine whether the page vertex contains in-edges.
     * \return True if it contains in-edges.
     */
	bool has_in_part() const {
		return in_array;
	}

    /** \brief Determine whether the page vertex contains out-edges.
     * \return True if it contains out-edges.
     */
	bool has_out_part() const {
		return out_array;
	}
};

/**
 * \brief This vertex class represents an undirected vertex in the page cache.
 */
class page_undirected_vertex: public page_vertex
{
	vertex_id_t id;
	vsize_t vertex_size;
	vsize_t num_edges;
	const page_byte_array &array;
public:
	page_undirected_vertex(const page_byte_array &arr): page_vertex(
			false), array(arr) {
		size_t size = arr.get_size();
		BOOST_VERIFY(size >= ext_mem_undirected_vertex::get_header_size());
		// We only want to know the header of the vertex, so we don't need to
		// know what data type an edge has.
		ext_mem_undirected_vertex v = arr.get<ext_mem_undirected_vertex>(0);
		BOOST_VERIFY((unsigned) size >= v.get_size());
		vertex_size = v.get_size();

		id = v.get_id();
		num_edges = v.get_num_edges();
	}

	size_t get_size() const {
		return vertex_size;
	}
    
    /**
     * \brief Get the number of edges of a specific `edge_type` associated with the vertex.
     * \param type The type of edge i.e `IN_EDGE`, `OUT_EDGE` are equivalent,
	 *             since it's an undirected vertex.
     * return The number of edges associated with the vertex.
     */
	size_t get_num_edges(edge_type type = edge_type::IN_EDGE) const {
		return num_edges;
	}

	/**
	 * \brief Get an STL-style const iterator pointing to the *first* element in the
	 *         neighbor list of a vertex.
	 * \param type The type of edge i.e `IN_EDGE`, `OUT_EDGE` are equivalent,
	 *             since it's an undirected vertex.
	 * \return A const iterator pointing to the *first* element in the
	 *         neighbor list of a vertex.
	 */
	page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			edge_type type) const {
		return array.begin<vertex_id_t>(
				ext_mem_undirected_vertex::get_header_size());
	}

    /**
     * \brief Get an STL-style const iterator pointing to the *end* of the
     *         neighbor list of a vertex.
     * \param type The type of edge i.e `IN_EDGE`, `OUT_EDGE` are equivalent,
	 *             since it's an undirected vertex.
     *  \return A const iterator pointing to the *end* of the
     *         neighbor list of a vertex.
     */
	page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			edge_type type) const {
		page_byte_array::const_iterator<vertex_id_t> it = get_neigh_begin(type);
		it += num_edges;
		return it;
	}
    
    /**
     * \brief Get a java-style sequential const iterator for the specified range
     * in a vertex's neighbor list.
     * \param type The type of edge i.e `IN_EDGE`, `OUT_EDGE` are equivalent,
	 *             since it's an undirected vertex.
	 * \param start The starting offset in the neighbor list iterated by
	 *              the sequential iterator.
	 * \param end The end offset in the neighbor list iterated by the sequential
	 *            iterator.
     * \return A sequential const iterator for the specified range
     * in a vertex's neighbor list.
     */
	page_byte_array::seq_const_iterator<vertex_id_t> get_neigh_seq_it(
			edge_type type, size_t start = 0, size_t end = -1) const {
		end = std::min(end, get_num_edges(type));
		assert(start <= end);
		assert(end <= get_num_edges(type));
		return array.get_seq_iterator<vertex_id_t>(
				ext_mem_undirected_vertex::get_header_size()
				+ start * sizeof(vertex_id_t),
				ext_mem_undirected_vertex::get_header_size()
				+ end * sizeof(vertex_id_t));
	}
    
    /**
     * \brief Read the edges of the specified type.
     * \param type The type of edge i.e `IN_EDGE`, `OUT_EDGE` are equivalent,
	 *             since it's an undirected vertex.
	 * \param edges The array of edges returned to a user.
	 * \param num The maximal number of edges read by a user.
     */
	virtual size_t read_edges(edge_type type, vertex_id_t edges[],
			size_t num) const {
		vsize_t num_edges = get_num_edges(type);
		assert(num_edges <= num);
		array.memcpy(ext_mem_undirected_vertex::get_header_size(),
				(char *) edges, sizeof(vertex_id_t) * num_edges);
		return num_edges;
	}
    
    /**
     * \brief Get the vertex ID.
     * \return The vertex ID.
     */
	vertex_id_t get_id() const {
		return id;
	}
};

/*
 * The offset of in- and out-edges in the edge list of a time-series vertex.
 */
struct edge_off
{
	// TODO vsize_t may not be enough. A graph with many timestamps may have
	// many edges.
	vsize_t in_off;
	vsize_t out_off;
};

class in_mem_vertex
{
public:
	virtual vertex_id_t get_id() const = 0;
	virtual bool has_edge_data() const = 0;
	virtual size_t get_edge_data_size() const = 0;
	virtual void serialize_edges(vertex_id_t ids[], edge_type type) const = 0;
	virtual void serialize_edge_data(char *data, edge_type type) const = 0;
	virtual size_t get_serialize_size(edge_type type) const = 0;
	virtual size_t get_num_edges(edge_type type) const = 0;
};

/*
 * This is the size of a page vertex (either directed or undirected).
 * It's mainly used for allocating a buffer from the stack for a page vertex.
 */
const size_t STACK_PAGE_VERTEX_SIZE = sizeof(page_directed_vertex);

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

	in_mem_directed_vertex(const page_directed_vertex &vertex,
			bool has_data) {
		this->id = vertex.get_id();
		this->has_data = has_data;
		in_edges.resize(vertex.get_num_edges(edge_type::IN_EDGE));
		vertex.read_edges(edge_type::IN_EDGE, in_edges.data(),
				vertex.get_num_edges(edge_type::IN_EDGE));
		out_edges.resize(vertex.get_num_edges(edge_type::OUT_EDGE));
		vertex.read_edges(edge_type::OUT_EDGE, out_edges.data(),
				vertex.get_num_edges(edge_type::OUT_EDGE));
		assert(!has_data);
	}

	vertex_id_t get_id() const {
		return id;
	}

	bool has_edge_data() const {
		return has_data;
	}

	size_t get_edge_data_size() const {
		return has_data ? sizeof(edge_data_type) : 0;
	}

	virtual void serialize_edges(vertex_id_t ids[], edge_type type) const {
		switch (type) {
			case edge_type::IN_EDGE:
				memcpy(ids, in_edges.data(), in_edges.size() * sizeof(ids[0]));
				break;
			case edge_type::OUT_EDGE:
				memcpy(ids, out_edges.data(), out_edges.size() * sizeof(ids[0]));
				break;
			default:
				abort();
		}
	}

	virtual void serialize_edge_data(char *data, edge_type type) const {
		assert(has_data);
		switch (type) {
			case edge_type::IN_EDGE:
				memcpy(data, in_data.data(),
						in_data.size() * sizeof(edge_data_type));
				break;
			case edge_type::OUT_EDGE:
				memcpy(data, out_data.data(),
						out_data.size() * sizeof(edge_data_type));
				break;
			default:
				abort();
		}
	}

	/*
	 * Add an in-edge to the vertex.
	 * We allow to have multiple same edges, but all edges must be added
	 * in a sorted order.
	 */
	void add_in_edge(const edge<edge_data_type> &e) {
		assert(e.get_to() == id);
		in_edges.push_back(e.get_from());
		if (has_edge_data())
			in_data.push_back(e.get_data());
	}

	/*
	 * Add an out-edge to the vertex.
	 * We allow to have multiple same edges, but all edges must be added
	 * in a sorted order.
	 */
	void add_out_edge(const edge<edge_data_type> &e) {
		assert(e.get_from() == id);
		out_edges.push_back(e.get_to());
		if (has_edge_data())
			out_data.push_back(e.get_data());
	}

	size_t get_num_edges(edge_type type) const {
		switch(type) {
			case edge_type::IN_EDGE:
				return get_num_in_edges();
			case edge_type::OUT_EDGE:
				return get_num_out_edges();
			case edge_type::BOTH_EDGES:
				return get_num_in_edges() + get_num_out_edges();
			default:
				abort();
		}
	}

	size_t get_num_in_edges() const {
		return in_edges.size();
	}

	size_t get_num_out_edges() const {
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

	size_t get_serialize_size(edge_type type) const {
		assert(type == edge_type::IN_EDGE || type == edge_type::OUT_EDGE);
		if (type == edge_type::IN_EDGE) {
			ext_mem_undirected_vertex v(0, in_edges.size(),
					has_data ? sizeof(edge_data_type) : 0);
			return v.get_size();
		}
		else {
			ext_mem_undirected_vertex v(0, out_edges.size(),
					has_data ? sizeof(edge_data_type) : 0);
			return v.get_size();
		}
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
};

template<class edge_data_type = empty_data>
class in_mem_undirected_vertex: public in_mem_vertex
{
	bool has_data;
	vertex_id_t id;
	std::vector<vertex_id_t> edges;
	std::vector<edge_data_type> data_arr;
public:
	in_mem_undirected_vertex(vertex_id_t id, bool has_data) {
		this->id = id;
		this->has_data = has_data;
	}

	in_mem_undirected_vertex(const page_undirected_vertex &vertex,
			bool has_data) {
		this->id = vertex.get_id();
		this->has_data = has_data;
		edges.resize(vertex.get_num_edges());
		vertex.read_edges(edge_type::IN_EDGE, edges.data(),
				vertex.get_num_edges());
		assert(!has_data);
	}

	vertex_id_t get_id() const {
		return id;
	}

	bool has_edge_data() const {
		return has_data;
	}

	size_t get_edge_data_size() const {
		return has_data ? sizeof(edge_data_type) : 0;
	}

	virtual void serialize_edges(vertex_id_t ids[], edge_type type) const {
		memcpy(ids, edges.data(), edges.size() * sizeof(ids[0]));
	}

	virtual void serialize_edge_data(char *data, edge_type type) const {
		memcpy(data, data_arr.data(), data_arr.size() * sizeof(edge_data_type));
	}

	/*
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

	size_t get_num_edges(edge_type type = edge_type::IN_EDGE) const {
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

	size_t get_serialize_size(edge_type type) const {
		ext_mem_undirected_vertex v(0, edges.size(),
				has_data ? sizeof(edge_data_type) : 0);
		return v.get_size();
	}
};

#endif
