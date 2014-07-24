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
};

class vertex_index;

/*
 * \brief This class contains the basic information about a vertex in the memory.
 *
 */
class in_mem_vertex_info
{
	vsize_t size;
	off_t off;
public:
	in_mem_vertex_info() {
		off = 0;
		size = 0;
	}

	in_mem_vertex_info(vertex_id_t id, off_t off, size_t size) {
		this->off = off;
		this->size = size;
	}

	off_t get_ext_mem_off() const {
		return off;
	}

	vsize_t get_ext_mem_size() const {
		return size;
	}
};

class in_mem_directed_vertex_info: public in_mem_vertex_info
{
	vsize_t num_in_edges;
	vsize_t num_out_edges;
public:
	in_mem_directed_vertex_info() {
		num_in_edges = 0;
		num_out_edges = 0;
	}

	in_mem_directed_vertex_info(vertex_id_t id, off_t off, size_t size,
			vsize_t num_in_edges, vsize_t num_out_edges): in_mem_vertex_info(id,
				off, size) {
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

/**
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

/**
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

class ts_ext_mem_directed_vertex;

template<class T>
struct delete_as_chararr
{
public:
	void operator()(T *obj) const {
		char *char_p = (char *) obj;
		delete [] char_p;
	}
};

/*
 * This vertex represents a directed vertex stored in the external memory.
 */
class ext_mem_directed_vertex
{
	vertex_id_t id;
	uint32_t edge_data_size;
	vsize_t num_in_edges;
	vsize_t num_out_edges;
	vertex_id_t neighbors[0];

	void set_id(vertex_id_t id) {
		this->id = id;
	}

	char *get_edge_data_addr() const {
		return (char *) ROUNDUP((long) (neighbors + num_in_edges
					+ num_out_edges), edge_data_size);
	}

	template<class edge_data_type = empty_data>
	const edge_data_type *get_edge_data_begin(edge_type type) const {
		return ((ext_mem_directed_vertex *) this)->get_edge_data_begin<edge_data_type>(type);
	}

	template<class edge_data_type = empty_data>
	edge_data_type *get_edge_data_begin(edge_type type) {
		assert(sizeof(edge_data_type) == edge_data_size);
		switch (type) {
			case edge_type::IN_EDGE:
				return (edge_data_type *) get_edge_data_addr();
			case edge_type::OUT_EDGE:
				return ((edge_data_type *) get_edge_data_addr()) + num_in_edges;
			default:
				assert(0);
		}
	}
public:
	typedef std::unique_ptr<ext_mem_directed_vertex,
			delete_as_chararr<ext_mem_directed_vertex> > unique_ptr;

	static size_t get_header_size() {
		return offsetof(ext_mem_directed_vertex, neighbors);
	}

	static off_t get_edge_data_offset(vsize_t num_in_edges,
			vsize_t num_out_edges, uint32_t edge_data_size) {
		ext_mem_directed_vertex v(0, num_in_edges, num_out_edges,
				edge_data_size);
		return v.get_edge_data_addr() - (char *) &v;
	}

	static ext_mem_directed_vertex *deserialize(char *buf, size_t size) {
		assert(size >= ext_mem_directed_vertex::get_header_size());
		ext_mem_directed_vertex *v = (ext_mem_directed_vertex *) buf;
		assert(size >= v->get_size());
		return v;
	}

	static unique_ptr merge(
			const std::vector<const ext_mem_directed_vertex *> &vertices);

	template<class edge_data_type = empty_data>
	static size_t serialize(const in_mem_directed_vertex<edge_data_type> &in_v,
			char *buf, size_t size) {
		assert(size >= sizeof(ext_mem_directed_vertex));

		ext_mem_directed_vertex *ext_v = (ext_mem_directed_vertex *) buf;
		ext_v->set_id(in_v.get_id());
		ext_v->num_in_edges = in_v.get_num_in_edges();
		ext_v->num_out_edges = in_v.get_num_out_edges();
		if (in_v.has_edge_data())
			ext_v->edge_data_size = sizeof(edge_data_type);
		else
			ext_v->edge_data_size = 0;
		size_t mem_size = ext_v->get_size();
		assert(mem_size <= MAX_VERTEX_SIZE);
		assert(size >= mem_size);

		vertex_id_t *neighbors = ext_v->neighbors;
		for (size_t i = 0; i < in_v.get_num_in_edges(); i++) {
			neighbors[i] = in_v.get_in_edge(i).get_from();
		}
		for (size_t i = 0; i < in_v.get_num_out_edges(); i++) {
			neighbors[i + in_v.get_num_in_edges()] = in_v.get_out_edge(i).get_to();
		}

		// serialize edge data
		if (in_v.has_edge_data()) {
			edge_data_type *in_data
				= ext_v->get_edge_data_begin<edge_data_type>(edge_type::IN_EDGE);
			for (size_t i = 0; i < in_v.get_num_in_edges(); i++) {
				edge<edge_data_type> e = in_v.get_in_edge(i);
				in_data[i] = e.get_data();
			}
			edge_data_type *out_data
				= ext_v->get_edge_data_begin<edge_data_type>(edge_type::OUT_EDGE);
			for (size_t i = 0; i < in_v.get_num_out_edges(); i++) {
				edge<edge_data_type> e = in_v.get_out_edge(i);
				out_data[i] = e.get_data();
			}
		}
		return mem_size;
	}

	ext_mem_directed_vertex() {
		this->id = 0;
		this->num_in_edges = 0;
		this->num_out_edges = 0;
		this->edge_data_size = 0;
	}

	ext_mem_directed_vertex(vertex_id_t id, vsize_t num_in_edges,
			vsize_t num_out_edges, uint32_t edge_data_size) {
		this->id = id;
		this->num_in_edges = num_in_edges;
		this->num_out_edges = num_out_edges;
		this->edge_data_size = edge_data_size;
	}

	bool has_edge_data() const {
		return edge_data_size > 0;
	}

	size_t get_edge_data_size() const {
		return edge_data_size;
	}

	size_t get_num_edges(edge_type type) const {
		if (type == IN_EDGE)
			return get_num_in_edges();
		else if (type == OUT_EDGE)
			return get_num_out_edges();
		else
			return get_num_in_edges() + get_num_out_edges();
	}

	size_t get_num_in_edges() const {
		return num_in_edges;
	}

	size_t get_num_out_edges() const {
		return num_out_edges;
	}

	const vertex_id_t get_id() const {
		return id;
	}

	template<class edge_data_type = empty_data>
	edge_const_iterator<edge_data_type> get_out_edge_begin() const {
		const edge_data_type *data = NULL;
		if (has_edge_data())
			data = get_edge_data_begin<edge_data_type>(edge_type::OUT_EDGE);
		return edge_const_iterator<edge_data_type>(get_id(),
				neighbors + num_in_edges, data, false);
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
		const edge_data_type *data = NULL;
		if (has_edge_data())
			data = get_edge_data_begin<edge_data_type>(edge_type::IN_EDGE);
		return edge_const_iterator<edge_data_type>(get_id(),
				neighbors, data, true);
	}

	template<class edge_data_type = empty_data>
	edge_const_iterator<edge_data_type> get_in_edge_end() const {
		edge_const_iterator<edge_data_type> it
			= get_in_edge_begin<edge_data_type>();
		it += get_num_in_edges();
		return it;
	}

	size_t get_size() const {
		if (has_edge_data())
			return ((size_t) get_edge_data_addr()) - ((size_t) this)
				+ (num_in_edges + num_out_edges) * edge_data_size;
		else
			return ext_mem_directed_vertex::get_header_size()
				+ (num_in_edges + num_out_edges) * sizeof(neighbors[0]);
	}

	friend class ts_ext_mem_directed_vertex;
};

/*
 * This vertex represents an undirected vertex in the external memory.
 */
class ext_mem_undirected_vertex
{
	vertex_id_t id;
	uint32_t edge_data_size;
	vsize_t num_edges;
	vertex_id_t neighbors[0];

	bool has_edge_data() const {
		return edge_data_size > 0;
	}

	void set_id(vertex_id_t id) {
		this->id = id;
	}

	char *get_edge_data_addr() const {
		return (char *) ROUNDUP((long) (neighbors + num_edges), edge_data_size);
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

	static ext_mem_undirected_vertex *deserialize(char *buf, size_t size) {
		assert(size >= ext_mem_undirected_vertex::get_header_size());
		ext_mem_undirected_vertex *v = (ext_mem_undirected_vertex *) buf;
		assert(size >= v->get_size());
		return v;
	}

	template<class edge_data_type = empty_data>
	static size_t serialize(const in_mem_undirected_vertex<edge_data_type> &v,
			char *buf, size_t size) {
		assert(ext_mem_undirected_vertex::get_header_size() <= size);
		ext_mem_undirected_vertex *ext_v = (ext_mem_undirected_vertex *) buf;
		ext_v->set_id(v.get_id());
		ext_v->num_edges = v.get_num_edges();
		if (v.has_edge_data())
			ext_v->edge_data_size = sizeof(edge_data_type);
		else
			ext_v->edge_data_size = 0;
		size_t mem_size = ext_v->get_size();
		assert(mem_size <= MAX_VERTEX_SIZE);
		assert(size >= mem_size);

		vertex_id_t *neighbors = ext_v->neighbors;
		for (size_t i = 0; i < v.get_num_edges(); i++) {
			neighbors[i] = v.get_edge(i).get_to();
		}

		// serialize edge data
		if (v.has_edge_data()) {
			edge_data_type *data = ext_v->get_edge_data_begin<edge_data_type>();
			for (size_t i = 0; i < v.get_num_edges(); i++) {
				edge<edge_data_type> e = v.get_edge(i);
				data[i] = e.get_data();
			}
		}

		return mem_size;
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
			return ((size_t) get_edge_data_addr()) - ((size_t) this)
				+ (num_edges) * edge_data_size;
		else
			return ext_mem_undirected_vertex::get_header_size()
				+ (num_edges) * sizeof(neighbors[0]);
	}

	size_t get_num_edges(edge_type type) const {
		return num_edges;
	}

	vertex_id_t get_neighbor(edge_type type, size_t idx) const {
		return neighbors[idx];
	}

	template<class edge_data_type = empty_data>
	edge_const_iterator<edge_data_type> get_edge_begin() const {
		const edge_data_type *data = NULL;
		if (has_edge_data())
			data = get_edge_data_begin<edge_data_type>();
		return edge_const_iterator<edge_data_type>(get_id(),
				neighbors, data, false);
	}

	template<class edge_data_type = empty_data>
	edge_const_iterator<edge_data_type> get_edge_end() const {
		edge_const_iterator<edge_data_type> it
			= get_edge_begin<edge_data_type>();
		it += num_edges;
		return it;
	}

	vertex_id_t get_id() const {
		return id;
	}
};

/**
 * \brief Vertex representation when in the page cache.
 */
class page_vertex
{
public:
    /**
     * \brief Get the number of edges connecting the vertex to othe vertices.
     * \return The number of edges conning the vertex.
     * \param type The type of edges a user wishes to iterate 
                over e.g `IN_EDGE`, `OUT_EDGE`, `BOTH_EDGES`.
     */
	virtual size_t get_num_edges(edge_type type) const = 0;
    
    /**
     * \brief Get an STL-style const iterator pointing to the *first* neighbor 
            in a vertex's neighbor list.
     * \return A const iterator pointing to the *first* neighbor in 
            a vertex's neighbor list.
     * \param type The type of edges a user wishes to iterate over
            e.g `IN_EDGE`, `OUT_EDGE`, `BOTH_EDGES`.
     */
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			edge_type type) const = 0;
    
    /**
     * \brief Get an STL-style const iterator pointing to the *end* of
                    a vertex's neighbor list.
     * \return A const iterator pointing to the *end* of a vertex's neighbor list.
     * \param type The type of edges a user wishes to iterate over
                e.g `IN_EDGE`, `OUT_EDGE`, `BOTH_EDGES`.
     */
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			edge_type type) const = 0;
    /**
     * \brief Get a java-style sequential const iterator that iterates
	 *        the neighbors in the specified range.
     * \return A sequential const iterator.
     * \param type The type of edges a user wishes to iterate over e.g `IN_EDGE`, 
     *          `OUT_EDGE`, `BOTH_EDGES`.
	 * \param start The starting offset in the neighbor list iterated by
	 *              the sequential iterator.
	 * \param end The end offset in the neighbor list iterated by the sequential
	 *            iterator.
     */
	virtual page_byte_array::seq_const_iterator<vertex_id_t> get_neigh_seq_it(
			edge_type type, size_t start, size_t end) const = 0;
    
    /**
     * \brief Get the vertex unique ID.
     * \return The vertex unique ID.
     */
	virtual vertex_id_t get_id() const = 0;
    
    /**
     * \brief Read the edges of the specified type.
     * \param type The type of edges a user wishes to read
     *      e.g `IN_EDGE`, `OUT_EDGE`, `BOTH_EDGES`.
	 * \param edges The array of edges returned to a user.
	 * \param num The maximal number of edges read by a user.
     */
	virtual size_t read_edges(edge_type type, vertex_id_t edges[],
			size_t num) const {
		assert(0);
		return 0;
	}
    
    /**
     * \brief Whether the vertex is partial or complete. A user may only
	 *        request the in-edge or out-edge list of a vertex.
	 * \return true if the page vertex contains the entire vertex, else false.
     */
	virtual bool is_complete() const {
		return true;
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
	using page_vertex::get_num_edges;
    
    /**
     * \brief Get the global number of edges associated with a vertex.
     * \return The number of edges associated with a vertex.
     */
	virtual size_t get_num_edges() const = 0;
    
    /**
     * \brief Get the number of edges associated with a vertex at a specific time point.
     * \param timestamp The specific time stamp where you want the vertex metadata evaluated.
     * \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `BOTH_EDGES`.
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
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `BOTH_EDGES`.
     *  \return A const iterator pointing to the *first* element in the
     *         neighbor list of a vertex.
     */
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			int timestamp, edge_type type) const = 0;
    
    /**
     * \brief Get an STL-style const iterator pointing to the *end* of the
     *         neighbor list of a vertex at a specific time point.
     *  \param timpstamp The time stamp of interest.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `BOTH_EDGES`.
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
	const page_byte_array &array;
	bool partial;
public:
    
    /**
	 * \internal
	 * The constructor for a directed vertex in the page cache.
     *  \param arr The byte array containing the directed vertex
	 *             in the page cache.
     */
	page_directed_vertex(const page_byte_array &arr): array(arr) {
		size_t size = arr.get_size();
		// We only want to know the header of the vertex, so we don't need to
		// know what data type an edge has.
		assert(size >= ext_mem_directed_vertex::get_header_size());
		ext_mem_directed_vertex v = arr.get<ext_mem_directed_vertex>(0);
		assert(size >= v.get_size());

		id = v.get_id();
		num_in_edges = v.get_num_in_edges();
		num_out_edges = v.get_num_out_edges();
		partial = false;
	}

	/** 
	 * \internal
     * This constructor is for partial directed vertex.
     * \param id The vertex ID.
     * \param num_in_edges The number of in-edges of the vertex brought
	 *                     to the page cache.
     * \param num_out_edges The number of out-edges of the vertex brought
	 *                      to the page cache.
     * \param arr The byte array containing the partial directed vertex.
     */
	page_directed_vertex(vertex_id_t id, vsize_t num_in_edges,
			vsize_t num_out_edges, const page_byte_array &arr): array(arr) {
		this->id = id;
		this->num_in_edges = num_in_edges;
		this->num_out_edges = num_out_edges;
		assert(arr.get_size() % sizeof(vertex_id_t) == 0);
		assert(arr.get_size() / sizeof(vertex_id_t) == num_in_edges
				|| arr.get_size() / sizeof(vertex_id_t) == num_out_edges);
		this->partial = true;
	}
    
    /*
     * \brief Get the number of edges associated with a vertex.
     * \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `BOTH_EDGES`.
     * \return The number of edges associated with a vertex.
     */
	size_t get_num_edges(edge_type type) const {
		if (type == IN_EDGE)
			return num_in_edges;
		else if (type == OUT_EDGE)
			return num_out_edges;
		else if (type == BOTH_EDGES)
			return num_in_edges + num_out_edges;
		else
			assert(0);
	}
    
    /*
     * \brief Get an STL-style const iterator pointing to the *first* element in the
     *         neighbor list of a vertex.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `BOTH_EDGES`.
     *  \return A const iterator pointing to the *first* element in the
     *         neighbor list of a vertex.
     */
	page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			edge_type type) const {
		if (partial)
			return array.begin<vertex_id_t>(0);
		else if (type == IN_EDGE || type == BOTH_EDGES)
			return array.begin<vertex_id_t>(
					ext_mem_directed_vertex::get_header_size());
		else if (type == OUT_EDGE)
			return array.begin<vertex_id_t>(
					ext_mem_directed_vertex::get_header_size()
				+ num_in_edges * sizeof(vertex_id_t));
		else
			assert(0);
	}
    
    /*
     * \brief Get an STL-style const iterator pointing to the *end* of the
     *         neighbor list of a vertex.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `BOTH_EDGES`.
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
     *          `OUT_EDGE`, `BOTH_EDGES`.
	 * \param start The starting offset in the neighbor list iterated by
	 *              the sequential iterator.
	 * \param end The end offset in the neighbor list iterated by the sequential
	 *            iterator.
     */
	page_byte_array::seq_const_iterator<vertex_id_t> get_neigh_seq_it(
			edge_type type, size_t start, size_t end) const {
		if (partial)
			return array.get_seq_iterator<vertex_id_t>(
					start * sizeof(vertex_id_t), end * sizeof(vertex_id_t));

		assert(start <= end);
		switch(type) {
			case IN_EDGE:
			case BOTH_EDGES:
				assert(end <= get_num_edges(type));
				return array.get_seq_iterator<vertex_id_t>(
						ext_mem_directed_vertex::get_header_size()
						+ start * sizeof(vertex_id_t),
						ext_mem_directed_vertex::get_header_size()
						+ end * sizeof(vertex_id_t));
			case OUT_EDGE:
				assert(end <= num_out_edges);
				return array.get_seq_iterator<vertex_id_t>(
						ext_mem_directed_vertex::get_header_size()
						+ (num_in_edges + start) * sizeof(vertex_id_t),
						ext_mem_directed_vertex::get_header_size()
						+ (num_in_edges + end) * sizeof(vertex_id_t));
			default:
				assert(0);
		}
	}
    
    /**
     * \brief Get an STL-style const iterator pointing to the *first* element in the
     *         edge data list of a vertex.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `BOTH_EDGES`.
     *  \return A const iterator pointing to the *first* element in the
     *         edge data list of a vertex.
     */
	template<class edge_data_type>
	page_byte_array::const_iterator<edge_data_type> get_data_begin(
			edge_type type) const {
		assert(!partial);
		off_t edge_end
			= ext_mem_directed_vertex::get_edge_data_offset(
					num_in_edges, num_out_edges, sizeof(edge_data_type));
		if (type == IN_EDGE || type == BOTH_EDGES)
			return array.begin<edge_data_type>(edge_end);
		else if (type == OUT_EDGE)
			return array.begin<edge_data_type>(edge_end
					+ num_in_edges * sizeof(edge_data_type));
		else
			assert(0);
	}

    /**
     * \brief Get an STL-style const iterator pointing to the *end* of the
     *         edge data list of a vertex.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `BOTH_EDGES`.
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
     *          `OUT_EDGE`, `BOTH_EDGES`.
	 * \param start The starting offset in the edge list.
	 * \param end The end offset in the edge list.
	 */
	template<class edge_data_type>
	page_byte_array::seq_const_iterator<edge_data_type> get_data_seq_it(
			edge_type type, size_t start, size_t end) const {
		// TODO we currently don't support to request a partial vertex.
		assert(!partial);
		off_t edge_end
			= ext_mem_directed_vertex::get_edge_data_offset(
					num_in_edges, num_out_edges, sizeof(edge_data_type));
		switch(type) {
			case IN_EDGE:
			case BOTH_EDGES:
				assert(end <= get_num_edges(type));
				return array.get_seq_iterator<edge_data_type>(
						edge_end + start * sizeof(edge_data_type),
						edge_end + end * sizeof(edge_data_type));
			case OUT_EDGE:
				assert(end <= num_out_edges);
				return array.get_seq_iterator<edge_data_type>(
						edge_end + (num_in_edges + start) * sizeof(edge_data_type),
						edge_end + (num_in_edges + end) * sizeof(edge_data_type));
			default:
				assert(0);
		}
	}

    /**
     * \brief Get a java-style sequential const iterator that iterates
	 *        the edge data list.
     * \return A sequential const iterator.
     * \param type The type of edges a user wishes to iterate over e.g `IN_EDGE`, 
     *          `OUT_EDGE`, `BOTH_EDGES`.
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
		vsize_t num_edges = get_num_edges(type);
		assert(num_edges <= num);
		if (partial)
			array.memcpy(0, (char *) edges, sizeof(vertex_id_t) * num_edges);
		else if (type == IN_EDGE || type == BOTH_EDGES)
			array.memcpy(ext_mem_directed_vertex::get_header_size(),
					(char *) edges, sizeof(vertex_id_t) * num_edges);
		else if (type == OUT_EDGE)
			array.memcpy(ext_mem_directed_vertex::get_header_size()
				+ num_in_edges * sizeof(vertex_id_t), (char *) edges,
				sizeof(vertex_id_t) * num_edges);
		else
			assert(0);

		return num_edges;
	}
    
    /** \brief Get the id of the vertex
     *  \return The vertex id
     */
	vertex_id_t get_id() const {
		return id;
	}
    
    /** \brief Determine whether the vertex is partial or complete.
     * \return True if complete else false
     */
	bool is_complete() const {
		return !partial;
	}
};

/**
 * \brief This vertex class represents an undirected vertex in the page cache.
 */
class page_undirected_vertex: public page_vertex
{
	vertex_id_t id;
	vsize_t num_edges;
	const page_byte_array &array;
public:
	page_undirected_vertex(const page_byte_array &arr): array(arr) {
		size_t size = arr.get_size();
		assert(size >= ext_mem_undirected_vertex::get_header_size());
		// We only want to know the header of the vertex, so we don't need to
		// know what data type an edge has.
		ext_mem_undirected_vertex v = arr.get<ext_mem_undirected_vertex>(0);
		assert((unsigned) size >= v.get_size());

		id = v.get_id();
		num_edges = v.get_num_edges(BOTH_EDGES);
	}
    
    /**
     * \brief Get the number of edges associated with the vertex.
     * \param type The type of edge i.e `IN_EDGE`, `OUT_EDGE`, `BOTH_EDGES` are equivalent,
	 *             since it's an undirected vertex.
     * return The number of edges associated with the vertex.
     */
	size_t get_num_edges(edge_type type) const {
		return num_edges;
	}

    /**
     * \brief Get an STL-style const iterator pointing to the *first* element in the
     *         neighbor list of a vertex.
     * \param type The type of edge i.e `IN_EDGE`, `OUT_EDGE`, `BOTH_EDGES` are equivalent,
	 *             since it's an undirected vertex.
     *  \return A const iterator pointing to the *first* element in the
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
     * \param type The type of edge i.e `IN_EDGE`, `OUT_EDGE`, `BOTH_EDGES` are equivalent,
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
     * \param type The type of edge i.e `IN_EDGE`, `OUT_EDGE`, `BOTH_EDGES` are equivalent,
	 *             since it's an undirected vertex.
	 * \param start The starting offset in the neighbor list iterated by
	 *              the sequential iterator.
	 * \param end The end offset in the neighbor list iterated by the sequential
	 *            iterator.
     * \return A sequential const iterator for the specified range
     * in a vertex's neighbor list.
     */
	page_byte_array::seq_const_iterator<vertex_id_t> get_neigh_seq_it(
			edge_type type, size_t start, size_t end) const {
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
     * \param type The type of edge i.e `IN_EDGE`, `OUT_EDGE`, `BOTH_EDGES` are equivalent,
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

template<class edge_data_type>
class ts_in_mem_directed_vertex;

/*
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
		vsize_t num_in_edges;
		vsize_t num_out_edges;
	public:
		edge_list(vertex_id_t *neighbors, size_t num_in_edges,
				size_t num_out_edges) {
			this->neighbors = neighbors;
			this->num_in_edges = num_in_edges;
			this->num_out_edges = num_out_edges;
		}

		size_t get_num_in_edges() const {
			return num_in_edges;
		}

		size_t get_num_out_edges() const {
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
		vsize_t num_in_edges;
		vsize_t num_out_edges;
	public:
		static edge_list to_edge_list(const_edge_list &list) {
			return edge_list((vertex_id_t *) list.neighbors, list.num_in_edges,
					list.num_out_edges);
		}

		const_edge_list(const vertex_id_t *neighbors, size_t num_in_edges,
				size_t num_out_edges) {
			this->neighbors = neighbors;
			this->num_in_edges = num_in_edges;
			this->num_out_edges = num_out_edges;
		}

		size_t get_num_in_edges() const {
			return num_in_edges;
		}

		size_t get_num_out_edges() const {
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
	uint32_t edge_data_size;
	// TODO it may overflow in a graph with many timestamps.
	vsize_t num_edges;
	int num_timestamps;

	static size_t get_timestamp_list_size(int num_timestamps) {
		// The edge off list need to align with the word size.
		int num = ROUNDUP(num_timestamps, 8);
		return sizeof(short) * num;
	}

	static size_t get_timestamp_table_size(int num_timestamps) {
		return get_timestamp_list_size(num_timestamps)
			+ sizeof(edge_off) * num_timestamps;
	}

	void set_edge_data_size(uint32_t size) {
		this->edge_data_size = size;
	}

	bool has_edge_data() const {
		return edge_data_size > 0;
	}

	void set_id(vertex_id_t id) {
		this->id = id;
	}

	/*
	 * Find the specified timestamp in the timestamp table.
	 * It returns the index of the specified timestamp in the table.
	 */
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
	size_t get_header_size() const {
		return sizeof(ts_ext_mem_directed_vertex)
			+ get_timestamp_table_size(num_timestamps);
	}

	/*
	 * Get the number of in-edges in the timestamp at the location
	 * specified by `idx'.
	 */
	size_t get_num_in_edges_idx(int idx) const {
		if (idx < 0)
			return 0;
		else
			return get_edge_off_begin()[idx].out_off
				- get_edge_off_begin()[idx].in_off;
	}

	size_t get_num_out_edges_idx(int idx) const {
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

	/*
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
		size_t num_prev_edges = get_edge_off_begin()[idx].in_off;
		size_t prev_edge_list_size = num_prev_edges * (sizeof(vertex_id_t)
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
			size_t prev_edge_list_size = num_edges * (sizeof(vertex_id_t)
					+ edge_data_size);
			return const_edge_list((vertex_id_t *) (edge_list_start
						+ prev_edge_list_size), 0, 0);
		}
	}

	/*
	 * This returns the offset (in bytes) of the edge list of the specified
	 * timestamp in the external memory.
	 * If the timestamp doesn't exist, return the offset of the end
	 * of the vertex.
	 */
	off_t get_edge_list_offset(int timestamp, edge_type type) const {
		const_edge_list edges = get_const_edge_list(timestamp);
		if (type == edge_type::IN_EDGE)
			return ((char *) edges.get_in_edge_begin()) - ((char *) this);
		else if (type == edge_type::OUT_EDGE)
			return ((char *) edges.get_out_edge_begin()) - ((char *) this);
		else
			assert(0);
	}

	/*
	 * This construct a header of the vertex based on the original vertex
	 * header passed to the method.
	 * \param edge_list_off: the relative location (in bytes) of the edge lists
	 * in the vertex.
	 * \param edge_list_size: the size (in bytes) of the edge lists.
	 */
	void construct_header(const ts_ext_mem_directed_vertex &header,
			off_t edge_list_off, size_t edge_list_size);

	/*
	 * This returns the offset (in bytes) of the edge data list of the specified
	 * timestamp in the external memory.
	 * If the timestamp doesn't exist, return the offset of the end
	 * of the vertex.
	 */
	template<class edge_data_type>
	off_t get_edge_data_offset(int timestamp, edge_type type) const {
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

	ts_ext_mem_directed_vertex(vertex_id_t id, size_t num_edges,
			int num_timestamps, uint32_t edge_data_size) {
		this->id = id;
		this->num_edges = num_edges;
		this->num_timestamps = num_timestamps;
		this->edge_data_size = edge_data_size;
	}
public:
	typedef std::unique_ptr<ts_ext_mem_directed_vertex,
			delete_as_chararr<ts_ext_mem_directed_vertex> > unique_ptr;

	static size_t get_vertex_size(int num_timestamps, vsize_t num_edges,
			int edge_data_size) {
		size_t size = sizeof(ts_ext_mem_directed_vertex)
			+ get_timestamp_table_size(num_timestamps)
			+ sizeof(vertex_id_t) * num_edges;
		if (edge_data_size > 0)
			size += edge_data_size * num_edges;
		return size;
	}

	template<class edge_data_type = empty_data>
	static size_t serialize(const ts_in_mem_directed_vertex<edge_data_type> &in_v,
			char *buf, size_t size) {
		ts_ext_mem_directed_vertex *v = new (buf) ts_ext_mem_directed_vertex(
				in_v.get_id(), in_v.get_num_edges(), in_v.get_num_timestamps(),
				in_v.has_edge_data() ? sizeof(edge_data_type) : 0);
		assert(v->get_size() <= MAX_VERTEX_SIZE);
		assert(v->get_size() <= size);

		// Generate the timestamp table.
		// It's likely that many vertices don't have edges in most of timestamps.
		std::vector<int> all_timestamps;
		in_v.get_all_timestamps(all_timestamps);
		assert((size_t) v->num_timestamps == all_timestamps.size());
		size_t off = 0;
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
			size_t num_in_edges = in_v.get_num_in_edges(timestamp);
			edge_list edges = v->get_edge_list(timestamp);
			edge_const_iterator<edge_data_type> it
				= in_v.get_in_edge_begin(timestamp);
			assert(edges.get_num_in_edges() == num_in_edges);
			for (size_t j = 0; j < num_in_edges; j++) {
				edge<edge_data_type> e = *it;
				++it;
				edges.get_in_edge_begin()[j] = e.get_from();
				if (v->has_edge_data())
					edges.get_in_edge_data<edge_data_type>()[j] = e.get_data();
			}

			size_t num_out_edges = in_v.get_num_out_edges(timestamp);
			it = in_v.get_out_edge_begin(timestamp);
			assert(edges.get_num_out_edges() == num_out_edges);
			for (size_t j = 0; j < num_out_edges; j++) {
				edge<edge_data_type> e = *it;
				++it;
				edges.get_out_edge_begin()[j] = e.get_to();
				if (v->has_edge_data())
					edges.get_out_edge_data<edge_data_type>()[j] = e.get_data();
			}
		}
		return v->get_size();
	}

	static unique_ptr merge(
			const std::vector<const ext_mem_directed_vertex *> &vertices);

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
		return get_vertex_size(num_timestamps, num_edges, edge_data_size);
	}

	size_t get_num_edges() const {
		return num_edges;
	}

	int get_num_timestamps() const {
		return num_timestamps;
	}

	uint32_t get_edge_data_size() const {
		return edge_data_size;
	}

	size_t get_num_in_edges(int timestamp) const {
		int idx = find_timestamp(timestamp);
		return get_num_in_edges_idx(idx);
	}

	size_t get_num_out_edges(int timestamp) const {
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

	size_t get_num_in_edges() const {
		size_t num = 0;
		// TODO I should use a simpler way to compute the result.
		for (int i = 0; i < num_timestamps; i++)
			num += get_num_in_edges(get_timestamps_begin()[i]);
		return num;
	}

	size_t get_num_out_edges() const {
		size_t num = 0;
		// TODO I should use a simpler way to compute the result.
		for (int i = 0; i < num_timestamps; i++)
			num += get_num_out_edges(get_timestamps_begin()[i]);
		return num;
	}

	void print() const {
		printf("v%ld has edge data: %d, # timestamps: %d, # edges: %ld\n",
				(unsigned long) get_id(), 0, num_timestamps, get_num_edges());
		for (int i = 0; i < num_timestamps; i++) {
			int timestamp = get_timestamps_begin()[i];
			// We need to skip the timestamps without edges.
			if (get_num_in_edges(timestamp) + get_num_out_edges(timestamp) == 0)
				continue;

			printf("timestamp %d\n", timestamp);
			size_t num_in_edges = get_num_in_edges(timestamp);
			printf("in-edges (%ld): ", num_in_edges);
			const vertex_id_t *in_edge_list = get_edge_list_begin()
				+ get_edge_off_begin()[i].in_off;
			for (size_t j = 0; j < num_in_edges; j++) {
				printf("%ld, ", (unsigned long) in_edge_list[j]);
			}
			printf("\n");
			size_t num_out_edges = get_num_out_edges(timestamp);
			printf("out-edges (%ld): ", num_out_edges);
			const vertex_id_t *out_edge_list = get_edge_list_begin()
				+ get_edge_off_begin()[i].out_off;
			for (size_t j = 0; j < num_out_edges; j++) {
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
 * \brief This is a data structure to interpret a time-series directed vertex
 * in the page cache. <br>
 *
 * The memory layout of a time-series directed vertex is as follows:
 *      - id
 *      - the total number of edges
 *      - the number of timestamps
 *      - timestamp table that contains the offsets to the edges for the timestamp
 *        within the edge list
 *      - edge offsets of existing timestamp (in-edge offset and out-edge offset)
 *      - edges
 *      - edge data <br>
 *
 * The size of this object is defined at the runtime, so it can only be
 * allocated from the heap.
 */
class TS_page_directed_vertex: public TS_page_vertex
{
	// The entire size of the vertex in the external memory.
	size_t entire_vertex_size;
	// The location of the vertex in the file.
	off_t vertex_begin_off;
	const page_byte_array *array;
	ts_ext_mem_directed_vertex ext_v;
    
    /**
	 * \internal
	 * The constructor for a time-series directed vertex in the page cache.
     *  \param arr The byte array containing the time-series directed vertex
	 *             in the page cache.
     */
	TS_page_directed_vertex(const page_byte_array &arr): array(&arr) {
		unsigned size = arr.get_size();
		assert((unsigned) size >= sizeof(ts_ext_mem_directed_vertex));
		ts_ext_mem_directed_vertex v = arr.get<ts_ext_mem_directed_vertex>(0);
		entire_vertex_size = v.get_size();
		vertex_begin_off = arr.get_offset();
		// We have to make sure the array size must be larger than the header
		// of the ts-vertex.
		size_t header_size = v.get_header_size();
		assert(arr.get_size() >= header_size);
		assert(header_size <= PAGE_SIZE);
		arr.memcpy(0, (char *) &ext_v, header_size);
		// If the vertex isn't complete, we probably only has the header
		// of the vertex. Don't give a uaser a chance to access it.
		if (!is_complete())
			array = NULL;
	}

	/**
	 * \internal
	 * \brief This is to create a vertex based on the header of the vertex
	 * and the edge list in the array. <br>
	 * The constructed vertex has only part of the vertex.
	 */
	TS_page_directed_vertex(const TS_page_directed_vertex *header,
			const page_byte_array &arr): array(&arr) {
		this->entire_vertex_size = header->entire_vertex_size;
		this->vertex_begin_off = header->vertex_begin_off;
		off_t rel_off = arr.get_offset() - vertex_begin_off;
		assert((size_t) rel_off >= header->ext_v.get_header_size());
		ext_v.construct_header(header->ext_v, rel_off, arr.get_size());
	}

	/**
	 * \brief **This object is not allowed to be copied** hence
     *      the copy constructor is disabled.
    */
	TS_page_directed_vertex(const TS_page_directed_vertex &);
    
    /** \brief **This object is not allowed to be copied** hence
     the copy constructor is disabled.
     */
	TS_page_directed_vertex &operator=(const TS_page_directed_vertex &);
    
    /**
     * \brief Get an STL-style const iterator pointing to the *end* of the
     *         neighbor list of a vertex.
     */
	page_byte_array::const_iterator<vertex_id_t> get_neigh_end() const {
		assert(array);
		return array->begin<vertex_id_t>(array->get_size());
	}
    
    /**
     * \brief Get an STL-style const iterator pointing to the *end*
     * of the edge data list of a vertex.
     */
	template<class edge_data_type>
	page_byte_array::const_iterator<edge_data_type> get_edge_data_end() const {
		assert(array);
		return array->begin<edge_data_type>(array->get_size());
	}
public:
	/** \brief The header size of the vertex.
     * \param num_timestamps The number of timestamps in the timestamp table of
	 * the vertex.
     */
	static size_t get_size(int num_timestamps) {
		return sizeof(TS_page_directed_vertex)
			+ ts_ext_mem_directed_vertex::get_timestamp_table_size(num_timestamps);
	}

	/** 
     * \brief We create the vertex object in the given 
     *           buffer in lieu of calling the constructor.
     *  \param arr A byte array in the page cache that contains the data of
	 *             the time-series vertex.
	 *  \param buf A buffer where we construct the object.
	 *  \param size The buffer size.
     */
	static TS_page_directed_vertex *create(const page_byte_array &arr,
			char *buf, size_t size) {
		return new (buf) TS_page_directed_vertex(arr);
	}
    
    /**
     * \brief We create the vertex object in the given
     *           buffer in lieu of calling the constructor.
     *  \param header A pointer to the time series page directed vertex header.
     *  \param arr A byte array in the page cache that contains the data of
	 *             the time-series vertex.
	 *  \param buf A buffer where we construct the object.
	 *  \param size The buffer size.
     */
	static TS_page_directed_vertex *create(const TS_page_directed_vertex *header,
			const page_byte_array &arr, char *buf, size_t size) {
		return new (buf) TS_page_directed_vertex(header, arr);
	}

	/**
	 * \brief This returns the relative offset of the edge lists specified
	 *          in the timestamp range in the vertex.
     *  \param range The timestamp range in the vertex.
     *  \return The relative offset of the edge lists specified
	 *          in the timestamp range in the vertex.
	 */
	virtual offset_pair get_edge_list_offset(const timestamp_pair &range) const;
    
    /**
     * \brief Check if the vertex is complete or partial.
     * \return True if complete, else false.
     */
	virtual bool is_complete() const {
		if (array == NULL)
			return false;
		else
			return array->get_size() >= entire_vertex_size;
	}
    
    /**
     * \brief Get the number of edges of a certain `edge_type` associated with
     *          the vertex.
     *  \param type The type of edge i.e `OUT_EDGE`, `IN_EDGE`, `BOTH_EDGES`.
     *  \return The number of edges.
     */
	virtual size_t get_num_edges(edge_type type) const {
		assert(0);
	}
    
    /**
     * \brief Get an STL-style const iterator pointing to the *first* element in the
     *         neighbor list of a vertex.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `BOTH_EDGES`.
     *  \return A const iterator pointing to the *first* element in the
     *         neighbor list of a vertex.
     */
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			edge_type type) const {
		assert(0);
	}
    
    /**
     * \brief Get an STL-style const iterator pointing to the *end* of the
     *         neighbor list of a vertex.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`, `BOTH_EDGES`.
     *  \return A const iterator pointing to the *end* of the
     *         neighbor list of a vertex.
     */
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_end(
			edge_type type) const {
		assert(0);
	}
    
    /**
     *  \brief Get a java-style sequential const iterator pointing to the *first*
            neighbor in a vertex's neighbor list.
     *
     *  \param type The type of edges a user wishes to iterate over e.g `IN_EDGE`,
     *          `OUT_EDGE`, `BOTH_EDGES`.
	 * \param start The starting offset in the neighbor list iterated by
	 *              the sequential iterator.
	 * \param end The end offset in the neighbor list iterated by the sequential
	 *            iterator.
     *  \return A sequential const iterator pointing to the *first* neighbor in a
     vertex's neighbor list.
     */
	virtual page_byte_array::seq_const_iterator<vertex_id_t> get_neigh_seq_it(
			edge_type type, size_t start, size_t end) const {
		assert(0);
	}
    
    /**
     * \brief Get the vertex unique ID.
     * \return The vertex unique ID.
     */
	virtual vertex_id_t get_id() const {
		return ext_v.get_id();
	}
    
    /**
     * \brief Get the **global** number of edges associated with the vertex.
     *  \return The number of edges associated with the vertex.
     */
	virtual size_t get_num_edges() const {
		return ext_v.get_num_edges();
	}
    
    /**
     * \brief Get the number of time stamps where this vertex participates in the graph.
     * \return The number of time stamps where this vertex participates in the graph.
     */
	virtual int get_num_timestamps() const {
		return ext_v.get_num_timestamps();
	}
    
    
    /**
     * \brief Get the number of edges associated of a single `edge_type` 
     *          with the vertex at a given timestamp.
     *  \param timestamp The timestamp of interest.
     *  \param The edge type i.e `OUT_EDGE`, `IN_EDGE`, `BOTH_EDGES`.
     */
	virtual size_t get_num_edges(int timestamp, edge_type type) const {
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
    
    /**
     * \brief Get an STL-style const iterator pointing to the *first* element in
     *		the neighbor list of a vertex at a given timestamp.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`,
     *			`BOTH_EDGES`.
     *  \param timestamp The timestamp of interest.
     *  \return A const iterator pointing to the *first* element in the
     *         neighbor list of a vertex.
     */
	virtual page_byte_array::const_iterator<vertex_id_t> get_neigh_begin(
			int timestamp, edge_type type) const {
		off_t offset = 0;
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
    
    /**
     * \brief Get an STL-style const iterator pointing to the *end* of
     *		the neighbor list of a vertex at a given timestamp.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`,
     *			`BOTH_EDGES`.
     *  \param timestamp The timestamp of interest.
     *  \return A const iterator pointing to the *end* of the
     *         neighbor list of a vertex.
     */
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
    
    /**
     * \brief Get an STL-style const iterator pointing to the *first* element in
     *		the edge data list of a vertex at a given timestamp.
     *  \param timestamp The timestamp of interest.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`,
     *			`BOTH_EDGES`.
     *  \return A const iterator pointing to the *first* element in the
     *         edge data of a vertex.
     */
	template<class edge_data_type>
	page_byte_array::const_iterator<edge_data_type> get_edge_data_begin(
			int timestamp, edge_type type) const {
		off_t offset = 0;
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
    
    /**
     * \brief Get an STL-style const iterator pointing to the *end* of
     *		the edge data of a vertex at a given timestamp.
     *  \param timestamp The timestamp of interest.
     *  \param type The type of edges a user wishes to evaluate e.g `IN_EDGE`,
     *			`BOTH_EDGES`.
     *  \return A const iterator pointing to the *end* of the
     *         edge data of a vertex.
     */
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
    
    /**
	 * \internal
     * \brief A visual representation of the vertex printed to stdout.
     *  \return The vertex's string representation.
     */
	virtual void print() const {
		printf("v%ld has edge data: %d, # timestamps: %d, # edges: %ld\n",
				(unsigned long) get_id(), 0, get_num_timestamps(), get_num_edges());
		for (int i = 0; i < get_num_timestamps(); i++) {
			int timestamp = ext_v.get_timestamps_begin()[i];
			// We need to skip the timestamps without edges.
			if (get_num_edges(timestamp, edge_type::IN_EDGE)
					+ get_num_edges(timestamp, edge_type::OUT_EDGE) == 0)
				continue;

			printf("timestamp %d\n", timestamp);
			size_t num_in_edges = get_num_edges(timestamp, edge_type::IN_EDGE);
			printf("in-edges (%ld): ", num_in_edges);
			page_byte_array::const_iterator<vertex_id_t> end_it
				= get_neigh_end(timestamp, edge_type::IN_EDGE);
			for (page_byte_array::const_iterator<vertex_id_t> it
					= get_neigh_begin(timestamp, edge_type::IN_EDGE);
					it != end_it; ++it) {
				printf("%ld, ", (unsigned long) *it);
			}
			printf("\n");
			size_t num_out_edges = get_num_edges(timestamp, edge_type::OUT_EDGE);
			printf("out-edges (%ld): ", num_out_edges);
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
	virtual size_t get_serialize_size() const = 0;
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

	in_mem_directed_vertex(const ext_mem_directed_vertex *v) {
		id = v->get_id();
		has_data = v->has_edge_data();

		edge_const_iterator<edge_data_type> in_it
			= v->get_in_edge_begin<edge_data_type>();
		edge_const_iterator<edge_data_type> in_end
			= v->get_in_edge_end<edge_data_type>();
		for (; in_it != in_end; ++in_it)
			this->add_in_edge(*in_it);

		edge_const_iterator<edge_data_type> out_it
			= v->get_out_edge_begin<edge_data_type>();
		edge_const_iterator<edge_data_type> out_end
			= v->get_out_edge_end<edge_data_type>();
		for (; out_it != out_end; ++out_it)
			this->add_out_edge(*out_it);
	}

	vertex_id_t get_id() const {
		return id;
	}

	bool has_edge_data() const {
		return has_data;
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
				assert(0);
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

	size_t get_serialize_size() const {
		ext_mem_directed_vertex v(0, get_num_in_edges(), get_num_out_edges(),
				has_data ? sizeof(edge_data_type) : 0);
		return v.get_size();
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
	in_mem_undirected_vertex(vertex_id_t id, bool has_data) {
		this->id = id;
		this->has_data = has_data;
	}

	in_mem_undirected_vertex(const page_undirected_vertex &vertex,
			bool has_data) {
		this->id = vertex.get_id();
		this->has_data = has_data;
		edges.resize(vertex.get_num_edges(edge_type::BOTH_EDGES));
		vertex.read_edges(edge_type::BOTH_EDGES, edges.data(),
				vertex.get_num_edges(edge_type::BOTH_EDGES));
		assert(!has_data);
	}

	vertex_id_t get_id() const {
		return id;
	}

	bool has_edge_data() const {
		return has_data;
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

	size_t get_num_edges(edge_type type = edge_type::BOTH_EDGES) const {
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

	size_t get_serialize_size() const {
		ext_mem_undirected_vertex v(0, edges.size(),
				has_data ? sizeof(edge_data_type) : 0);
		return v.get_size();
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

	size_t get_num_in_edges(int timestamp) const {
		typename std::map<int, ts_edge_pair>::const_iterator it
			= ts_edges.find(timestamp);
		if (it == ts_edges.end())
			return 0;
		else
			return it->second.in_edges.size();
	}

	size_t get_num_out_edges(int timestamp) const {
		typename std::map<int, ts_edge_pair>::const_iterator it
			= ts_edges.find(timestamp);
		if (it == ts_edges.end())
			return 0;
		else
			return it->second.out_edges.size();
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
				assert(0);
		}
	}

	size_t get_num_in_edges() const {
		size_t num_edges = 0;
		for (typename std::map<int, ts_edge_pair>::const_iterator it
				= ts_edges.begin(); it != ts_edges.end(); it++)
			num_edges += it->second.in_edges.size();
		return num_edges;
	}

	size_t get_num_out_edges() const {
		size_t num_edges = 0;
		for (typename std::map<int, ts_edge_pair>::const_iterator it
				= ts_edges.begin(); it != ts_edges.end(); it++)
			num_edges += it->second.out_edges.size();
		return num_edges;
	}

	size_t get_num_edges() const {
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

	size_t get_serialize_size() const {
		ts_ext_mem_directed_vertex v(get_id(), get_num_edges(),
				get_num_timestamps(),
				has_edge_data() ? sizeof(edge_data_type) : 0);
		return v.get_size();
	}

	void print() const {
		printf("v%ld has edge data: %d, # timestamps: %d, # edges: %ld\n",
				(unsigned long) get_id(), has_edge_data(),
				get_num_timestamps(), get_num_edges());
		for (typename std::map<int, ts_edge_pair>::const_iterator it
				= ts_edges.begin(); it != ts_edges.end(); it++) {
			printf("timestamp %d\n", it->first);
			printf("in-edges (%ld): ", get_num_in_edges(it->first));
			for (size_t i = 0; i < it->second.in_edges.size(); i++) {
				printf("%ld, ", (unsigned long) it->second.in_edges[i]);
			}
			printf("\n");
			printf("out-edges (%ld): ", get_num_out_edges(it->first));
			for (size_t i = 0; i < it->second.out_edges.size(); i++) {
				printf("%ld, ",  (unsigned long)it->second.out_edges[i]);
			}
			printf("\n");
		}
	}
};

#endif
