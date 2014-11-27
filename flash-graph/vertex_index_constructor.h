#ifndef __VERTEX_INDEX_CONSTRUCTOR_H__
#define __VERTEX_INDEX_CONSTRUCTOR_H__

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
#include "vertex_index.h"

namespace fg
{

/*
 * These are the in-mem counterparts of vertex index above.
 * These in-memory data structures are used to construct vertex indices.
 */

class in_mem_vertex_index
{
public:
	virtual ~in_mem_vertex_index() {
	}

	virtual void add_vertex(const in_mem_vertex &) = 0;
	virtual void dump(const std::string &file, const graph_header &header,
			bool compressed) = 0;
	virtual vertex_index::ptr dump(const graph_header &header, bool compressed) = 0;
};

/*
 * Compressed in-memory vertex index for an undirected graph.
 */
class cdefault_in_mem_vertex_index: public in_mem_vertex_index
{
	static const size_t ENTRY_SIZE = compressed_vertex_entry::ENTRY_SIZE;
	size_t edge_data_size;
	std::vector<compressed_undirected_vertex_entry> centries;
	// We temporarily keep the vertex entries in an uncompressed form.
	// Once the number of uncompressed vertex entries is larger than
	// ENTRY_SIZE, we convert them to the compressed form and keep them
	// in `centries'.
	std::vector<vertex_offset> last_entries;
	std::vector<large_vertex_t> large_vertices;

	void finalize();
	void add_large_vertices(vertex_id_t start_vertex_id,
			const compressed_undirected_vertex_entry entry,
			const std::vector<vertex_offset> &vertices);
	size_t get_num_vertices() const {
		return centries.size() * ENTRY_SIZE + last_entries.size() - 1;
	}
public:
	cdefault_in_mem_vertex_index(size_t edge_data_size) {
		this->edge_data_size = edge_data_size;
		last_entries.push_back(vertex_offset(sizeof(graph_header)));
	}

	virtual void add_vertex(const in_mem_vertex &v);

	virtual void dump(const std::string &file, const graph_header &header,
			bool compressed);
	virtual vertex_index::ptr dump(const graph_header &header, bool compressed);

	friend class cundirected_vertex_index;
};

class default_in_mem_vertex_index: public in_mem_vertex_index
{
	std::vector<vertex_offset> vertices;
public:
	default_in_mem_vertex_index() {
		vertices.push_back(vertex_offset(sizeof(graph_header)));
	}

	virtual void add_vertex(const in_mem_vertex &v) {
		assert(v.get_id() + 1 == vertices.size());
		vertex_offset off;
		off.init(vertices.back(), v);
		vertices.push_back(off);
	}

	virtual void dump(const std::string &file, const graph_header &header,
			bool compressed) {
		if (compressed)
			dump(header, compressed)->dump(file);
		else
			default_vertex_index::dump(file, header, vertices);
	}

	virtual vertex_index::ptr dump(const graph_header &header, bool compressed) {
		if (compressed) {
			vertex_index::ptr index = default_vertex_index::create(header,
					vertices);
			default_vertex_index::ptr def_index
				= std::static_pointer_cast<default_vertex_index, vertex_index>(
						index);
			cundirected_vertex_index::ptr cindex
				= cundirected_vertex_index::construct(*def_index);
			return std::static_pointer_cast<vertex_index,
				   cundirected_vertex_index>(cindex);
		}
		else
			return default_vertex_index::create(header, vertices);
	}
};

typedef default_in_mem_vertex_index undirected_in_mem_vertex_index;
typedef default_in_mem_vertex_index ts_directed_in_mem_vertex_index;

/*
 * Compressed in-memory vertex index for a directed graph.
 */
class cdirected_in_mem_vertex_index: public in_mem_vertex_index
{
	static const size_t ENTRY_SIZE = compressed_vertex_entry::ENTRY_SIZE;
	size_t edge_data_size;
	std::vector<compressed_directed_vertex_entry> centries;
	// We temporarily keep the vertex entries in an uncompressed form.
	// Once the number of uncompressed vertex entries is larger than
	// ENTRY_SIZE, we convert them to the compressed form and keep them
	// in `centries'.
	std::vector<directed_vertex_entry> last_entries;
	std::vector<large_vertex_t> large_in_vertices;
	std::vector<large_vertex_t> large_out_vertices;

	void finalize();
	void add_large_vertices(vertex_id_t start_vertex_id,
			const compressed_directed_vertex_entry entry,
			const std::vector<directed_vertex_entry> &vertices);
	size_t get_num_vertices() const {
		return centries.size() * ENTRY_SIZE + last_entries.size() - 1;
	}
public:
	cdirected_in_mem_vertex_index(size_t edge_data_size) {
		this->edge_data_size = edge_data_size;
		last_entries.push_back(directed_vertex_entry(sizeof(graph_header), 0));
	}

	virtual void add_vertex(const in_mem_vertex &v);

	virtual void dump(const std::string &file, const graph_header &header,
			bool compressed);
	virtual vertex_index::ptr dump(const graph_header &header, bool compressed);

	friend class cdirected_vertex_index;
};

class directed_in_mem_vertex_index: public in_mem_vertex_index
{
	std::vector<directed_vertex_entry> vertices;

	void finalize() {
		size_t in_part_size = vertices.back().get_in_off();
		// If all out-edge lists have been moved to behind in-edge lists,
		// ignore it.
		if ((size_t) vertices.front().get_out_off() >= in_part_size)
			return;
		for (size_t i = 0; i < vertices.size(); i++)
			vertices[i] = directed_vertex_entry(vertices[i].get_in_off(),
					vertices[i].get_out_off() + in_part_size);
	}
public:
	directed_in_mem_vertex_index() {
		vertices.push_back(directed_vertex_entry(sizeof(graph_header), 0));
	}

	virtual void add_vertex(const in_mem_vertex &v) {
		assert(v.get_id() + 1 == vertices.size());
		directed_vertex_entry entry;
		entry.init(vertices.back(), v);
		vertices.push_back(entry);
	}

	virtual void dump(const std::string &file, const graph_header &header,
			bool compressed) {
		finalize();
		if (compressed)
			dump(header, compressed)->dump(file);
		else
			directed_vertex_index::dump(file, header, vertices);
	}

	virtual vertex_index::ptr dump(const graph_header &header, bool compressed) {
		finalize();
		if (compressed) {
			vertex_index::ptr index = directed_vertex_index::create(header,
					vertices);
			directed_vertex_index::ptr dir_index
				= std::static_pointer_cast<directed_vertex_index, vertex_index>(
						index);
			cdirected_vertex_index::ptr cindex
				= cdirected_vertex_index::construct(*dir_index);
			return std::static_pointer_cast<vertex_index,
				   cdirected_vertex_index>(cindex);
		}
		else
			return directed_vertex_index::create(header, vertices);
	}
};

}

#endif
