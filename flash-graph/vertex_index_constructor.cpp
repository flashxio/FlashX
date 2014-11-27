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

#include "exception.h"

#include "vertex_index_constructor.h"

namespace fg
{

/*
 * Compressed in-memory vertex index for an undirected graph.
 */
class cundirected_vertex_index_construct: public vertex_index_construct
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
	cundirected_vertex_index_construct(size_t edge_data_size) {
		this->edge_data_size = edge_data_size;
		last_entries.push_back(vertex_offset(sizeof(graph_header)));
	}

	virtual void add_vertex(const in_mem_vertex &v);

	virtual void dump(const std::string &file, const graph_header &header,
			bool compressed);
	virtual vertex_index::ptr dump(const graph_header &header, bool compressed);
};

class undirected_vertex_index_construct: public vertex_index_construct
{
	std::vector<vertex_offset> vertices;
public:
	undirected_vertex_index_construct() {
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
			undirected_vertex_index::dump(file, header, vertices);
	}

	virtual vertex_index::ptr dump(const graph_header &header, bool compressed) {
		if (compressed) {
			vertex_index::ptr index = undirected_vertex_index::create(header,
					vertices);
			undirected_vertex_index::ptr def_index
				= std::static_pointer_cast<undirected_vertex_index, vertex_index>(
						index);
			cundirected_vertex_index::ptr cindex
				= cundirected_vertex_index::construct(*def_index);
			return std::static_pointer_cast<vertex_index,
				   cundirected_vertex_index>(cindex);
		}
		else
			return undirected_vertex_index::create(header, vertices);
	}
};

/*
 * Compressed in-memory vertex index for a directed graph.
 */
class cdirected_vertex_index_construct: public vertex_index_construct
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
	cdirected_vertex_index_construct(size_t edge_data_size) {
		this->edge_data_size = edge_data_size;
		last_entries.push_back(directed_vertex_entry(sizeof(graph_header), 0));
	}

	virtual void add_vertex(const in_mem_vertex &v);

	virtual void dump(const std::string &file, const graph_header &header,
			bool compressed);
	virtual vertex_index::ptr dump(const graph_header &header, bool compressed);
};

class directed_vertex_index_construct: public vertex_index_construct
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
	directed_vertex_index_construct() {
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

void cundirected_vertex_index_construct::finalize()
{
	if (last_entries.empty())
		return;

	compressed_undirected_vertex_entry centry(last_entries.data(),
			edge_data_size, last_entries.size());
	add_large_vertices(centries.size() * ENTRY_SIZE, centry, last_entries);
	centries.push_back(centry);
	last_entries.clear();
}

void cundirected_vertex_index_construct::add_large_vertices(vertex_id_t start_vertex_id,
		const compressed_undirected_vertex_entry entry,
		const std::vector<vertex_offset> &entries)
{
	int num_entries = entries.size() - 1;
	for (int i = 0;i < num_entries; i++) {
		if (entry.is_large_vertex(i)) {
			size_t size = entries[i + 1].get_off() - entries[i].get_off();
			large_vertices.push_back(large_vertex_t(start_vertex_id + i,
					ext_mem_undirected_vertex::vsize2num_edges(
						size, edge_data_size)));
		}
	}
}

void cundirected_vertex_index_construct::add_vertex(const in_mem_vertex &v)
{
	assert(get_num_vertices() == v.get_id());
	assert(last_entries.size() <= ENTRY_SIZE + 1);
	assert(!last_entries.empty());
	vertex_offset last_entry = last_entries.back();
	if (last_entries.size() == ENTRY_SIZE + 1) {
		compressed_undirected_vertex_entry centry(last_entries.data(),
				edge_data_size, last_entries.size());
		add_large_vertices(centries.size() * ENTRY_SIZE, centry, last_entries);
		centries.push_back(centry);
		last_entries.clear();
		last_entries.push_back(last_entry);
	}
	vertex_offset entry;
	entry.init(last_entry, v);
	last_entries.push_back(entry);
}

void cundirected_vertex_index_construct::dump(const std::string &file,
		const graph_header &header, bool compressed)
{
	throw safs::unsupported_exception();
}

vertex_index::ptr cundirected_vertex_index_construct::dump(const graph_header &header,
		bool compressed)
{
	finalize();
	return std::static_pointer_cast<vertex_index>(
			cundirected_vertex_index::construct(centries, large_vertices, header));
}

void cdirected_vertex_index_construct::add_large_vertices(
		vertex_id_t start_vertex_id, const compressed_directed_vertex_entry entry,
		const std::vector<directed_vertex_entry> &entries)
{
	int num_entries = entries.size() - 1;
	for (int i = 0;i < num_entries; i++) {
		if (entry.is_large_in_vertex(i)) {
			size_t size = entries[i + 1].get_in_off() - entries[i].get_in_off();
			large_in_vertices.push_back(large_vertex_t(start_vertex_id + i,
					ext_mem_undirected_vertex::vsize2num_edges(
						size, edge_data_size)));
		}
		if (entry.is_large_out_vertex(i)) {
			size_t size = entries[i + 1].get_out_off() - entries[i].get_out_off();
			large_out_vertices.push_back(large_vertex_t(start_vertex_id + i,
					ext_mem_undirected_vertex::vsize2num_edges(
						size, edge_data_size)));
		}
	}
}

void cdirected_vertex_index_construct::add_vertex(const in_mem_vertex &v)
{
	assert(get_num_vertices() == v.get_id());
	assert(last_entries.size() <= ENTRY_SIZE + 1);
	assert(!last_entries.empty());
	directed_vertex_entry last_entry = last_entries.back();
	if (last_entries.size() == ENTRY_SIZE + 1) {
		compressed_directed_vertex_entry centry(last_entries.data(),
				edge_data_size, last_entries.size());
		add_large_vertices(centries.size() * ENTRY_SIZE, centry, last_entries);
		centries.push_back(centry);
		last_entries.clear();
		last_entries.push_back(last_entry);
	}
	directed_vertex_entry entry;
	entry.init(last_entry, v);
	last_entries.push_back(entry);
}

void cdirected_vertex_index_construct::finalize()
{
	if (last_entries.empty())
		return;

	compressed_directed_vertex_entry centry(last_entries.data(),
			edge_data_size, last_entries.size());
	add_large_vertices(centries.size() * ENTRY_SIZE, centry, last_entries);
	centries.push_back(centry);
	directed_vertex_entry last_entry = last_entries.back();
	last_entries.clear();

	size_t in_part_size = last_entry.get_in_off();
	// If all out-edge lists have been moved to behind in-edge lists,
	// ignore it.
	assert((size_t) centries.front().get_start_out_off() < in_part_size);
	for (size_t i = 0; i < centries.size(); i++) {
		off_t start_in_off = centries[i].get_start_in_off();
		off_t start_out_off = centries[i].get_start_out_off();
		centries[i].reset_start_offs(start_in_off, start_out_off + in_part_size);
	}
}

void cdirected_vertex_index_construct::dump(const std::string &file,
		const graph_header &header, bool compressed)
{
	throw safs::unsupported_exception();
}

vertex_index::ptr cdirected_vertex_index_construct::dump(const graph_header &header,
		bool compressed)
{
	finalize();
	return std::static_pointer_cast<vertex_index>(
			cdirected_vertex_index::construct(centries, large_in_vertices,
				large_out_vertices, header));
}

vertex_index_construct::ptr vertex_index_construct::create(bool directed)
{
	if (directed)
		return vertex_index_construct::ptr(new undirected_vertex_index_construct());
	else
		return vertex_index_construct::ptr(new directed_vertex_index_construct());
}

vertex_index_construct::ptr vertex_index_construct::create_compressed(bool directed,
		size_t edge_data_size)
{
	if (directed)
		return vertex_index_construct::ptr(new cundirected_vertex_index_construct(
					edge_data_size));
	else
		return vertex_index_construct::ptr(new cdirected_vertex_index_construct(
					edge_data_size));
}

}
