#ifndef __VERTEX_INDEX_H__
#define __VERTEX_INDEX_H__

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

#include <string.h>

#include <string>

#include "native_file.h"

#include "vertex.h"
#include "graph_file_header.h"

class vertex_index
{
protected:
	graph_header header;
	// The total size of the graph in the form of adjacency list
	// in the external memory.
	size_t tot_size;
	// The size of the entire index
	size_t index_size;

	vertex_index() {
		tot_size = 0;
		index_size = 0;
	}

	size_t get_serialize_size() const {
		return index_size;
	}

public:
	static vertex_index *load(const std::string &index_file) {
		native_file local_f(index_file);
		ssize_t size = local_f.get_size();
		assert(size > 0);
		assert((size_t) size >= sizeof(vertex_index));
		char *buf = (char *) malloc(size);
		assert(buf);
		FILE *fd = fopen(index_file.c_str(), "r");
		size_t ret = fread(buf, size, 1, fd);
		assert(ret == 1);
		fclose(fd);

		vertex_index *idx = (vertex_index *) buf;
		assert((size_t) size >= idx->index_size);
		idx->header.verify();

		return idx;
	}

	static void destroy(vertex_index *index) {
		free(index);
	}

	void dump(const std::string &file) {
		FILE *f = fopen(file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			assert(0);
		}

		ssize_t ret = fwrite(this, get_serialize_size(), 1, f);
		assert(ret);

		fclose(f);
	}

	const graph_header &get_graph_header() const {
		return header;
	}

	size_t get_num_vertices() const {
		return header.get_num_vertices();
	}

	vertex_id_t get_max_id() const {
		return get_num_vertices() - 1;
	}

	size_t get_graph_size() const {
		return tot_size;
	}
};

/**
 * This vertex index maps a vertex id to the location of the vertex in a file.
 */
template <class vertex_entry_type>
class vertex_index_temp: public vertex_index
{
	vertex_entry_type vertices[0];

protected:
	vertex_index_temp() {
	}

	vertex_index_temp(size_t num) {
		memset(vertices, 0, sizeof(vertices[0]) * num);
	}
public:
	template<class vertex_type>
	static vertex_index_temp<vertex_entry_type> *create(const graph_header &header,
			const std::vector<vertex_type> &vertices) {
		void *addr = malloc(sizeof(vertex_index_temp<vertex_entry_type>)
				+ sizeof(vertex_entry_type) * vertices.size());
		assert(addr);
		assert(header.get_num_vertices() == vertices.size());
		vertex_index_temp<vertex_entry_type> *index
			= new (addr) vertex_index_temp<vertex_entry_type>(vertices.size());
		index->header = header;
		assert(sizeof(index->header) == PAGE_SIZE);
		// All data of adjacency lists are stored after the header.
		size_t tot_size = sizeof(header);
		for (size_t i = 0; i < vertices.size(); i++) {
			assert(i == vertices[i].get_id());
			index->vertices[i].init(tot_size, vertices[i]);
			tot_size += vertices[i].get_serialize_size();
		}
		index->tot_size = tot_size;
		index->index_size = sizeof(vertex_index_temp<vertex_entry_type>)
			+ index->get_num_vertices() * sizeof(index->vertices[0]);

		return index;
	}

	const vertex_entry_type &get_vertex(vertex_id_t id) const {
		assert(id < get_num_vertices());
		return vertices[id];
	}

	off_t get_vertex_off(vertex_id_t id) const {
		assert(id < get_num_vertices());
		return vertices[id].get_off();
	}

	size_t get_vertex_size(vertex_id_t id) const {
		assert(id < get_num_vertices());
		if (id < get_num_vertices() - 1)
			return vertices[id + 1].get_off() - vertices[id].get_off();
		else
			return tot_size - vertices[id].get_off();
	}

	static int get_header_size() {
		return sizeof(vertex_index);
	}
};

class vertex_offset
{
	off_t off;
public:
	vertex_offset() {
		off = 0;
	}

	void init(off_t off, const in_mem_vertex &v) {
		this->off = off;
	}

	off_t get_off() const {
		return off;
	}
};

class default_vertex_index: public vertex_index_temp<vertex_offset>
{
public:
	static default_vertex_index *load(const std::string &index_file) {
		default_vertex_index *ret
			= (default_vertex_index *) vertex_index_temp<vertex_offset>::load(
					index_file);
		ret->verify();
		return ret;
	}

	void verify() const {
		// Right now all other types of graphs use the default vertex index.
		assert(get_graph_header().get_graph_type() != graph_type::DIRECTED);
		assert(index_size == sizeof(default_vertex_index)
			+ get_num_vertices() * sizeof(vertex_offset));
	}
};

class directed_vertex_entry
{
	off_t off;
	vsize_t in_edges;
	vsize_t out_edges;
public:
	directed_vertex_entry() {
		off = 0;
		in_edges = 0;
		out_edges = 0;
	}

	directed_vertex_entry(off_t off, int in_edges, int out_edges) {
		this->off = off;
		this->in_edges = in_edges;
		this->out_edges = out_edges;
	}

	void init(off_t off, const in_mem_vertex &v) {
		this->off = off;
		this->in_edges = v.get_num_edges(edge_type::IN_EDGE);
		this->out_edges = v.get_num_edges(edge_type::OUT_EDGE);
	}

	off_t get_off() const {
		return off;
	}

	size_t get_num_in_edges() const {
		return in_edges;
	}

	size_t get_num_out_edges() const {
		return out_edges;
	}
};

class directed_in_mem_vertex_index;

class directed_vertex_index: public vertex_index_temp<directed_vertex_entry>
{
	directed_vertex_index() {
	}
public:
	static directed_vertex_index *load(const std::string &index_file) {
		directed_vertex_index *ret
			= (directed_vertex_index *) vertex_index_temp<directed_vertex_entry>::load(
					index_file);
		ret->verify();
		return ret;
	}

	void verify() const {
		assert(get_graph_header().get_graph_type() == graph_type::DIRECTED);
		assert(index_size == sizeof(directed_vertex_index)
			+ get_num_vertices() * sizeof(directed_vertex_entry));
	}

	size_t get_num_in_edges(vertex_id_t id) const {
		return this->get_vertex(id).get_num_in_edges();
	}

	size_t get_num_out_edges(vertex_id_t id) const {
		return this->get_vertex(id).get_num_out_edges();
	}

	friend class directed_in_mem_vertex_index;
};

static inline off_t get_vertex_off(vertex_index *index, vertex_id_t id)
{
	if (index->get_graph_header().get_graph_type() == graph_type::DIRECTED)
		return ((directed_vertex_index *) index)->get_vertex_off(id);
	else
		return ((default_vertex_index *) index)->get_vertex_off(id);
}

static inline size_t get_vertex_size(vertex_index *index, vertex_id_t id)
{
	if (index->get_graph_header().get_graph_type() == graph_type::DIRECTED)
		return ((directed_vertex_index *) index)->get_vertex_size(id);
	else
		return ((default_vertex_index *) index)->get_vertex_size(id);
}

/**
 * These are the in-mem counterparts of vertex index above.
 */

class in_mem_vertex_index
{
public:
	virtual void add_vertex(char *ext_mem_vertex) = 0;
	virtual void dump(const std::string &file, const graph_header &header) = 0;
};

class directed_in_mem_vertex_index: public in_mem_vertex_index
{
	size_t tot_size;
	std::vector<directed_vertex_entry> vertices;
public:
	directed_in_mem_vertex_index() {
		tot_size = sizeof(graph_header);
	}

	virtual void add_vertex(char *ext_mem_vertex) {
		ext_mem_directed_vertex *v = (ext_mem_directed_vertex *) ext_mem_vertex;
		vertices.push_back(directed_vertex_entry(tot_size, v->get_num_in_edges(),
					v->get_num_out_edges()));
		assert(v->get_id() + 1 == vertices.size());
		tot_size += v->get_size();
	}

	virtual void dump(const std::string &file, const graph_header &header) {
		directed_vertex_index index;
		index.header = header;
		index.tot_size = tot_size;
		index.index_size = directed_vertex_index::get_header_size()
			+ sizeof(directed_vertex_entry) * vertices.size();

		FILE *f = fopen(file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			assert(0);
		}
		ssize_t ret = fwrite(&index, directed_vertex_index::get_header_size(),
				1, f);
		assert(ret == 1);
		ret = fwrite(vertices.data(),
				sizeof(directed_vertex_entry) * vertices.size(), 1, f);
		assert(ret == 1);
		fclose(f);
	}
};

template<class entry_type>
class vertex_index_iterator
{
	static const int BUF_SIZE = 1024 * 1024;
	size_t index_file_size;
	size_t graph_file_size;
	FILE *f;

	entry_type entry_buf[BUF_SIZE];
	off_t idx;
	size_t num_entries;
	vertex_id_t curr_id;

	vertex_index_iterator(const vertex_index_iterator &);
	vertex_index_iterator &operator=(const vertex_index_iterator &);

	void read_block() {
		size_t read_size = min(index_file_size - ftell(f), sizeof(entry_buf));
		assert(read_size % sizeof(entry_type) == 0);
		size_t ret = fread(entry_buf, read_size, 1, f);
		assert(ret == 1);
		idx = 0;
		num_entries = read_size / sizeof(entry_type);
	}

	vertex_index_iterator(const std::string &index_file) {
		curr_id = 0;
		memset(entry_buf, 0, sizeof(entry_buf));
		native_file nf(index_file);
		index_file_size = nf.get_size();
		f = fopen(index_file.c_str(), "r");
		assert(f);
		int header_size = vertex_index_temp<entry_type>::get_header_size();
		char header_buf[header_size];
		size_t ret = fread(header_buf, header_size, 1, f);
		assert(ret == 1);
		vertex_index_temp<entry_type> *header
			= (vertex_index_temp<entry_type> *) header_buf;
		graph_file_size = header->get_graph_size();
		read_block();
	}

	~vertex_index_iterator() {
		fclose(f);
	}
public:
	static vertex_index_iterator<entry_type> *create(
			const std::string &index_file) {
		return new vertex_index_iterator(index_file);
	}

	static void destroy(vertex_index_iterator<entry_type> *it) {
		delete it;
	}

	in_mem_vertex_info next() {
		assert((size_t) idx < num_entries);
		if ((size_t) idx == num_entries - 1) {
			// End of the file
			if (index_file_size - ftell(f) == 0) {
				entry_type entry = entry_buf[idx++];
				vertex_id_t id = curr_id++;
				return in_mem_vertex_info(id, entry.get_off(),
						graph_file_size - entry.get_off());
			}
			else {
				// Read a block and compute the size of the vertx.
				entry_type entry1 = entry_buf[idx];
				read_block();
				entry_type entry2 = entry_buf[idx];
				vertex_id_t id = curr_id++;
				return in_mem_vertex_info(id, entry1.get_off(),
						entry2.get_off() - entry1.get_off());
			}
		}
		else {
			entry_type entry1 = entry_buf[idx++];
			entry_type entry2 = entry_buf[idx];
			vertex_id_t id = curr_id++;
			return in_mem_vertex_info(id, entry1.get_off(),
					entry2.get_off() - entry1.get_off());
		}
	}

	bool has_next() const {
		assert((size_t) idx <= num_entries);
		if ((size_t) idx == num_entries)
			return index_file_size - ftell(f) > 0;
		else
			return true;
	}
};

typedef vertex_index_iterator<directed_vertex_entry> directed_vertex_index_iterator;

#endif
