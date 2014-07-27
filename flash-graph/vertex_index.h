#ifndef __VERTEX_INDEX_H__
#define __VERTEX_INDEX_H__

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

	class destroy_index
	{
	public:
		void operator()(vertex_index *index) {
			free(index);
		}
	};

public:
	typedef std::shared_ptr<vertex_index> ptr;

	/*
	 * Load the vertex index from SAFS.
	 */
	static vertex_index::ptr safs_load(const std::string &index_file);
	/*
	 * Load the vertex index from the Linux filesystem.
	 */
	static vertex_index::ptr load(const std::string &index_file);

	static size_t get_header_size() {
		return sizeof(vertex_index);
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

	size_t get_index_size() const {
		return index_size;
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
	typedef std::shared_ptr<vertex_index_temp<vertex_entry_type> > ptr;

	static typename vertex_index_temp<vertex_entry_type>::ptr cast(
			vertex_index::ptr index) {
		return std::static_pointer_cast<vertex_index_temp<vertex_entry_type>,
			   vertex_index>(index);
	}

	template<class vertex_type>
	static typename vertex_index_temp<vertex_entry_type>::ptr create(
			const graph_header &header,
			const std::vector<vertex_type> &vertices) {
		void *addr = malloc(sizeof(vertex_index_temp<vertex_entry_type>)
				+ sizeof(vertex_entry_type) * vertices.size());
		assert(addr);
		assert(header.get_num_vertices() == vertices.size());
		vertex_index_temp<vertex_entry_type>::ptr index(
				new (addr) vertex_index_temp<vertex_entry_type>(vertices.size()),
				destroy_index());
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

	template<class vertex_type>
	static typename vertex_index_temp<vertex_entry_type>::ptr create(
			const graph_header &header,
			const std::map<vertex_id_t, vertex_type> &vertices) {
		assert(!vertices.empty());
		vsize_t num_vertices = vertices.rbegin()->first + 1;
		bool has_edge_data = vertices.rbegin()->second.has_edge_data();
		void *addr = malloc(sizeof(vertex_index_temp<vertex_entry_type>)
				+ sizeof(vertex_entry_type) * num_vertices);
		assert(addr);
		assert(header.get_num_vertices() == num_vertices);
		vertex_index_temp<vertex_entry_type>::ptr index(
				new (addr) vertex_index_temp<vertex_entry_type>(num_vertices),
				destroy_index());
		index->header = header;
		assert(sizeof(index->header) == PAGE_SIZE);
		// All data of adjacency lists are stored after the header.
		size_t tot_size = sizeof(header);
		vertex_id_t id = 0;
		for (typename std::map<vertex_id_t, vertex_type>::const_iterator it
				= vertices.begin(); it != vertices.end(); it++) {
			const vertex_type &v = it->second;
			assert(v.has_edge_data() == has_edge_data);
			while (id < v.get_id()) {
				vertex_type empty_v(id, has_edge_data);
				index->vertices[id].init(tot_size, empty_v);
				tot_size += empty_v.get_serialize_size();
				id++;
			}
			assert(id == v.get_id());
			index->vertices[id].init(tot_size, v);
			tot_size += v.get_serialize_size();
			id++;
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

	vertex_entry_type *get_data() {
		return vertices;
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

	vertex_offset(off_t off) {
		this->off = off;
	}

	void init(off_t off) {
		this->off = off;
	}

	void init(off_t off, const in_mem_vertex &v) {
		this->off = off;
	}

	off_t get_off() const {
		return off;
	}
};

class default_in_mem_vertex_index;

class default_vertex_index: public vertex_index_temp<vertex_offset>
{
public:
	typedef std::shared_ptr<default_vertex_index> ptr;

	static default_vertex_index::ptr cast(vertex_index::ptr index) {
		return std::static_pointer_cast<default_vertex_index, vertex_index>(
				index);
	}

	static default_vertex_index::ptr load(const std::string &index_file) {
		default_vertex_index::ptr ret
			= cast(vertex_index_temp<vertex_offset>::load(index_file));
		ret->verify();
		return ret;
	}

	void verify() const {
		// Right now all other types of graphs use the default vertex index.
		assert(get_graph_header().get_graph_type() != graph_type::DIRECTED);
		assert(index_size == sizeof(default_vertex_index)
			+ get_num_vertices() * sizeof(vertex_offset));
	}

	friend class default_in_mem_vertex_index;
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

	directed_vertex_entry(off_t off) {
		this->off = off;
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
	typedef std::shared_ptr<directed_vertex_index> ptr;

	static directed_vertex_index::ptr cast(vertex_index::ptr index) {
		return std::static_pointer_cast<directed_vertex_index, vertex_index>(
				index);
	}

	static directed_vertex_index::ptr load(const std::string &index_file) {
		directed_vertex_index::ptr ret
			= cast(vertex_index_temp<directed_vertex_entry>::load(index_file));
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
	virtual ~in_mem_vertex_index() {
	}

	virtual void add_vertex(char *ext_mem_vertex) = 0;
	virtual void dump(const std::string &file, const graph_header &header) = 0;
};

class default_in_mem_vertex_index: public in_mem_vertex_index
{
	size_t tot_size;
	std::vector<vertex_offset> vertices;
public:
	default_in_mem_vertex_index() {
		tot_size = sizeof(graph_header);
	}

	virtual vsize_t get_vertex_size(char *ext_mem_vertex) = 0;
	virtual vertex_id_t get_vertex_id(char *ext_mem_vertex) = 0;

	virtual void add_vertex(char *ext_mem_vertex) {
		vertex_offset off;
		off.init(tot_size);
		vertices.push_back(off);
		assert(get_vertex_id(ext_mem_vertex) + 1 == vertices.size());
		tot_size += get_vertex_size(ext_mem_vertex);
	}

	virtual void dump(const std::string &file, const graph_header &header) {
		default_vertex_index index;
		index.header = header;
		index.tot_size = tot_size;
		index.index_size = default_vertex_index::get_header_size()
			+ sizeof(vertex_offset) * vertices.size();

		FILE *f = fopen(file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			assert(0);
		}
		ssize_t ret = fwrite(&index, default_vertex_index::get_header_size(),
				1, f);
		assert(ret == 1);
		if (!vertices.empty()) {
			ret = fwrite(vertices.data(),
					sizeof(vertex_offset) * vertices.size(), 1, f);
			assert(ret == 1);
		}
		fclose(f);
	}
};

class undirected_in_mem_vertex_index: public default_in_mem_vertex_index
{
public:
	virtual vsize_t get_vertex_size(char *ext_mem_vertex) {
		ext_mem_undirected_vertex *v
			= (ext_mem_undirected_vertex *) ext_mem_vertex;
		return v->get_size();
	}

	virtual vertex_id_t get_vertex_id(char *ext_mem_vertex) {
		ext_mem_undirected_vertex *v
			= (ext_mem_undirected_vertex *) ext_mem_vertex;
		return v->get_id();
	}
};

class ts_directed_in_mem_vertex_index: public default_in_mem_vertex_index
{
public:
	virtual vsize_t get_vertex_size(char *ext_mem_vertex) {
		ts_ext_mem_directed_vertex *v
			= (ts_ext_mem_directed_vertex *) ext_mem_vertex;
		return v->get_size();
	}

	virtual vertex_id_t get_vertex_id(char *ext_mem_vertex) {
		ts_ext_mem_directed_vertex *v
			= (ts_ext_mem_directed_vertex *) ext_mem_vertex;
		return v->get_id();
	}
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
		if (!vertices.empty()) {
			ret = fwrite(vertices.data(),
					sizeof(directed_vertex_entry) * vertices.size(), 1, f);
			assert(ret == 1);
		}
		fclose(f);
	}
};

class vertex_index_iterator
{
public:
	virtual ~vertex_index_iterator() {
	}

	virtual in_mem_vertex_info next() = 0;
	virtual bool has_next() const = 0;
};

template<class entry_type>
class vertex_index_iterator_impl: public vertex_index_iterator
{
	static const int BUF_SIZE = 1024 * 1024;
	size_t index_file_size;
	size_t graph_file_size;
	FILE *f;

	entry_type entry_buf[BUF_SIZE];
	off_t idx;
	size_t num_entries;
	vertex_id_t curr_id;

	vertex_index_iterator_impl(const vertex_index_iterator_impl<entry_type> &);
	vertex_index_iterator_impl &operator=(
			const vertex_index_iterator_impl<entry_type> &);

	void read_block() {
		size_t read_size = min(index_file_size - ftell(f), sizeof(entry_buf));
		assert(read_size % sizeof(entry_type) == 0);
		if (read_size > 0) {
			size_t ret = fread(entry_buf, read_size, 1, f);
			assert(ret == 1);
		}
		idx = 0;
		num_entries = read_size / sizeof(entry_type);
	}

	vertex_index_iterator_impl(const std::string &index_file) {
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

	~vertex_index_iterator_impl() {
		fclose(f);
	}
public:
	static vertex_index_iterator_impl<entry_type> *create(
			const std::string &index_file) {
		return new vertex_index_iterator_impl<entry_type>(index_file);
	}

	static void destroy(vertex_index_iterator_impl<entry_type> *it) {
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

typedef vertex_index_iterator_impl<directed_vertex_entry> directed_vertex_index_iterator;
typedef vertex_index_iterator_impl<vertex_offset> default_vertex_index_iterator;

#endif
