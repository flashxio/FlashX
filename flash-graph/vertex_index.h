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
#include <unordered_map>
#include <boost/format.hpp>

#include "native_file.h"

#include "vertex.h"
#include "graph_file_header.h"

class vertex_index
{
protected:
	union {
		struct {
			struct graph_header_struct header;
			size_t entry_size;
			size_t num_entries;
			off_t out_part_loc;

			// These are used for compressed vertx index.
			bool compressed;
			size_t num_large_in_vertices;
			size_t num_large_out_vertices;
		} data;
		char page[PAGE_SIZE];
	} h;

	vertex_index(size_t entry_size) {
		assert(sizeof(*this) == PAGE_SIZE);
		memset(this, 0, sizeof(*this));
		graph_header::init(h.data.header);
		h.data.entry_size = entry_size;
		h.data.num_entries = 0;
		h.data.out_part_loc = 0;

		h.data.compressed = false;
		h.data.num_large_in_vertices = 0;
		h.data.num_large_out_vertices = 0;
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

	const graph_header &get_graph_header() const {
		return (const graph_header &) *this;
	}

	size_t get_num_vertices() const {
		return h.data.header.num_vertices;
	}

	size_t get_num_entries() const {
		return h.data.num_entries;
	}

	vertex_id_t get_max_id() const {
		return get_num_vertices() - 1;
	}

	size_t get_index_size() const;

	off_t get_out_part_loc() const {
		return h.data.out_part_loc;
	}

	bool is_compressed() const {
		return h.data.compressed;
	}

	void dump(const std::string &file) const {
		FILE *f = fopen(file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			abort();
		}
		BOOST_VERIFY(fwrite(this, vertex_index::get_index_size(), 1, f));

		fclose(f);
	}
};

/**
 * This vertex index maps a vertex id to the location of the vertex in a file.
 */
template <class vertex_entry_type>
class vertex_index_temp: public vertex_index
{
	vertex_entry_type vertices[1];

protected:
	vertex_index_temp(const graph_header &header): vertex_index(
			sizeof(vertex_entry_type)) {
		memcpy(this, &header, sizeof(h.data.header));
		h.data.num_entries = 1;
		vertices[0] = vertex_entry_type(sizeof(graph_header));
	}
public:
	typedef std::shared_ptr<vertex_index_temp<vertex_entry_type> > ptr;

	static typename vertex_index_temp<vertex_entry_type>::ptr cast(
			vertex_index::ptr index) {
		return std::static_pointer_cast<vertex_index_temp<vertex_entry_type>,
			   vertex_index>(index);
	}

	static vertex_index::ptr create(const graph_header &header,
			const std::vector<vertex_entry_type> &vertices) {
		char *buf = (char *) malloc(vertex_index::get_header_size()
				+ vertices.size() * sizeof(vertices[0]));
		vertex_index_temp<vertex_entry_type> *index
			= new (buf) vertex_index_temp<vertex_entry_type>(header);
		index->h.data.num_entries = vertices.size();
		assert(header.get_num_vertices() + 1 == vertices.size());
		memcpy(buf + vertex_index::get_header_size(), vertices.data(),
				vertices.size() * sizeof(vertices[0]));
		return vertex_index::ptr(index);
	}

	static void dump(const std::string &file, const graph_header &header,
			const std::vector<vertex_entry_type> &vertices) {
		vertex_index_temp<vertex_entry_type> index(header);
		index.h.data.num_entries = vertices.size();
		assert(header.get_num_vertices() + 1 == vertices.size());
		FILE *f = fopen(file.c_str(), "w");
		if (f == NULL)
			ABORT_MSG(boost::format("fail to open %1%: %2%")
					% file % strerror(errno));

		BOOST_VERIFY(fwrite(&index, vertex_index::get_header_size(), 1, f));
		BOOST_VERIFY(fwrite(vertices.data(),
					vertices.size() * sizeof(vertices[0]), 1, f));

		fclose(f);
	}

	const vertex_entry_type &get_vertex(vertex_id_t id) const {
		assert(id < h.data.num_entries);
		return vertices[id];
	}

	vertex_entry_type *get_data() {
		return vertices;
	}

	const vertex_entry_type *get_data() const {
		return vertices;
	}

	size_t cal_index_size() const {
		return sizeof(vertex_index)
			+ h.data.num_entries * h.data.entry_size;
	}

	void verify() const {
		assert(h.data.entry_size == sizeof(vertex_entry_type));
		assert(h.data.header.num_vertices + 1 == h.data.num_entries);
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

	void init(const vertex_offset &prev, const in_mem_vertex &v) {
		this->off = prev.off + v.get_serialize_size(edge_type::IN_EDGE);
	}

	off_t get_off() const {
		return off;
	}
};

class default_vertex_index: public vertex_index_temp<vertex_offset>
{
public:
	typedef std::shared_ptr<default_vertex_index> ptr;

	static default_vertex_index::ptr cast(vertex_index::ptr index) {
		assert(!index->is_compressed());
		assert(!index->get_graph_header().is_directed_graph());
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
		vertex_index_temp<vertex_offset>::verify();
		// Right now all other types of graphs use the default vertex index.
		assert(get_graph_header().get_graph_type() != graph_type::DIRECTED);
		assert(get_vertex(0).get_off() == sizeof(graph_header));
		assert(h.data.out_part_loc == 0);
	}

	size_t get_graph_size() const {
		off_t last_idx = h.data.num_entries - 1;
		assert((size_t) last_idx == h.data.header.num_vertices);
		return get_vertex(last_idx).get_off();
	}

	ext_mem_vertex_info get_vertex_info(vertex_id_t id) const {
		off_t next_off = get_vertex(id + 1).get_off();
		off_t off = get_vertex(id).get_off();
		return ext_mem_vertex_info(id, off, next_off - off);
	}
};

class directed_vertex_entry
{
	off_t in_off;
	off_t out_off;
public:
	directed_vertex_entry() {
		in_off = 0;
		out_off = 0;
	}

	// When the vertex index is constructed, we assume in-part and out-part
	// of vertices are stored in two separate files.
	directed_vertex_entry(off_t off) {
		in_off = off;
		out_off = off;
	}

	directed_vertex_entry(off_t in_off, off_t out_off) {
		this->in_off = in_off;
		this->out_off = out_off;
	}

	void init(const directed_vertex_entry &prev, const in_mem_vertex &v) {
		this->in_off = prev.in_off + v.get_serialize_size(IN_EDGE);
		this->out_off = prev.out_off + v.get_serialize_size(OUT_EDGE);
	}

	off_t get_in_off() const {
		return in_off;
	}

	off_t get_out_off() const {
		return out_off;
	}
};

class directed_vertex_index: public vertex_index_temp<directed_vertex_entry>
{
	directed_vertex_index(const graph_header &header): vertex_index_temp<directed_vertex_entry>(header) {
	}
public:
	typedef std::shared_ptr<directed_vertex_index> ptr;

	static directed_vertex_index::ptr cast(vertex_index::ptr index) {
		assert(!index->is_compressed());
		assert(index->get_graph_header().is_directed_graph());
		return std::static_pointer_cast<directed_vertex_index, vertex_index>(
				index);
	}

	static directed_vertex_index::ptr load(const std::string &index_file) {
		directed_vertex_index::ptr ret
			= cast(vertex_index_temp<directed_vertex_entry>::load(index_file));
		ret->verify();
		return ret;
	}

	static vertex_index::ptr create(const graph_header &header,
			const std::vector<directed_vertex_entry> &vertices) {
		char *buf = (char *) malloc(vertex_index::get_header_size()
				+ vertices.size() * sizeof(vertices[0]));
		directed_vertex_index *index = new (buf) directed_vertex_index(header);
		index->h.data.num_entries = vertices.size();
		index->h.data.out_part_loc = vertices.front().get_out_off();
		assert(header.get_num_vertices() + 1 == vertices.size());
		memcpy(buf + vertex_index::get_header_size(), vertices.data(),
				vertices.size() * sizeof(vertices[0]));
		return vertex_index::ptr(index);
	}

	static void dump(const std::string &file, const graph_header &header,
			const std::vector<directed_vertex_entry> &vertices) {
		directed_vertex_index index(header);
		index.h.data.num_entries = vertices.size();
		index.h.data.out_part_loc = vertices.front().get_out_off();
		assert(header.get_num_vertices() + 1 == vertices.size());
		FILE *f = fopen(file.c_str(), "w");
		if (f == NULL)
			ABORT_MSG(boost::format("fail to open %1%: %2%")
					% file % strerror(errno));

		BOOST_VERIFY(fwrite(&index, vertex_index::get_header_size(), 1, f));
		BOOST_VERIFY(fwrite(vertices.data(),
					vertices.size() * sizeof(vertices[0]), 1, f));

		fclose(f);
	}

	void verify() const {
		vertex_index_temp<directed_vertex_entry>::verify();
		TEST(get_graph_header().get_graph_type() == graph_type::DIRECTED);
		TEST(get_vertex(0).get_in_off() == sizeof(graph_header));
		// All out-part of vertices are stored behind the in-part of vertices.
		TEST(get_vertex(0).get_out_off()
				== get_vertex(get_num_vertices()).get_in_off());
		TEST(h.data.out_part_loc == get_vertex(0).get_out_off());
	}

	size_t get_graph_size() const {
		off_t last_idx = h.data.num_entries - 1;
		assert((size_t) last_idx == h.data.header.num_vertices);
		return get_vertex(last_idx).get_out_off();
	}

	ext_mem_vertex_info get_vertex_info_in(vertex_id_t id) const {
		off_t next_off = get_vertex(id + 1).get_in_off();
		off_t off = get_vertex(id).get_in_off();
		return ext_mem_vertex_info(id, off, next_off - off);
	}

	ext_mem_vertex_info get_vertex_info_out(vertex_id_t id) const {
		off_t next_off = get_vertex(id + 1).get_out_off();
		off_t off = get_vertex(id).get_out_off();
		return ext_mem_vertex_info(id, off, next_off - off);
	}
};

/*
 * A compressed vertex index groups multiple vertices into a single data
 * structure, called compressed vertex entry. Each entry contains the location
 * of the first vertex in the entry and the number of edges of each vertex
 * in the entry.
 */

class compressed_vertex_entry
{
public:
	typedef unsigned char edge_size_t;
	static const size_t ENTRY_SIZE = 32;
	static const size_t LARGE_VERTEX_SIZE
		= std::numeric_limits<edge_size_t>::max();
};

/*
 * This defines the compressed vertex entry for a directed graph.
 */
class compressed_directed_vertex_entry: public compressed_vertex_entry
{
	typedef std::pair<edge_size_t, edge_size_t> edge_pair_t;
	directed_vertex_entry start_offs;
	edge_pair_t edges[ENTRY_SIZE];
public:
	compressed_directed_vertex_entry() {
		memset(this, 0, sizeof(*this));
	}

	compressed_directed_vertex_entry(const directed_vertex_entry offs[],
			size_t edge_data_size, size_t num);

	vsize_t get_num_in_edges(int idx) const {
		return edges[idx].first;
	}

	vsize_t get_num_out_edges(int idx) const {
		return edges[idx].second;
	}

	bool is_large_in_vertex(int idx) const {
		return edges[idx].first == LARGE_VERTEX_SIZE;
	}

	bool is_large_out_vertex(int idx) const {
		return edges[idx].second == LARGE_VERTEX_SIZE;
	}

	off_t get_start_in_off() const {
		return start_offs.get_in_off();
	}

	off_t get_start_out_off() const {
		return start_offs.get_out_off();
	}

	const directed_vertex_entry &get_start_offs() const {
		return start_offs;
	}
};

/*
 * This defines the compressed vertex entry for an undirected graph.
 */
class compressed_undirected_vertex_entry: public compressed_vertex_entry
{
	vertex_offset start_offs;
	edge_size_t edges[ENTRY_SIZE];
public:
	compressed_undirected_vertex_entry() {
		memset(this, 0, sizeof(*this));
	}

	compressed_undirected_vertex_entry(const vertex_offset offs[],
			size_t edge_data_size, size_t num);

	vsize_t get_num_edges(int idx) const {
		return edges[idx];
	}

	bool is_large_vertex(int idx) const {
		return edges[idx] == LARGE_VERTEX_SIZE;
	}

	off_t get_start_off() const {
		return start_offs.get_off();
	}

	const vertex_offset &get_start_offs() const {
		return start_offs;
	}
};

class default_in_mem_vertex_index;

/*
 * This class defines the data layout of a compressed vertex index for
 * an undirected graph in the external memory. It's not used to answer
 * queries in memory. The data layout of the index is
 *	header,
 *	vertex entries,
 *	large vertices.
 */
typedef std::pair<vertex_id_t, vsize_t> large_vertex_t;
class cundirected_vertex_index: public vertex_index
{
	static const size_t ENTRY_SIZE = compressed_undirected_vertex_entry::ENTRY_SIZE;
	compressed_undirected_vertex_entry entries[0];

	large_vertex_t *get_large_vertices() {
		return (large_vertex_t *) &entries[h.data.num_entries];
	}
public:
	typedef std::shared_ptr<cundirected_vertex_index> ptr;

	static typename cundirected_vertex_index::ptr cast(vertex_index::ptr index) {
		assert(index->is_compressed());
		assert(!index->get_graph_header().is_directed_graph());
		return std::static_pointer_cast<cundirected_vertex_index, vertex_index>(
				index);
	}

	static ptr construct(default_vertex_index &index);

	const compressed_undirected_vertex_entry *get_entries() const {
		return entries;
	}

	const large_vertex_t *get_large_vertices() const {
		return (const large_vertex_t *) &entries[h.data.num_entries];
	}

	size_t cal_index_size() const {
		return sizeof(cundirected_vertex_index)
		+ sizeof(entries[0]) * h.data.num_entries
		// We use num_large_in_vertices to store the number of large
		// undirected vertices.
		+ sizeof(large_vertex_t) * h.data.num_large_in_vertices;
	}

	size_t get_num_large_vertices() const {
		return h.data.num_large_in_vertices;
	}

	void verify() const {
		// We use num_large_in_vertices to store the number of large
		// undirected vertices and set num_large_out_vertices to 0.
		assert(h.data.num_large_out_vertices == 0);
		assert(h.data.entry_size == sizeof(compressed_undirected_vertex_entry));
		assert(ROUNDUP(h.data.header.num_vertices, ENTRY_SIZE) / ENTRY_SIZE
				== h.data.num_entries);
	}
};

class directed_in_mem_vertex_index;

/*
 * This class defines the data layout of a compressed vertex index for
 * a directed graph in the external memory. It's not used to answer
 * queries in memory. The data layout of the index is
 *	header,
 *	vertex entries,
 *	large in-vertices,
 *	large out-vertices.
 */
class cdirected_vertex_index: public vertex_index
{
	static const size_t ENTRY_SIZE = compressed_directed_vertex_entry::ENTRY_SIZE;
	compressed_directed_vertex_entry entries[0];

	large_vertex_t *get_large_in_vertices() {
		return (large_vertex_t *) &entries[h.data.num_entries];
	}

	large_vertex_t *get_large_out_vertices() {
		return get_large_in_vertices() + h.data.num_large_in_vertices;
	}
public:
	typedef std::shared_ptr<cdirected_vertex_index> ptr;

	static typename cdirected_vertex_index::ptr cast(vertex_index::ptr index) {
		assert(index->is_compressed());
		assert(index->get_graph_header().is_directed_graph());
		return std::static_pointer_cast<cdirected_vertex_index, vertex_index>(
				index);
	}

	static ptr construct(directed_vertex_index &index);

	const compressed_directed_vertex_entry *get_entries() const {
		return entries;
	}

	const large_vertex_t *get_large_in_vertices() const {
		return (const large_vertex_t *) &entries[h.data.num_entries];
	}

	const large_vertex_t *get_large_out_vertices() const {
		return get_large_in_vertices() + h.data.num_large_in_vertices;
	}

	size_t cal_index_size() const {
		return sizeof(cdirected_vertex_index)
		+ sizeof(entries[0]) * h.data.num_entries
		+ sizeof(large_vertex_t) * h.data.num_large_in_vertices
		+ sizeof(large_vertex_t) * h.data.num_large_out_vertices;
	}

	size_t get_num_large_in_vertices() const {
		return h.data.num_large_in_vertices;
	}

	size_t get_num_large_out_vertices() const {
		return h.data.num_large_out_vertices;
	}

	void verify() const {
		assert(h.data.entry_size == sizeof(compressed_directed_vertex_entry));
		assert(ROUNDUP(h.data.header.num_vertices, ENTRY_SIZE) / ENTRY_SIZE
				== h.data.num_entries);
	}
};

/*
 * This defines the interface that other parts of FlashGraph can query
 * the in-memory vertex index.
 */
class in_mem_query_vertex_index
{
	bool directed;
	bool compressed;
protected:
	in_mem_query_vertex_index(bool directed, bool compressed) {
		this->directed = directed;
		this->compressed = compressed;
	}
public:
	typedef std::shared_ptr<in_mem_query_vertex_index> ptr;
	static ptr create(vertex_index::ptr index, bool compress);

	bool is_directed() const {
		return directed;
	}

	bool is_compressed() const {
		return compressed;
	}

	virtual vsize_t get_num_edges(vertex_id_t id, edge_type type) const = 0;
	virtual vertex_index::ptr get_raw_index() const = 0;
};

/*
 * This is the in-memory compressed vertex index for an undirected graph.
 * This is used to answer queries from other parts of FlashGraph.
 */
class in_mem_cundirected_vertex_index: public in_mem_query_vertex_index
{
	static const size_t ENTRY_SIZE = compressed_undirected_vertex_entry::ENTRY_SIZE;
	static const size_t ENTRY_MASK = ENTRY_SIZE - 1;

	typedef std::unordered_map<vertex_id_t, vsize_t> vertex_map_t;

	size_t num_vertices;
	size_t edge_data_size;
	vertex_map_t large_vmap;
	std::vector<compressed_undirected_vertex_entry> entries;

	in_mem_cundirected_vertex_index(vertex_index &index);

	void init(const default_vertex_index &index);
	void init(const cundirected_vertex_index &index);
public:
	typedef std::shared_ptr<in_mem_cundirected_vertex_index> ptr;

	static ptr cast(in_mem_query_vertex_index::ptr index) {
		assert(index->is_compressed());
		assert(!index->is_directed());
		return std::static_pointer_cast<in_mem_cundirected_vertex_index,
			   in_mem_query_vertex_index>(index);
	}

	static ptr create(vertex_index &index) {
		return ptr(new in_mem_cundirected_vertex_index(index));
	}

	size_t get_num_vertices() const {
		return num_vertices;
	}

	vsize_t get_num_edges(vertex_id_t id, edge_type type) const {
		size_t entry_idx = id / ENTRY_SIZE;
		vsize_t num_edges = entries[entry_idx].get_num_edges(id % ENTRY_SIZE);
		if (num_edges >= compressed_undirected_vertex_entry::LARGE_VERTEX_SIZE) {
			vertex_map_t::const_iterator it = large_vmap.find(id);
			assert(it != large_vmap.end());
			num_edges = it->second;
		}
		return num_edges;
	}

	vertex_index::ptr get_raw_index() const {
		return vertex_index::ptr();
	}

	size_t get_size(vertex_id_t id) const {
		vsize_t num_edges = get_num_edges(id, edge_type::IN_EDGE);
		return ext_mem_undirected_vertex::num_edges2vsize(num_edges,
				edge_data_size);
	}

	vertex_offset get_vertex(vertex_id_t id) const;

	void verify_against(default_vertex_index &index);
};

/*
 * This is the in-memory compressed vertex index for a directed graph.
 * This is used to answer queries from other parts of FlashGraph.
 */
class in_mem_cdirected_vertex_index: public in_mem_query_vertex_index
{
	static const size_t ENTRY_SIZE = compressed_directed_vertex_entry::ENTRY_SIZE;
	static const size_t ENTRY_MASK = ENTRY_SIZE - 1;

	typedef std::unordered_map<vertex_id_t, vsize_t> vertex_map_t;

	size_t num_vertices;
	size_t edge_data_size;
	vertex_map_t large_in_vmap;
	vertex_map_t large_out_vmap;
	std::vector<compressed_directed_vertex_entry> entries;

	in_mem_cdirected_vertex_index(vertex_index &index);

	void init(const directed_vertex_index &index);
	void init(const cdirected_vertex_index &index);
public:
	typedef std::shared_ptr<in_mem_cdirected_vertex_index> ptr;

	static ptr cast(in_mem_query_vertex_index::ptr index) {
		assert(index->is_compressed());
		assert(index->is_directed());
		return std::static_pointer_cast<in_mem_cdirected_vertex_index,
			   in_mem_query_vertex_index>(index);
	}

	static ptr create(vertex_index &index) {
		return ptr(new in_mem_cdirected_vertex_index(index));
	}

	size_t get_num_vertices() const {
		return num_vertices;
	}

	vsize_t get_num_in_edges(vertex_id_t id) const {
		size_t entry_idx = id / ENTRY_SIZE;
		vsize_t num_edges = entries[entry_idx].get_num_in_edges(id % ENTRY_SIZE);
		if (num_edges >= compressed_directed_vertex_entry::LARGE_VERTEX_SIZE) {
			vertex_map_t::const_iterator it = large_in_vmap.find(id);
			assert(it != large_in_vmap.end());
			num_edges = it->second;
		}
		return num_edges;
	}

	vsize_t get_num_out_edges(vertex_id_t id) const {
		size_t entry_idx = id / ENTRY_SIZE;
		vsize_t num_edges = entries[entry_idx].get_num_out_edges(id % ENTRY_SIZE);
		if (num_edges >= compressed_directed_vertex_entry::LARGE_VERTEX_SIZE) {
			vertex_map_t::const_iterator it = large_out_vmap.find(id);
			assert(it != large_out_vmap.end());
			num_edges = it->second;
		}
		return num_edges;
	}

	virtual vsize_t get_num_edges(vertex_id_t id, edge_type type) const {
		switch (type) {
			case edge_type::IN_EDGE:
				return get_num_in_edges(id);
			case edge_type::OUT_EDGE:
				return get_num_out_edges(id);
			case edge_type::BOTH_EDGES:
				return get_num_in_edges(id) + get_num_out_edges(id);
			default:
				ABORT_MSG("wrong edge type");
		}
	}

	vertex_index::ptr get_raw_index() const {
		return vertex_index::ptr();
	}

	size_t get_in_size(vertex_id_t id) const {
		vsize_t num_edges = get_num_in_edges(id);
		return ext_mem_undirected_vertex::num_edges2vsize(num_edges,
				edge_data_size);
	}

	size_t get_out_size(vertex_id_t id) const {
		vsize_t num_edges = get_num_out_edges(id);
		return ext_mem_undirected_vertex::num_edges2vsize(num_edges,
				edge_data_size);
	}

	directed_vertex_entry get_vertex(vertex_id_t id) const;

	void verify_against(directed_vertex_index &index);
};

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

static inline int get_index_entry_size(graph_type type)
{
	switch (type) {
		case graph_type::UNDIRECTED:
			return sizeof(vertex_offset);
		case graph_type::DIRECTED:
			return sizeof(directed_vertex_entry);
		default:
			ABORT_MSG("wrong graph type");
	}
}

#endif
