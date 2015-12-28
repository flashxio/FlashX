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

#include <unistd.h>
#ifdef USE_GZIP
#include <zlib.h>
#endif

#include <memory>
#include <algorithm>
#include <atomic>

#include <boost/foreach.hpp>
#if defined(_OPENMP)
#include <parallel/algorithm>
#endif

#include "thread.h"
#include "native_file.h"
#include "safs_exception.h"
#include "comm_exception.h"

#include "graph.h"
//#include "edge_type.h"
#include "vertex.h"
#include "in_mem_storage.h"
#include "utils.h"
#include "vertex_index_constructor.h"

using namespace safs;

namespace fg
{

namespace utils
{

static const int EDGE_LIST_BLOCK_SIZE = 16 * 1024 * 1024;
static const vsize_t VERTEX_TASK_SIZE = 1024 * 128;

void serial_graph::add_vertex(const in_mem_vertex &v)
{
	num_vertices++;
	// To get the total number of edges, I only accumulate on in-edges
	// or out-edges.
	num_edges += v.get_num_edges(edge_type::IN_EDGE);
	if (v.get_num_edges(edge_type::BOTH_EDGES) > 0)
		num_non_empty++;
	index->add_vertex(v);
}

vertex_index::ptr serial_graph::dump_index(bool compressed) const
{
	graph_header header(get_graph_type(), this->get_num_vertices(),
			this->get_num_edges(), this->get_edge_data_size());
	return index->dump(header, compressed);
}

serial_graph::~serial_graph()
{
}

class directed_vertex_info: public in_mem_vertex
{
	vertex_id_t id;
	int edge_data_size;
	size_t in_size;
	size_t out_size;
	size_t num_in_edges;
	size_t num_out_edges;
public:
	directed_vertex_info(const in_mem_vertex &v) {
		id = v.get_id();
		if (v.has_edge_data())
			edge_data_size = v.get_edge_data_size();
		else
			edge_data_size = 0;
		in_size = v.get_serialize_size(IN_EDGE);
		out_size = v.get_serialize_size(OUT_EDGE);
		num_in_edges = v.get_num_edges(IN_EDGE);
		num_out_edges = v.get_num_edges(OUT_EDGE);
	}

	virtual vertex_id_t get_id() const {
		return id;
	}
	virtual bool has_edge_data() const {
		return edge_data_size > 0;
	}
	virtual size_t get_edge_data_size() const {
		return edge_data_size;
	}
	virtual void serialize_edges(vertex_id_t ids[], edge_type type) const {
		ABORT_MSG("serialize_edges isn't implemented");
	}
	virtual void serialize_edge_data(char *data, edge_type type) const {
		ABORT_MSG("serialize_edge_data isn't implemented");
	}
	virtual size_t get_serialize_size(edge_type type) const {
		switch(type) {
			case IN_EDGE:
				return in_size;
			case OUT_EDGE:
				return out_size;
			case BOTH_EDGES:
				return in_size + out_size;
			default:
				ABORT_MSG("wrong edge type");
		}
	}
	virtual size_t get_num_edges(edge_type type) const {
		switch(type) {
			case IN_EDGE:
				return num_in_edges;
			case OUT_EDGE:
				return num_out_edges;
			case BOTH_EDGES:
				return num_in_edges + num_out_edges;
			default:
				ABORT_MSG("wrong edge type");
		}
	}

	in_mem_vertex::ptr create_remapped_vertex(
			const std::unordered_map<vertex_id_t, vertex_id_t> &map) const {
		ABORT_MSG("create_remapped_vertex isn't implemented");
	}

	void remap(const std::unordered_map<vertex_id_t, vertex_id_t> &map) {
		ABORT_MSG("remapped isn't implemented");
	}
};

class undirected_vertex_info: public in_mem_vertex
{
	vertex_id_t id;
	int edge_data_size;
	size_t size;
	size_t num_edges;
public:
	undirected_vertex_info(const in_mem_vertex &v) {
		id = v.get_id();
		if (v.has_edge_data())
			edge_data_size = v.get_edge_data_size();
		else
			edge_data_size = 0;
		size = v.get_serialize_size(OUT_EDGE);
		num_edges = v.get_num_edges(OUT_EDGE);
	}

	virtual vertex_id_t get_id() const {
		return id;
	}
	virtual bool has_edge_data() const {
		return edge_data_size > 0;
	}
	virtual size_t get_edge_data_size() const {
		return edge_data_size;
	}
	virtual void serialize_edges(vertex_id_t ids[], edge_type type) const {
		ABORT_MSG("serialize_edges isn't implemented");
	}
	virtual void serialize_edge_data(char *data, edge_type type) const {
		ABORT_MSG("serialize_edge_data isn't implemented");
	}
	virtual size_t get_serialize_size(edge_type type) const {
		return size;
	}
	virtual size_t get_num_edges(edge_type type) const {
		return num_edges;
	}

	in_mem_vertex::ptr create_remapped_vertex(
			const std::unordered_map<vertex_id_t, vertex_id_t> &map) const {
		ABORT_MSG("create_remapped_vertex isn't implemented");
	}

	void remap(const std::unordered_map<vertex_id_t, vertex_id_t> &map) {
		ABORT_MSG("remapped isn't implemented");
	}
};

class serial_subgraph
{
public:
	virtual ~serial_subgraph() {
	}

	virtual size_t get_num_vertices() const = 0;
	virtual vertex_id_t get_start_id() const = 0;
	virtual vertex_id_t get_end_id() const = 0;
	virtual size_t get_size() const = 0;
};

class mem_graph_store
{
	struct deleter {
		void operator()(char *buf) {
			delete [] buf;
		}
	};

	size_t buf_cap;
	size_t buf_bytes;
	char *buf;

	void expand_buf(size_t least_size) {
		while (buf_cap < least_size)
			buf_cap *= 2;
		char *tmp = new char[buf_cap];
		memcpy(tmp, buf, buf_bytes);
		delete [] buf;
		buf = tmp;
	}
public:
	mem_graph_store() {
		this->buf_cap = 1024 * 1024;
		buf_bytes = 0;
		buf = new char[buf_cap];
	}

	/*
	 * This constructor reserves some space in the memory buffer.
	 */
	mem_graph_store(size_t reserve) {
		this->buf_cap = 1024 * 1024;
		assert(reserve <= buf_cap);
		buf_bytes = reserve;
		buf = new char[buf_cap];
	}

	~mem_graph_store() {
		if (buf)
			delete [] buf;
	}

	void add_vertex(const in_mem_vertex &v, edge_type type) {
		int size = v.get_serialize_size(type);
		if (buf_bytes + size > buf_cap)
			expand_buf(buf_bytes + size);
		assert(buf_bytes + size <= buf_cap);
		ext_mem_undirected_vertex::serialize(v, buf + buf_bytes, size, type);
		buf_bytes += size;
	}

	size_t get_size() const {
		return buf_bytes;
	}

	const char *get_buf() const {
		return buf;
	}

	char *get_buf() {
		return buf;
	}

	void merge(const mem_graph_store &store) {
		if (buf_bytes + store.get_size() > buf_cap)
			expand_buf(buf_bytes + store.get_size());
		assert(buf_bytes + store.get_size() <= buf_cap);
		memcpy(buf + buf_bytes, store.get_buf(), store.get_size());
		buf_bytes += store.get_size();
	}

	std::shared_ptr<char> reset() {
		char *tmp = buf;
		buf = NULL;
		buf_cap = 0;
		buf_bytes = 0;
		return std::shared_ptr<char>(tmp, deleter());
	}
};

class directed_serial_subgraph: public serial_subgraph
{
	mem_graph_store in_store;
	mem_graph_store out_store;
	std::vector<directed_vertex_info> vertices;
public:
	void add_vertex(const in_mem_vertex &v) {
		if (!vertices.empty())
			assert(vertices.back().get_id() + 1 == v.get_id());
		vertices.push_back(directed_vertex_info(v));
		in_store.add_vertex(v, IN_EDGE);
		out_store.add_vertex(v, OUT_EDGE);
	}

	const directed_vertex_info &get_vertex_info(off_t idx) const {
		return vertices[idx];
	}

	size_t get_num_vertices() const {
		return vertices.size();
	}

	const char *get_in_buf() const {
		return in_store.get_buf();
	}

	size_t get_in_size() const {
		return in_store.get_size();
	}

	const mem_graph_store &get_in_store() const {
		return in_store;
	}

	const char *get_out_buf() const {
		return out_store.get_buf();
	}

	size_t get_out_size() const {
		return out_store.get_size();
	}

	size_t get_size() const {
		return get_in_size() + get_out_size();
	}

	const mem_graph_store &get_out_store() const {
		return out_store;
	}

	vertex_id_t get_start_id() const {
		assert(!vertices.empty());
		return vertices.front().get_id();
	}

	vertex_id_t get_end_id() const {
		assert(!vertices.empty());
		return vertices.back().get_id() + 1;
	}
};

class undirected_serial_subgraph: public serial_subgraph
{
	mem_graph_store store;
	std::vector<undirected_vertex_info> vertices;
public:
	void add_vertex(const in_mem_vertex &v) {
		if (!vertices.empty())
			assert(vertices.back().get_id() + 1 == v.get_id());
		vertices.push_back(undirected_vertex_info(v));
		store.add_vertex(v, OUT_EDGE);
	}

	const undirected_vertex_info &get_vertex_info(off_t idx) const {
		return vertices[idx];
	}

	size_t get_num_vertices() const {
		return vertices.size();
	}

	const char *get_buf() const {
		return store.get_buf();
	}

	size_t get_size() const {
		return store.get_size();
	}

	const mem_graph_store &get_store() const {
		return store;
	}

	vertex_id_t get_start_id() const {
		assert(!vertices.empty());
		return vertices.front().get_id();
	}

	vertex_id_t get_end_id() const {
		assert(!vertices.empty());
		return vertices.back().get_id() + 1;
	}
};

class mem_directed_graph: public mem_serial_graph
{
	mem_graph_store in_store;
	mem_graph_store out_store;
public:
	mem_directed_graph(size_t edge_data_size): mem_serial_graph(
			vertex_index_construct::create_compressed(true, edge_data_size),
			edge_data_size), in_store(graph_header::get_header_size()) {
	}

	virtual bool is_directed() const {
		return true;
	}

	virtual void add_vertex(const in_mem_vertex &v) {
		serial_graph::add_vertex(v);
		in_store.add_vertex(v, IN_EDGE);
		out_store.add_vertex(v, OUT_EDGE);
	}

	virtual void add_empty_vertex(vertex_id_t id) {
		in_mem_directed_vertex<> v(id, false);
		serial_graph::add_vertex(v);
		in_store.add_vertex(v, IN_EDGE);
		out_store.add_vertex(v, OUT_EDGE);
	}

	void add_vertices(const serial_subgraph &subg) {
		const directed_serial_subgraph &d_subg = (const directed_serial_subgraph &) subg;
		for (size_t i = 0; i < d_subg.get_num_vertices(); i++)
			serial_graph::add_vertex(d_subg.get_vertex_info(i));
		in_store.merge(d_subg.get_in_store());
		out_store.merge(d_subg.get_out_store());
	}

	virtual graph_type get_graph_type() const {
		return graph_type::DIRECTED;
	}

	in_mem_graph::ptr dump_graph(const std::string &graph_name) {
		graph_header header(get_graph_type(), this->get_num_vertices(),
				this->get_num_edges(), this->get_edge_data_size());
		memcpy(in_store.get_buf(), &header, graph_header::get_header_size());
		in_store.merge(out_store);
		size_t graph_size = in_store.get_size();
		in_mem_graph::ptr ret = in_mem_graph::create(graph_name,
				in_store.reset(), graph_size);
		out_store.reset();
		return ret;
	}
};

class mem_undirected_graph: public mem_serial_graph
{
	mem_graph_store store;
public:
	mem_undirected_graph(size_t edge_data_size): mem_serial_graph(
			vertex_index_construct::create_compressed(false, edge_data_size),
			edge_data_size), store(graph_header::get_header_size()) {
	}

	virtual bool is_directed() const {
		return false;
	}

	virtual size_t get_num_edges() const {
		return serial_graph::get_num_edges() / 2;
	}

	virtual void add_vertex(const in_mem_vertex &v) {
		serial_graph::add_vertex(v);
		store.add_vertex(v, OUT_EDGE);
	}

	virtual void add_empty_vertex(vertex_id_t id) {
		in_mem_undirected_vertex<> v(id, false);
		serial_graph::add_vertex(v);
		store.add_vertex(v, IN_EDGE);
	}

	void add_vertices(const serial_subgraph &subg) {
		const undirected_serial_subgraph &u_subg
			= (const undirected_serial_subgraph &) subg;
		for (size_t i = 0; i < u_subg.get_num_vertices(); i++)
			serial_graph::add_vertex(u_subg.get_vertex_info(i));
		store.merge(u_subg.get_store());
	}

	virtual graph_type get_graph_type() const {
		return graph_type::UNDIRECTED;
	}

	in_mem_graph::ptr dump_graph(const std::string &graph_name) {
		graph_header header(get_graph_type(), this->get_num_vertices(),
				this->get_num_edges(), this->get_edge_data_size());
		memcpy(store.get_buf(), &header, graph_header::get_header_size());
		size_t graph_size = store.get_size();
		in_mem_graph::ptr ret = in_mem_graph::create(graph_name,
				store.reset(), graph_size);
		return ret;
	}
};

mem_serial_graph::ptr mem_serial_graph::create(bool directed,
		size_t edge_data_size)
{
	if (directed)
		return mem_serial_graph::ptr(new mem_directed_graph(edge_data_size));
	else
		return mem_serial_graph::ptr(new mem_undirected_graph(edge_data_size));
}

}

}
