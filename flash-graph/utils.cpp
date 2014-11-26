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

#include <unistd.h>
#ifdef USE_GZIP
#include <zlib.h>
#endif

#include <memory>
#include <algorithm>

#include <boost/foreach.hpp>
#ifdef USE_STXXL
#define STXXL_PARALLEL_MODE_EXPLICIT
#include <stxxl.h>
#endif
#if defined(_OPENMP)
#include <parallel/algorithm>
#endif

#include "thread.h"
#include "native_file.h"

#include "graph.h"
//#include "edge_type.h"
#include "vertex.h"
#include "in_mem_storage.h"
#include "utils.h"

using namespace safs;

namespace fg
{

static const int EDGE_LIST_BLOCK_SIZE = 16 * 1024 * 1024;
static const size_t SORT_BUF_SIZE = 1024L * 1024 * 1024 * 2;
static const vsize_t VERTEX_TASK_SIZE = 1024 * 128;

static int num_threads = 1;

class format_error: public std::exception
{
	std::string msg;
public:
	format_error(const std::string &msg) {
		this->msg = msg;
	}

	~format_error() throw() {
	}

	const char* what() const throw() {
		return msg.c_str();
	}
};

static bool is_compressed(const std::string &file_name)
{
	size_t pos = file_name.rfind(".gz");
	if (pos == std::string::npos)
		return false;
	return pos + 3 == file_name.length();
}

template<class edge_data_type>
struct comp_edge {
	bool operator() (const edge<edge_data_type> &e1, const edge<edge_data_type> &e2) const {
		if (e1.get_from() == e2.get_from())
			return e1.get_to() < e2.get_to();
		else
			return e1.get_from() < e2.get_from();
	}

	static edge<edge_data_type> min_value() {
		vertex_id_t min_id = std::numeric_limits<vertex_id_t>::min();
		assert(min_id == 0);
		return edge<edge_data_type>(min_id, min_id);
	}

	static edge<edge_data_type> max_value() {
		vertex_id_t max_id = std::numeric_limits<vertex_id_t>::max();
		assert(max_id == INVALID_VERTEX_ID);
		return edge<edge_data_type>(max_id, max_id);
	}
};

template<>
struct comp_edge<ts_edge_data> {
	comp_edge() {
		printf("compare timestamp edge\n");
	}

	bool operator() (const edge<ts_edge_data> &e1, const edge<ts_edge_data> &e2) const {
		if (e1.get_from() != e2.get_from())
			return e1.get_from() < e2.get_from();
		else if (e1.get_data().get_timestamp() != e2.get_data().get_timestamp())
			return e1.get_data().get_timestamp() < e2.get_data().get_timestamp();
		else
			return e1.get_to() < e2.get_to();
	}

	static edge<ts_edge_data> min_value() {
		vertex_id_t min_id = std::numeric_limits<vertex_id_t>::min();
		time_t min_time = std::numeric_limits<time_t>::min();
		return edge<ts_edge_data>(min_id, min_id, ts_edge_data(min_time));
	}

	static edge<ts_edge_data> max_value() {
		vertex_id_t max_id = std::numeric_limits<vertex_id_t>::max();
		time_t max_time = std::numeric_limits<time_t>::max();
		return edge<ts_edge_data>(max_id, max_id, ts_edge_data(max_time));
	}
};

template<class edge_data_type>
struct comp_in_edge {
	bool operator() (const edge<edge_data_type> &e1, const edge<edge_data_type> &e2) const {
		if (e1.get_to() == e2.get_to())
			return e1.get_from() < e2.get_from();
		else
			return e1.get_to() < e2.get_to();
	}

	static edge<edge_data_type> min_value() {
		vertex_id_t min_id = std::numeric_limits<vertex_id_t>::min();
		assert(min_id == 0);
		return edge<edge_data_type>(min_id, min_id);
	}

	static edge<edge_data_type> max_value() {
		vertex_id_t max_id = std::numeric_limits<vertex_id_t>::max();
		assert(max_id == INVALID_VERTEX_ID);
		return edge<edge_data_type>(max_id, max_id);
	}
};

template<>
struct comp_in_edge<ts_edge_data> {
	comp_in_edge() {
		printf("compare timestamp in-edge\n");
	}

	bool operator() (const edge<ts_edge_data> &e1, const edge<ts_edge_data> &e2) const {
		if (e1.get_to() != e2.get_to())
			return e1.get_to() < e2.get_to();
		else if (e1.get_data().get_timestamp() != e2.get_data().get_timestamp())
			return e1.get_data().get_timestamp() < e2.get_data().get_timestamp();
		else
			return e1.get_from() < e2.get_from();
	}

	static edge<ts_edge_data> min_value() {
		vertex_id_t min_id = std::numeric_limits<vertex_id_t>::min();
		time_t min_time = std::numeric_limits<time_t>::min();
		return edge<ts_edge_data>(min_id, min_id, ts_edge_data(min_time));
	}

	static edge<ts_edge_data> max_value() {
		vertex_id_t max_id = std::numeric_limits<vertex_id_t>::max();
		time_t max_time = std::numeric_limits<time_t>::max();
		return edge<ts_edge_data>(max_id, max_id, ts_edge_data(max_time));
	}
};

template<class edge_data_type>
class edge_vector
{
protected:
	class bulk_iterator {
	public:
		typedef std::shared_ptr<bulk_iterator> ptr;

		virtual bool has_next() const = 0;
		virtual int fetch(std::vector<edge<edge_data_type> > &edges) = 0;
	};
public:
	class edge_stream {
		typedef std::vector<edge<edge_data_type> > edge_buf_t;
		edge_buf_t buf;
		typename edge_buf_t::iterator it;
		typename bulk_iterator::ptr bulk_it;
	public:
		edge_stream(typename bulk_iterator::ptr bulk_it) {
			this->bulk_it = bulk_it;
			it = buf.end();
		}

		bool empty() {
			// If we have edges in the buffer, it's not empty.
			if (it != buf.end())
				return false;
			// If the buffer is empty and the bulk iterator is empty,
			// the edge stream is empty.
			else if (!bulk_it->has_next())
				return true;
			// Otherwise, the edge stream isn't empty, and we should fetch
			// some edges from the bulk iterator and keep them in the buffer.
			else {
				buf.clear();
				bulk_it->fetch(buf);
				it = buf.begin();
				return false;
			}
		}

		edge<edge_data_type> operator*() const {
			return *it;
		}

		edge_stream &operator++() {
			++it;
			return *this;
		}

		const edge<edge_data_type> *operator->() const {
			return it.operator->();
		}
	};

	typedef std::shared_ptr<edge_vector<edge_data_type> > ptr;

	virtual void push_back(const edge<edge_data_type> &e) = 0;
	virtual void append(const std::vector<edge<edge_data_type> > &vec) = 0;
	virtual void sort(bool out_edge) = 0;
	virtual edge_stream get_stream() const = 0;
	virtual ptr clone() const = 0;
	virtual size_t size() const = 0;
	virtual bool empty() const = 0;
	virtual const edge<edge_data_type> &back() const = 0;
};

#ifdef USE_STXXL

template<class edge_data_type>
class stxxl_edge_vector: public edge_vector<edge_data_type>
{
	typename stxxl::VECTOR_GENERATOR<edge<edge_data_type>, 4, 2, 32 * 1024 * 1024>::result data;
	typedef typename stxxl::VECTOR_GENERATOR<edge<edge_data_type>, 4, 2,
			32 * 1024 * 1024>::result::const_iterator const_iterator;

	class stxxl_iterator: public edge_vector<edge_data_type>::bulk_iterator
	{
		stxxl::stream::vector_iterator2stream<typename stxxl_edge_vector<edge_data_type>::const_iterator> strm;
	public:
		stxxl_iterator(typename stxxl_edge_vector<edge_data_type>::const_iterator begin,
				typename stxxl_edge_vector<edge_data_type>::const_iterator end): strm(
					begin, end) {
		}

		virtual bool has_next() const {
			return !strm.empty();
		}

		virtual int fetch(std::vector<edge<edge_data_type> > &edges) {
			int i = 0;
			for (; i < 1024 && !strm.empty(); i++, ++strm)
				edges.push_back(*strm);
			return i;
		}
	};
public:
	stxxl_edge_vector() {
	}

	stxxl_edge_vector(const stxxl_edge_vector<edge_data_type> &vec): data(vec.data) {
	}

	virtual void push_back(const edge<edge_data_type> &e) {
		data.push_back(e);
	}

	void append(const std::vector<edge<edge_data_type> > &vec) {
		typename std::vector<edge<edge_data_type> >::const_iterator it = vec.begin();
		typename std::vector<edge<edge_data_type> >::const_iterator end = vec.end();
		for (; it != end; it++)
			this->push_back(*it);
	}

	void sort(bool out_edge) {
		if (out_edge) {
			comp_edge<edge_data_type> edge_comparator;
			stxxl::sort(data.begin(), data.end(), edge_comparator, SORT_BUF_SIZE);
		}
		else {
			comp_in_edge<edge_data_type> in_edge_comparator;
			stxxl::sort(data.begin(), data.end(), in_edge_comparator, SORT_BUF_SIZE);
		}
	}

	virtual typename edge_vector<edge_data_type>::edge_stream get_stream() const {
		typename edge_vector<edge_data_type>::bulk_iterator::ptr it(
				new stxxl_iterator(data.begin(), data.end()));
		return typename edge_vector<edge_data_type>::edge_stream(it);
	}

	virtual typename edge_vector<edge_data_type>::ptr clone() const {
		return typename edge_vector<edge_data_type>::ptr(
				new stxxl_edge_vector<edge_data_type>(*this));
	}

	virtual size_t size() const {
		return data.size();
	}

	virtual bool empty() const {
		return data.empty();
	}

	virtual const edge<edge_data_type> &back() const {
		return data.back();
	}
};

#endif

template<class edge_data_type>
class std_edge_vector: public edge_vector<edge_data_type>
{
	std::vector<edge<edge_data_type> > data;

	class std_iterator: public edge_vector<edge_data_type>::bulk_iterator
	{
		typedef typename std::vector<edge<edge_data_type> >::const_iterator std_edge_iterator;
		std_edge_iterator it;
		std_edge_iterator end;
	public:
		std_iterator(std_edge_iterator begin, std_edge_iterator end) {
			this->it = begin;
			this->end = end;
		}

		virtual bool has_next() const {
			return it != end;
		}

		virtual int fetch(std::vector<edge<edge_data_type> > &edges) {
			int i = 0;
			for (; i < 1024 && it != end; i++, it++)
				edges.push_back(*it);
			return i;
		}
	};
public:
	typedef typename std::vector<edge<edge_data_type> >::const_iterator const_iterator;

	std_edge_vector() {
	}

	std_edge_vector(const std_edge_vector<edge_data_type> &vec): data(vec.data) {
	}

	virtual void push_back(const edge<edge_data_type> &e) {
		data.push_back(e);
	}

	void append(const edge_vector<edge_data_type> &vec) {
		typename edge_vector<edge_data_type>::edge_stream strm
			= vec.get_stream();
		while (!strm.empty()) {
			data.push_back(*strm);
			++strm;
		}
	}

	void append(const std::vector<edge<edge_data_type> > &vec) {
		data.insert(data.end(), vec.begin(), vec.end());
	}

	void sort(bool out_edge) {
		if (out_edge) {
			comp_edge<edge_data_type> edge_comparator;
#if defined(_OPENMP)
			__gnu_parallel::sort(data.begin(), data.end(), edge_comparator);
#else
			std::sort(data.begin(), data.end(), edge_comparator);
#endif
		}
		else {
			comp_in_edge<edge_data_type> in_edge_comparator;
#if defined(_OPENMP)
			__gnu_parallel::sort(data.begin(), data.end(), in_edge_comparator);
#else
			std::sort(data.begin(), data.end(), in_edge_comparator);
#endif
		}
	}

	virtual typename edge_vector<edge_data_type>::edge_stream get_stream() const {
		typename edge_vector<edge_data_type>::bulk_iterator::ptr it(
				new std_iterator(data.begin(), data.end()));
		return typename edge_vector<edge_data_type>::edge_stream(it);
	}

	virtual typename edge_vector<edge_data_type>::ptr clone() const {
		return typename edge_vector<edge_data_type>::ptr(
				new std_edge_vector<edge_data_type>(*this));
	}

	virtual size_t size() const {
		return data.size();
	}

	virtual bool empty() const {
		return data.empty();
	}

	virtual const edge<edge_data_type> &back() const {
		return data.back();
	}

	const_iterator cbegin() const {
		return data.cbegin();
	}

	const_iterator cend() const {
		return data.cend();
	}
};

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
	delete index;
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
};

class mem_graph_store
{
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

	char *reset() {
		char *tmp = buf;
		buf = NULL;
		buf_cap = 0;
		buf_bytes = 0;
		return tmp;
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

class disk_directed_graph: public disk_serial_graph
{
	large_writer::ptr in_f;
	large_writer::ptr out_f;
	embedded_array<char> buf;
	std::string tmp_in_graph_file;
	std::string tmp_out_graph_file;

	void check_ext_graph(const edge_graph &edge_g,
			const in_mem_cdirected_vertex_index &idx,
			large_reader::ptr reader) const;
	void check_ext_graph(const edge_graph &edge_g,
			const directed_vertex_index &idx, large_reader::ptr reader) const;
public:
	disk_directed_graph(const edge_graph &g,
			large_io_creator::ptr creator): disk_serial_graph(
			new directed_in_mem_vertex_index(), g.get_edge_data_size(), creator) {
		tmp_in_graph_file = basename(tempnam(".", "in-directed"));
		in_f = creator->create_writer(tmp_in_graph_file);
		BOOST_VERIFY(in_f->seek(sizeof(graph_header), SEEK_SET) == sizeof(graph_header));
		tmp_out_graph_file = basename(tempnam(".", "out-directed"));
		out_f = creator->create_writer(tmp_out_graph_file);
	}

	~disk_directed_graph() {
		if (in_f) {
			in_f->delete_file();
			in_f = NULL;
		}
		if (out_f) {
			out_f->delete_file();
			out_f = NULL;
		}
	}

	virtual bool is_directed() const {
		return true;
	}

	virtual void check_ext_graph(const edge_graph &edge_g,
			const std::string &index_file, large_reader::ptr reader) const;

	virtual void add_vertex(const in_mem_vertex &v) {
		throw unsupported_exception();
#if 0
		serial_graph::add_vertex(v);

		assert(in_f);
		assert(out_f);
		int size = v.get_serialize_size(IN_EDGE);
		buf.resize(size);
		ext_mem_undirected_vertex::serialize(v, buf.data(), size, IN_EDGE);
		BOOST_VERIFY(fwrite(buf.data(), size, 1, in_f) == 1);

		size = v.get_serialize_size(OUT_EDGE);
		buf.resize(size);
		ext_mem_undirected_vertex::serialize(v, buf.data(), size, OUT_EDGE);
		BOOST_VERIFY(fwrite(buf.data(), size, 1, out_f) == 1);
#endif
	}

	void add_vertices(const serial_subgraph &subg) {
		const directed_serial_subgraph &d_subg = (const directed_serial_subgraph &) subg;
		for (size_t i = 0; i < d_subg.get_num_vertices(); i++)
			serial_graph::add_vertex(d_subg.get_vertex_info(i));
		printf("write %ld vertices (%ld, %ld)\n", d_subg.get_num_vertices(),
				d_subg.get_in_size(), d_subg.get_out_size());
		BOOST_VERIFY(in_f->write(d_subg.get_in_buf(), d_subg.get_in_size())
				== (ssize_t) d_subg.get_in_size());
		BOOST_VERIFY(out_f->write(d_subg.get_out_buf(), d_subg.get_out_size())
				== (ssize_t) d_subg.get_out_size());
	}

	void copy_file(large_reader::ptr reader, size_t from_size, large_writer::ptr to) {
		const size_t BUF_SIZE = 128 * 1024 * 1024;
		std::unique_ptr<char[]> buf = std::unique_ptr<char[]>(new char[BUF_SIZE]);
		size_t remain_size = from_size;
		size_t read_size = std::min(remain_size, BUF_SIZE);
		while (read_size > 0) {
			size_t ret = reader->read(buf.get(), read_size);
			BOOST_VERIFY(ret == read_size);
			ret = to->write(buf.get(), read_size);
			BOOST_VERIFY(ret == read_size);
			remain_size -= read_size;
			read_size = std::min(remain_size, BUF_SIZE);
		}
	}

	virtual void finalize_graph_file(const std::string &adj_file) {
		size_t out_size = out_f->get_write_bytes();
		out_f = NULL;

		large_reader::ptr reader = get_creator()->create_reader(tmp_out_graph_file);
		assert(reader);
		copy_file(reader, out_size, in_f);
		reader = NULL;
		large_writer::ptr writer = get_creator()->create_writer(tmp_out_graph_file);
		writer->delete_file();
		writer = NULL;

		// Write the real graph header.
		graph_header header(get_graph_type(), this->get_num_vertices(),
				this->get_num_edges(), this->get_edge_data_size());
		BOOST_VERIFY(in_f->seek(0, SEEK_SET) == 0);
		BOOST_VERIFY(in_f->write((char *) &header, sizeof(header)) == sizeof(header));
		BOOST_VERIFY(in_f->rename2(adj_file) == 0);
		in_f = NULL;
	}

	virtual graph_type get_graph_type() const {
		return graph_type::DIRECTED;
	}
};

class disk_undirected_graph: public disk_serial_graph
{
	large_writer::ptr f;
	embedded_array<char> buf;
	std::string tmp_graph_file;

	void check_ext_graph(const edge_graph &edge_g,
			const in_mem_cundirected_vertex_index &idx,
			large_reader::ptr reader) const;
	void check_ext_graph(const edge_graph &edge_g,
			const default_vertex_index &idx, large_reader::ptr reader) const;
public:
	disk_undirected_graph(const edge_graph &g,
			large_io_creator::ptr creator): disk_serial_graph(
			new undirected_in_mem_vertex_index(), g.get_edge_data_size(), creator) {
		tmp_graph_file = basename(tempnam(".", "undirected"));
		f = creator->create_writer(tmp_graph_file);
		BOOST_VERIFY(f->seek(sizeof(graph_header), SEEK_SET) == sizeof(graph_header));
	}

	~disk_undirected_graph() {
		if (f) {
			f->delete_file();
			f = NULL;
		}
	}

	virtual bool is_directed() const {
		return false;
	}

	virtual void check_ext_graph(const edge_graph &edge_g,
			const std::string &index_file, large_reader::ptr reader) const;

	virtual size_t get_num_edges() const {
		return serial_graph::get_num_edges() / 2;
	}

	virtual void add_vertex(const in_mem_vertex &v) {
		throw unsupported_exception();
#if 0
		serial_graph::add_vertex(v);
		assert(f);
		int size = v.get_serialize_size(IN_EDGE);
		buf.resize(size);
		ext_mem_undirected_vertex::serialize(v, buf.data(), size, IN_EDGE);
		BOOST_VERIFY(fwrite(buf.data(), size, 1, f) == 1);
#endif
	}

	virtual void finalize_graph_file(const std::string &adj_file) {
		// Write the real graph header.
		graph_header header(get_graph_type(), this->get_num_vertices(),
				this->get_num_edges(), this->get_edge_data_size());
		BOOST_VERIFY(f->seek(0, SEEK_SET) == 0);
		BOOST_VERIFY(f->write((char *) &header, sizeof(header)) == sizeof(header));
		if (f->rename2(adj_file) != 0) {
			fprintf(stderr, "can't rename %s to %s: %s\n",
					tmp_graph_file.c_str(), adj_file.c_str(), strerror(errno));
			exit(1);
		}
		f = NULL;
	}

	void add_vertices(const serial_subgraph &subg) {
		const undirected_serial_subgraph &u_subg = (const undirected_serial_subgraph &) subg;
		for (size_t i = 0; i < u_subg.get_num_vertices(); i++)
			serial_graph::add_vertex(u_subg.get_vertex_info(i));
		if (f->write(u_subg.get_buf(), u_subg.get_size()) != (ssize_t) u_subg.get_size()) {
			fprintf(stderr, "fail to write %ld bytes for %ld vertices: %s\n",
					u_subg.get_size(), u_subg.get_num_vertices(), strerror(errno));
			exit(1);
		}
	}

	virtual graph_type get_graph_type() const {
		return graph_type::UNDIRECTED;
	}
};

class mem_directed_graph: public mem_serial_graph
{
	mem_graph_store in_store;
	mem_graph_store out_store;
public:
	mem_directed_graph(size_t edge_data_size): mem_serial_graph(
			new directed_in_mem_vertex_index(), edge_data_size), in_store(
			graph_header::get_header_size()) {
	}

	virtual bool is_directed() const {
		return true;
	}

	virtual void add_vertex(const in_mem_vertex &v) {
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
			new undirected_in_mem_vertex_index(), edge_data_size), store(
			graph_header::get_header_size()) {
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

template<class edge_data_type = empty_data>
class undirected_edge_graph: public edge_graph
{
	typedef std::vector<edge<edge_data_type> > edge_list_t;
	typedef typename edge_vector<edge_data_type>::edge_stream edge_stream_t;

	std::vector<std::shared_ptr<edge_vector<edge_data_type> > > edge_lists;

	off_t add_edges(const edge_vector<edge_data_type> &edges, off_t idx,
			vertex_id_t id, std::vector<edge<edge_data_type> > &v_edges) const;

	void read_edges(edge_stream_t &, vertex_id_t until_id, edge_list_t &v_edges) const;

	vertex_id_t get_max_vertex_id() const {
		vertex_id_t max_id = 0;
		for (size_t i = 0; i < edge_lists.size(); i++)
			max_id = std::max(edge_lists[i]->back().get_from(), max_id);
		return max_id;
	}

	serial_graph::ptr create_serial_graph(large_io_creator::ptr creator) const {
		if (creator == NULL)
			return serial_graph::ptr(new mem_undirected_graph(
						this->get_edge_data_size()));
		else
			return serial_graph::ptr(new disk_undirected_graph(*this, creator));
	}
public:
	/**
	 * num_edges tells the edge graph that there will be num_edges
	 * edges added to the graph.
	 */
	undirected_edge_graph(
			std::vector<std::shared_ptr<edge_vector<edge_data_type> > > &edge_lists,
			bool has_data): edge_graph(has_data) {
		this->edge_lists = edge_lists;
	}

	void sort_edges() {
		for (size_t i = 0; i < edge_lists.size(); i++)
			edge_lists[i]->sort(true);
	}

	size_t get_num_edges() const {
		size_t num_edges = 0;
		for (size_t i = 0; i < edge_lists.size(); i++)
			num_edges += edge_lists[i]->size();
		return num_edges / 2;
	}

	void check_vertices(
			const std::vector<ext_mem_undirected_vertex *> &vertices,
			bool in_part) const;
	virtual std::shared_ptr<serial_graph> serialize_graph(
			large_io_creator::ptr creator) const;
};

/**
 * This represents a directed graph in the form of edge list.
 * It maintains a sorted list of out-edges (sorted on the from vertices)
 * and a sorted list of in-edges (sorted on the to vertices).
 */
template<class edge_data_type = empty_data>
class directed_edge_graph: public edge_graph
{
	typedef std::vector<edge<edge_data_type> > edge_list_t;
	typedef typename edge_vector<edge_data_type>::edge_stream edge_stream_t;

	std::vector<std::shared_ptr<edge_vector<edge_data_type> > > in_edge_lists;
	std::vector<std::shared_ptr<edge_vector<edge_data_type> > > out_edge_lists;

	void read_out_edges(edge_stream_t &, vertex_id_t until_id, edge_list_t &v_edges) const;
	void read_in_edges(edge_stream_t &, vertex_id_t until_id, edge_list_t &v_edges) const;

	vertex_id_t get_max_vertex_id() const {
		vertex_id_t max_id = 0;
		for (size_t i = 0; i < out_edge_lists.size(); i++) {
			if (!out_edge_lists[i]->empty())
				max_id = std::max(out_edge_lists[i]->back().get_from(), max_id);
			if (!in_edge_lists[i]->empty())
				max_id = std::max(in_edge_lists[i]->back().get_to(), max_id);
		}
		return max_id;
	}

	serial_graph::ptr create_serial_graph(large_io_creator::ptr creator) const {
		if (creator == NULL)
			return serial_graph::ptr(new mem_directed_graph(
						this->get_edge_data_size()));
		else
			return serial_graph::ptr(new disk_directed_graph(*this, creator));
	}
public:
	/**
	 * num_edges tells the edge graph that there will be num_edges
	 * edges added to the graph.
	 */
	directed_edge_graph(
			std::vector<std::shared_ptr<edge_vector<edge_data_type> > > &edge_lists,
			bool has_data): edge_graph(has_data) {
		this->in_edge_lists = edge_lists;
		this->out_edge_lists.resize(edge_lists.size());
		for (size_t i = 0; i < edge_lists.size(); i++)
			this->out_edge_lists[i] = edge_lists[i]->clone();
	}

	void sort_edges() {
		for (size_t i = 0; i < in_edge_lists.size(); i++) {
			out_edge_lists[i]->sort(true);
			in_edge_lists[i]->sort(false);
		}
	}

	void check_vertices(
			const std::vector<ext_mem_undirected_vertex *> &vertices,
			bool in_part) const;
	virtual std::shared_ptr<serial_graph> serialize_graph(
			large_io_creator::ptr creator) const;

	size_t get_num_edges() const {
		size_t num_edges = 0;
		for (size_t i = 0; i < in_edge_lists.size(); i++)
			num_edges += in_edge_lists[i]->size();
		return num_edges;
	}
};

template<class edge_data_type>
off_t undirected_edge_graph<edge_data_type>::add_edges(
		const edge_vector<edge_data_type> &edges, off_t idx,
		vertex_id_t id, std::vector<edge<edge_data_type> > &v_edges) const
{
	if ((size_t) idx >= edges.size())
		return idx;

	assert(edges[idx].get_from() >= id);
	off_t num_edges = edges.size();
	while (idx < num_edges && edges[idx].get_from() == id) {
		v_edges.push_back(edges[idx++]);
	}
	return idx;
}

template<class edge_data_type>
void get_all_edges(
		const std::vector<std::shared_ptr<edge_vector<edge_data_type> > > &edge_lists,
		std_edge_vector<edge_data_type> &edges)
{
	BOOST_FOREACH(std::shared_ptr<edge_vector<edge_data_type> > vec, edge_lists) {
		edges.append(*vec);
	}
}

template<class edge_data_type>
void undirected_edge_graph<edge_data_type>::check_vertices(
		const std::vector<ext_mem_undirected_vertex *> &vertices, bool) const
{
	assert(!vertices.empty());
	std_edge_vector<edge_data_type> edges;
	get_all_edges(edge_lists, edges);
	edges.sort(true);
	typename std_edge_vector<edge_data_type>::const_iterator it
		= std::lower_bound(edges.cbegin(), edges.cend(),
				edge<edge_data_type>(vertices[0]->get_id(), 0),
				comp_edge<edge_data_type>());

	for (size_t i = 0; i < vertices.size(); i++) {
		ext_mem_undirected_vertex *v = vertices[i];
		for (size_t j = 0; j < v->get_num_edges(); j++, it++) {
			assert(it != edges.cend());
			edge<edge_data_type> e = *it;
			assert(v->get_neighbor(j) == e.get_to());
			assert(v->get_id() == e.get_from());
			if (v->has_edge_data())
				assert(v->get_edge_data<edge_data_type>(j) == e.get_data());
		}
	}
}

template<class edge_data_type>
void directed_edge_graph<edge_data_type>::check_vertices(
		const std::vector<ext_mem_undirected_vertex *> &vertices, bool in_part) const
{
	assert(!vertices.empty());
	typename std_edge_vector<edge_data_type>::const_iterator it;
	std_edge_vector<edge_data_type> edges;
	if (in_part) {
		get_all_edges(in_edge_lists, edges);
		edges.sort(false);
	}
	else {
		get_all_edges(out_edge_lists, edges);
		edges.sort(true);
	}

	if (in_part)
		it = std::lower_bound(edges.cbegin(), edges.cend(),
				edge<edge_data_type>(0, vertices[0]->get_id()),
				comp_in_edge<edge_data_type>());
	else
		it = std::lower_bound(edges.cbegin(), edges.cend(),
				edge<edge_data_type>(vertices[0]->get_id(), 0),
				comp_edge<edge_data_type>());

	for (size_t i = 0; i < vertices.size(); i++) {
		// Check in-edges
		if (in_part) {
			ext_mem_undirected_vertex *v = vertices[i];
			for (size_t j = 0; j < v->get_num_edges(); j++, it++) {
				assert(it != edges.cend());
				edge<edge_data_type> e = *it;
				assert(v->get_neighbor(j) == e.get_from());
				assert(v->get_id() == e.get_to());
				if (v->has_edge_data())
					assert(v->get_edge_data<edge_data_type>(j) == e.get_data());
			}
		}
		else {
			// Check out-edges
			ext_mem_undirected_vertex *v = vertices[i];
			for (size_t j = 0; j < v->get_num_edges(); j++, it++) {
				assert(it != edges.cend());
				edge<edge_data_type> e = *it;
				assert(v->get_id() == e.get_from());
				assert(v->get_neighbor(j) == e.get_to());
				if (v->has_edge_data())
					assert(v->get_edge_data<edge_data_type>(j) == e.get_data());
			}
		}
	}
}

static const size_t BUF_SIZE = 1024L * 1024 * 1024 * 32;

size_t cal_vertex_size(const std::vector<ext_mem_vertex_info> &infos)
{
	assert(!infos.empty());
	return infos.back().get_off() + infos.back().get_size()
		- infos.front().get_off();
}

std::unique_ptr<char[]> read_vertices(large_reader::ptr reader,
		const std::vector<ext_mem_vertex_info> &infos,
		std::vector<ext_mem_undirected_vertex *> &vertices)
{
	size_t size = cal_vertex_size(infos);
	std::unique_ptr<char[]> buf = std::unique_ptr<char[]>(new char[size]);
	off_t off_begin = infos.front().get_off();
	BOOST_VERIFY(reader->seek(off_begin, SEEK_SET) == off_begin);
	BOOST_VERIFY(reader->read(buf.get(), size) == (ssize_t) size);
	BOOST_FOREACH(ext_mem_vertex_info info, infos) {
		off_t rel_off = info.get_off() - off_begin;
		vertices.push_back((ext_mem_undirected_vertex *) (buf.get() + rel_off));
	}
	return buf;
}

template<class VertexIndexType, class GetInfoFunc>
size_t check_all_vertices(large_reader::ptr reader, const VertexIndexType &idx,
		GetInfoFunc func, const edge_graph &edge_g, bool in_part)
{
	size_t num_vertices = 0;
	std::vector<ext_mem_vertex_info> infos;
	infos.push_back(func(idx, 0));
	while (num_vertices < idx.get_num_vertices()) {
		while (cal_vertex_size(infos) < BUF_SIZE
				&& infos.back().get_id() < idx.get_num_vertices() - 1) {
			infos.push_back(func(idx, infos.back().get_id() + 1));
		}
		std::vector<ext_mem_undirected_vertex *> vertices;
		std::unique_ptr<char[]> buf = read_vertices(reader, infos, vertices);
		num_vertices += vertices.size();
		edge_g.check_vertices(vertices, in_part);
		vertex_id_t last_id = infos.back().get_id();
		infos.clear();
		if (last_id < idx.get_num_vertices() - 1) {
			infos.push_back(func(idx, last_id + 1));
			assert(num_vertices < idx.get_num_vertices());
		}
	}
	return num_vertices;
}

void disk_undirected_graph::check_ext_graph(const edge_graph &edge_g,
		const std::string &index_file, large_reader::ptr adj_reader) const
{
	BOOST_LOG_TRIVIAL(info) << "check the graph in the external memory";
	vertex_index::ptr idx;
	// TODO is there a better option?
	if (adj_reader->is_safs())
		idx = vertex_index::safs_load(index_file);
	else
		idx = vertex_index::load(index_file);
	if (idx->is_compressed()) {
		in_mem_cundirected_vertex_index::ptr cidx
			= in_mem_cundirected_vertex_index::create(*idx);
		check_ext_graph(edge_g, *cidx, adj_reader);
	}
	else
		check_ext_graph(edge_g, *default_vertex_index::cast(idx), adj_reader);
}

void disk_undirected_graph::check_ext_graph(const edge_graph &edge_g,
		const in_mem_cundirected_vertex_index &idx,
		large_reader::ptr reader) const
{
	class get_undirected_info_func {
		ext_mem_vertex_info buf_info;
	public:
		ext_mem_vertex_info operator()(const in_mem_cundirected_vertex_index &idx,
				vertex_id_t id) {
			if (!buf_info.is_valid()) {
				ext_mem_vertex_info info(id, idx.get_vertex(id).get_off(),
						idx.get_size(id));
				buf_info = info;
				return info;
			}
			else if (buf_info.get_id() == id)
				return buf_info;
			else if (buf_info.get_id() + 1 == id) {
				ext_mem_vertex_info info(id, buf_info.get_off()
						+ idx.get_size(id - 1), idx.get_size(id));
				buf_info = info;
				return info;
			}
			else {
				ext_mem_vertex_info info(id, idx.get_vertex(id).get_off(),
						idx.get_size(id));
				buf_info = info;
				return info;
			}
		}
	};

	size_t num_vertices = check_all_vertices(reader, idx,
			get_undirected_info_func(), edge_g, true);
	BOOST_LOG_TRIVIAL(info) << boost::format("%1% vertices are checked")
		% num_vertices;
}

void disk_undirected_graph::check_ext_graph(const edge_graph &edge_g,
		const default_vertex_index &idx, large_reader::ptr reader) const
{
	struct get_undirected_info_func {
		ext_mem_vertex_info operator()(const default_vertex_index &idx,
				vertex_id_t id) {
			return idx.get_vertex_info(id);
		}
	};
	size_t num_vertices = check_all_vertices(reader, idx,
			get_undirected_info_func(), edge_g, true);
	BOOST_LOG_TRIVIAL(info) << boost::format("%1% vertices are checked")
		% num_vertices;
}

void disk_directed_graph::check_ext_graph(const edge_graph &edge_g,
		const std::string &index_file, large_reader::ptr adj_reader) const
{
	BOOST_LOG_TRIVIAL(info) << "check the graph in the external memory";
	vertex_index::ptr idx;
	// TODO is there a better option?
	if (adj_reader->is_safs())
		idx = vertex_index::safs_load(index_file);
	else
		idx = vertex_index::load(index_file);
	if (idx->is_compressed()) {
		in_mem_cdirected_vertex_index::ptr cidx
			= in_mem_cdirected_vertex_index::create(*idx);
		check_ext_graph(edge_g, *cidx, adj_reader);
	}
	else
		check_ext_graph(edge_g, *directed_vertex_index::cast(idx), adj_reader);
}

void disk_directed_graph::check_ext_graph(const edge_graph &edge_g,
		const in_mem_cdirected_vertex_index &idx,
		large_reader::ptr reader) const
{
	class get_in_part_info_func {
		ext_mem_vertex_info buf_info;
	public:
		ext_mem_vertex_info operator()(const in_mem_cdirected_vertex_index &idx,
				vertex_id_t id) {
			if (!buf_info.is_valid()) {
				ext_mem_vertex_info info(id, idx.get_vertex(id).get_in_off(),
						idx.get_in_size(id));
				buf_info = info;
				return info;
			}
			else if (buf_info.get_id() == id)
				return buf_info;
			else if (buf_info.get_id() + 1 == id) {
				ext_mem_vertex_info info(id, buf_info.get_off()
						+ idx.get_in_size(id - 1), idx.get_in_size(id));
				buf_info = info;
				return info;
			}
			else {
				ext_mem_vertex_info info(id, idx.get_vertex(id).get_in_off(),
						idx.get_in_size(id));
				buf_info = info;
				return info;
			}
		}
	};

	class get_out_part_info_func {
		ext_mem_vertex_info buf_info;
	public:
		ext_mem_vertex_info operator()(const in_mem_cdirected_vertex_index &idx,
				vertex_id_t id) {
			if (!buf_info.is_valid()) {
				ext_mem_vertex_info info(id, idx.get_vertex(id).get_out_off(),
						idx.get_out_size(id));
				buf_info = info;
				return info;
			}
			else if (buf_info.get_id() == id)
				return buf_info;
			else if (buf_info.get_id() + 1 == id) {
				ext_mem_vertex_info info(id, buf_info.get_off()
						+ idx.get_out_size(id - 1), idx.get_out_size(id));
				buf_info = info;
				return info;
			}
			else {
				ext_mem_vertex_info info(id, idx.get_vertex(id).get_out_off(),
						idx.get_out_size(id));
				buf_info = info;
				return info;
			}
		}
	};

	size_t num_vertices = check_all_vertices(reader, idx, get_in_part_info_func(),
			edge_g, true);
	size_t num_vertices1 = check_all_vertices(reader, idx, get_out_part_info_func(),
			edge_g, false);
	assert(num_vertices == num_vertices1);
	BOOST_LOG_TRIVIAL(info) << boost::format("%1% vertices are checked")
		% num_vertices;
}

void disk_directed_graph::check_ext_graph(const edge_graph &edge_g,
		const directed_vertex_index &idx, large_reader::ptr reader) const
{
	struct get_in_part_info_func {
		ext_mem_vertex_info operator()(const directed_vertex_index &idx,
				vertex_id_t id) {
			return idx.get_vertex_info_in(id);
		}
	};
	struct get_out_part_info_func {
		ext_mem_vertex_info operator()(const directed_vertex_index &idx,
				vertex_id_t id) {
			return idx.get_vertex_info_out(id);
		}
	};
	size_t num_vertices = check_all_vertices(reader, idx, get_in_part_info_func(),
			edge_g, true);
	size_t num_vertices1 = check_all_vertices(reader, idx, get_out_part_info_func(),
			edge_g, false);
	assert(num_vertices == num_vertices1);
	BOOST_LOG_TRIVIAL(info) << boost::format("%1% vertices are checked")
		% num_vertices;
}

void disk_serial_graph::dump(const std::string &index_file,
		const std::string &graph_file, bool compressed_index)
{
	struct timeval start, end;
	gettimeofday(&start, NULL);

	// Write the adjacency lists to the graph file.
	finalize_graph_file(graph_file);
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info) << boost::format(
			"It takes %1% seconds to dump the graph") % time_diff(start, end);

	start = end;
	graph_header header(get_graph_type(), this->get_num_vertices(),
			this->get_num_edges(), this->get_edge_data_size());
	vertex_index::ptr index = get_index().dump(header, compressed_index);
	creator->create_writer(index_file)->write((const char *) index.get(),
			index->get_index_size());
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info) << boost::format(
			"It takes %1% seconds to dump the index") % time_diff(start, end);
}

template<class edge_data_type = empty_data>
edge_graph::ptr par_load_edge_list_text(
		const std::vector<std::string> &files, bool has_edge_data,
		bool directed);

class graph_file_io
{
public:
	typedef std::shared_ptr<graph_file_io> ptr;

	virtual ~graph_file_io() {
	}

	/**
	 * It read a text of an edge list roughly the size of the wanted bytes.
	 * The returned text may be a little more than the wanted bytes, but
	 * it's guaranteed that all lines are complete.
	 * The returned string ends with '\0'.
	 */
	virtual std::unique_ptr<char[]> read_edge_list_text(const size_t wanted_bytes,
			size_t &read_bytes) = 0;

	virtual bool eof() const = 0;
};

class text_graph_file_io: public graph_file_io
{
	FILE *f;
	ssize_t file_size;
public:
	text_graph_file_io(const std::string file) {
		f = fopen(file.c_str(), "r");
		if (f == NULL)
			ABORT_MSG(boost::format("fail to open %1%: %2%")
					% file % strerror(errno));
		native_file local_f(file);
		file_size = local_f.get_size();
	}

	~text_graph_file_io() {
		if (f)
			fclose(f);
	}

	std::unique_ptr<char[]> read_edge_list_text(const size_t wanted_bytes,
			size_t &read_bytes);

	bool eof() const {
		off_t curr_off = ftell(f);
		return file_size - curr_off == 0;
	}
};

#ifdef USE_GZIP
class gz_graph_file_io: public graph_file_io
{
	std::vector<char> prev_buf;
	size_t prev_buf_bytes;

	gzFile f;
public:
	gz_graph_file_io(const std::string file) {
		BOOST_LOG_TRIVIAL(info) << (std::string("read gz file: ") + file);
		f = gzopen(file.c_str(), "rb");
		prev_buf_bytes = 0;
		prev_buf.resize(PAGE_SIZE);
	}

	~gz_graph_file_io() {
		gzclose(f);
	}

	std::unique_ptr<char[]> read_edge_list_text(const size_t wanted_bytes,
			size_t &read_bytes);

	bool eof() const {
		return gzeof(f) && prev_buf_bytes == 0;
	}
};

std::unique_ptr<char[]> gz_graph_file_io::read_edge_list_text(
		const size_t wanted_bytes1, size_t &read_bytes)
{
	read_bytes = 0;
	size_t wanted_bytes = wanted_bytes1;
	size_t buf_size = wanted_bytes + PAGE_SIZE;
	char *buf = new char[buf_size];
	std::unique_ptr<char[]> ret_buf(buf);
	if (prev_buf_bytes > 0) {
		memcpy(buf, prev_buf.data(), prev_buf_bytes);
		buf += prev_buf_bytes;
		read_bytes += prev_buf_bytes;
		wanted_bytes -= prev_buf_bytes;
		prev_buf_bytes = 0;
	}

	if (!gzeof(f)) {
		int ret = gzread(f, buf, wanted_bytes + PAGE_SIZE);
		if (ret <= 0) {
			if (ret < 0 || !gzeof(f)) {
				BOOST_LOG_TRIVIAL(fatal) << gzerror(f, &ret);
				exit(1);
			}
		}

		if ((size_t) ret > wanted_bytes) {
			int i = 0;
			int over_read = ret - wanted_bytes;
			for (; i < over_read; i++) {
				if (buf[wanted_bytes + i] == '\n') {
					i++;
					break;
				}
			}
			read_bytes += wanted_bytes + i;
			buf += wanted_bytes + i;

			prev_buf_bytes = over_read - i;
			assert(prev_buf_bytes <= PAGE_SIZE);
			memcpy(prev_buf.data(), buf, prev_buf_bytes);
		}
		else
			read_bytes += ret;
	}
	// The line buffer must end with '\0'.
	assert(read_bytes < buf_size);
	ret_buf[read_bytes] = 0;
	return ret_buf;
}
#endif

std::unique_ptr<char[]> text_graph_file_io::read_edge_list_text(
		const size_t wanted_bytes, size_t &read_bytes)
{
	off_t curr_off = ftell(f);
	off_t off = curr_off + wanted_bytes;
	// After we just to the new location, we need to further read another
	// page to search for the end of a line. If there isn't enough data,
	// we can just read all remaining data.
	if (off + PAGE_SIZE < file_size) {
		int ret = fseek(f, off, SEEK_SET);
		if (ret < 0) {
			perror("fseek");
			return NULL;
		}

		char buf[PAGE_SIZE];
		ret = fread(buf, sizeof(buf), 1, f);
		if (ret != 1) {
			perror("fread");
			return NULL;
		}
		unsigned i;
		for (i = 0; i < sizeof(buf); i++)
			if (buf[i] == '\n')
				break;
		// A line shouldn't be longer than a page.
		assert(i != sizeof(buf));

		// We read a little more than asked to make sure that we read
		// the entire line.
		read_bytes = wanted_bytes + i + 1;

		// Go back to the original offset in the file.
		ret = fseek(f, curr_off, SEEK_SET);
		assert(ret == 0);
	}
	else {
		read_bytes = file_size - curr_off;
	}

	// The line buffer must end with '\0'.
	char *line_buf = new char[read_bytes + 1];
	BOOST_VERIFY(fread(line_buf, read_bytes, 1, f) == 1);
	line_buf[read_bytes] = 0;

	return std::unique_ptr<char[]>(line_buf);
}

struct edge_line
{
	vertex_id_t from;
	vertex_id_t to;
	std::string data;

	edge_line(vertex_id_t from, vertex_id_t to, std::string data) {
		this->from = from;
		this->to = to;
		this->data = data;
	}

	edge_line(vertex_id_t from, vertex_id_t to) {
		this->from = from;
		this->to = to;
	}
};

struct edge_line parse_line(char *line)
{
	int len = strlen(line);

	char *first = line;
	for (; isspace(*first); first++);
	if (!isdigit(*first))
		throw format_error(
				std::string("the first entry isn't a number: ") + first);

	char *second = first;
	for (; isdigit(*second); second++);
	*second = 0;
	long from = atol(first);
	assert(from >= 0 && from < MAX_VERTEX_ID);

	if (second - line == len)
		throw format_error(std::string("there isn't second entry: ") + line);
	second++;
	if (!isdigit(*second))
		throw format_error(
				std::string("the second entry isn't a number: ") + second);
	char *third = second;
	for (; isdigit(*third); third++);
	*third = 0;
	long to = atol(second);
	assert(to >= 0 && to < MAX_VERTEX_ID);

	if (third - line == len)
		return edge_line(from, to);
	else {
		third++;
		return edge_line(from, to, third);
	}
}

size_t parse_edge_list_line(char *line, edge<ts_edge_data> &e)
{
	if (line[0] == '#')
		return 0;
	struct edge_line res = parse_line(line);
	if (!isdigit(res.data[0]))
		throw format_error(std::string("the third entry isn't a number: ")
				+ res.data);
	time_t timestamp = atol(res.data.c_str());
	ts_edge_data data(timestamp);
	e = edge<ts_edge_data>(res.from, res.to, data);
	return 1;
}

int parse_edge_list_line(char *line, edge<edge_count> &e)
{
	if (line[0] == '#')
		return 0;
	struct edge_line res = parse_line(line);
	if (!isdigit(res.data[0]))
		throw format_error(std::string("the third entry isn't a number: ")
				+ res.data);
	edge_count c(atol(res.data.c_str()));
	e = edge<edge_count>(res.from, res.to, c);
	return 1;
}

int parse_edge_list_line(char *line, edge<> &e)
{
	if (line[0] == '#')
		return 0;
	struct edge_line res = parse_line(line);
	e = edge<>(res.from, res.to);
	return 1;
}

/**
 * Parse the edge list in the character buffer.
 * `size' doesn't include '\0'.
 */
template<class edge_data_type>
size_t parse_edge_list_text(char *line_buf, size_t size,
		std::vector<edge<edge_data_type> > &edges)
{
	char *line_end;
	char *line = line_buf;
	size_t num_edges = 0;
	while ((line_end = strchr(line, '\n'))) {
		assert(line_end - line_buf <= (ssize_t) size);
		*line_end = 0;
		if (*(line_end - 1) == '\r')
			*(line_end - 1) = 0;
		edge<edge_data_type> e;
		int num = parse_edge_list_line(line, e);
		if (num > 0)
			edges.push_back(e);
		num_edges += num;
		line = line_end + 1;
	}
	if (line - line_buf < (ssize_t) size) {
		edge<edge_data_type> e;
		int num = parse_edge_list_line(line, e);
		if (num > 0)
			edges.push_back(e);
		num_edges += num;
	}
	return num_edges;
}

template<class edge_data_type>
class text_edge_task: public thread_task
{
	std::unique_ptr<char[]> line_buf;
	size_t size;
	bool directed;
public:
	text_edge_task(std::unique_ptr<char[]> line_buf, size_t size, bool directed) {
		this->line_buf = std::move(line_buf);
		this->size = size;
		this->directed = directed;
	}

	void run() {
		std::vector<edge<edge_data_type> > edges;
		parse_edge_list_text(line_buf.get(), size, edges);
		edge_vector<edge_data_type> *local_edge_buf
			= (edge_vector<edge_data_type> *) thread::get_curr_thread()->get_user_data();
		local_edge_buf->append(edges);

		// For an undirected graph, we need to store each edge twice
		// and each copy is the reverse of the original edge.
		if (!directed) {
			BOOST_FOREACH(edge<edge_data_type> e, edges) {
				e.reverse_dir();
				local_edge_buf->push_back(e);
			}
		}
	}
};

template<class edge_data_type>
class text_edge_file_task: public thread_task
{
	std::string file_name;
public:
	text_edge_file_task(const std::string file_name) {
		this->file_name = file_name;
	}

	void run();
};

template<class edge_data_type>
void text_edge_file_task<edge_data_type>::run()
{
	graph_file_io::ptr io;
	if (is_compressed(file_name)) {
#ifdef USE_GZIP
		io = graph_file_io::ptr(new gz_graph_file_io(file_name));
#else
		BOOST_LOG_TRIVIAL(error) << "Doesn't support reading gz file";
		BOOST_LOG_TRIVIAL(error)
			<< "zlib is required to support reading gz file";
		exit(1);
#endif
	}
	else
		io = graph_file_io::ptr(new text_graph_file_io(file_name));

	edge_vector<edge_data_type> *local_edge_buf
		= (edge_vector<edge_data_type> *) thread::get_curr_thread()->get_user_data();
	while (!io->eof()) {
		size_t size = 0;
		std::unique_ptr<char[]> data = io->read_edge_list_text(
				EDGE_LIST_BLOCK_SIZE, size);

		std::vector<edge<edge_data_type> > edges;
		parse_edge_list_text(data.get(), size, edges);
		local_edge_buf->append(edges);
	}

	std::cout << boost::format("There are %1% edges in thread %2%\n")
		% local_edge_buf->size() % thread::get_curr_thread()->get_id();
}

template<class edge_data_type>
void directed_edge_graph<edge_data_type>::read_out_edges(edge_stream_t &stream,
		vertex_id_t until_id, std::vector<edge<edge_data_type> > &v_edges) const
{
	if (stream.empty())
		return;

	while (!stream.empty() && stream->get_from() < until_id) {
		v_edges.push_back(*stream);
		++stream;
	}
}

template<class edge_data_type>
void directed_edge_graph<edge_data_type>::read_in_edges(edge_stream_t &stream,
		vertex_id_t until_id, std::vector<edge<edge_data_type> > &v_edges) const
{
	if (stream.empty())
		return;

	while (!stream.empty() && stream->get_to() < until_id) {
		v_edges.push_back(*stream);
		++stream;
	}
}

template<class edge_data_type>
void undirected_edge_graph<edge_data_type>::read_edges(edge_stream_t &stream,
		vertex_id_t until_id, std::vector<edge<edge_data_type> > &v_edges) const
{
	if (stream.empty())
		return;

	while (!stream.empty() && stream->get_from() < until_id) {
		v_edges.push_back(*stream);
		++stream;
	}
}

class write_graph_thread: public thread
{
	typedef std::shared_ptr<serial_subgraph> subgraph_ptr;
	struct subgraph_comp {
		bool operator()(const subgraph_ptr &g1, const subgraph_ptr &g2) {
			return g1->get_start_id() > g2->get_start_id();
		}
	};

	std::vector<subgraph_ptr> added_subgraphs;
	std::priority_queue<subgraph_ptr, std::vector<subgraph_ptr>, subgraph_comp> subgraphs;
	pthread_spinlock_t lock;
	serial_graph &g;
	vertex_id_t curr_id;
	vertex_id_t max_id;
public:
	write_graph_thread(serial_graph &_g,
			vertex_id_t max_id): thread("write-thread", -1), g(_g) {
		curr_id = 0;
		this->max_id = max_id;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	void add_vertices(subgraph_ptr subg) {
		pthread_spin_lock(&lock);
		added_subgraphs.push_back(subg);
		pthread_spin_unlock(&lock);
		activate();
	}

	void run();
};

void write_graph_thread::run()
{
	do {
		std::vector<subgraph_ptr> copy;
		pthread_spin_lock(&lock);
		copy = added_subgraphs;
		added_subgraphs.clear();
		pthread_spin_unlock(&lock);
		if (copy.empty() && (subgraphs.empty() || subgraphs.top()->get_start_id() > curr_id)) {
			usleep(10000);
		}

		BOOST_FOREACH(subgraph_ptr subg, copy) {
			subgraphs.push(subg);
		}

		if (!subgraphs.empty())
			assert(subgraphs.top()->get_start_id() >= curr_id);
		while (!subgraphs.empty() && subgraphs.top()->get_start_id() == curr_id) {
			subgraph_ptr subg = subgraphs.top();
			g.add_vertices(*subg);
			subgraphs.pop();
			curr_id = subg->get_end_id();
		}
	} while (curr_id <= max_id);
	BOOST_LOG_TRIVIAL(info) << boost::format("write %1% vertices") % curr_id;
	stop();
}

template<class edge_data_type>
class construct_directed_vertex_task: public thread_task
{
	typedef std::vector<edge<edge_data_type> > edge_list_t;
	std::shared_ptr<edge_list_t> in_edges;
	std::shared_ptr<edge_list_t> out_edges;
	vertex_id_t start_id;
	vertex_id_t end_id;
	write_graph_thread &write_thread;
	bool has_edge_data;

	typename edge_list_t::const_iterator add_in_edges(
			typename edge_list_t::const_iterator it,
			typename edge_list_t::const_iterator end, vertex_id_t id,
			in_mem_directed_vertex<edge_data_type> &v) {
		if (it == end)
			return it;
		assert(it->get_to() >= id);
		while (it != end && it->get_to() == id) {
			v.add_in_edge(*it);
			it++;
		}
		return it;
	}

	typename edge_list_t::const_iterator add_out_edges(
			typename edge_list_t::const_iterator it,
			typename edge_list_t::const_iterator end, vertex_id_t id,
			in_mem_directed_vertex<edge_data_type> &v) {
		if (it == end)
			return it;
		assert(it->get_from() >= id);
		while (it != end && it->get_from() == id) {
			v.add_out_edge(*it);
			it++;
		}
		return it;
	}
public:
	construct_directed_vertex_task(write_graph_thread &_write_thread,
			bool has_edge_data, vertex_id_t start_id, vertex_id_t end_id,
			std::shared_ptr<edge_list_t> in_edges,
			std::shared_ptr<edge_list_t> out_edges): write_thread(_write_thread) {
		this->in_edges = in_edges;
		this->out_edges = out_edges;
		this->start_id = start_id;
		this->end_id = end_id;
		this->has_edge_data = has_edge_data;
	}

	void run() {
		comp_edge<edge_data_type> edge_comparator;
		comp_in_edge<edge_data_type> in_edge_comparator;
		std::sort(in_edges->begin(), in_edges->end(), in_edge_comparator);
		std::sort(out_edges->begin(), out_edges->end(), edge_comparator);

		std::shared_ptr<directed_serial_subgraph> subg
			= std::shared_ptr<directed_serial_subgraph>(new directed_serial_subgraph());
		typename edge_list_t::const_iterator in_it = in_edges->begin();
		typename edge_list_t::const_iterator out_it = out_edges->begin();
		for (vertex_id_t id = start_id; id < end_id; id++) {
			in_mem_directed_vertex<edge_data_type> v(id, has_edge_data);
			in_it = add_in_edges(in_it, in_edges->end(), id, v);
			out_it = add_out_edges(out_it, out_edges->end(), id, v);
			subg->add_vertex(v);
		}
		write_thread.add_vertices(subg);
	}
};

template<class edge_data_type>
class construct_undirected_vertex_task: public thread_task
{
	typedef std::vector<edge<edge_data_type> > edge_list_t;
	std::shared_ptr<edge_list_t> edges;
	vertex_id_t start_id;
	vertex_id_t end_id;
	write_graph_thread &write_thread;
	bool has_edge_data;

	typename edge_list_t::const_iterator add_edges(
			typename edge_list_t::const_iterator it,
			typename edge_list_t::const_iterator end, vertex_id_t id,
			in_mem_undirected_vertex<edge_data_type> &v) {
		if (it == end)
			return it;
		assert(it->get_from() >= id);
		while (it != end && it->get_from() == id) {
			v.add_edge(*it);
			it++;
		}
		return it;
	}
public:
	construct_undirected_vertex_task(write_graph_thread &_write_thread,
			bool has_edge_data, vertex_id_t start_id, vertex_id_t end_id,
			std::shared_ptr<edge_list_t> edges): write_thread(_write_thread) {
		this->edges = edges;
		this->start_id = start_id;
		this->end_id = end_id;
		this->has_edge_data = has_edge_data;
	}

	void run() {
		comp_edge<edge_data_type> edge_comparator;
		std::sort(edges->begin(), edges->end(), edge_comparator);

		std::shared_ptr<undirected_serial_subgraph> subg
			= std::shared_ptr<undirected_serial_subgraph>(new undirected_serial_subgraph());
		typename edge_list_t::const_iterator it = edges->begin();
		for (vertex_id_t id = start_id; id < end_id; id++) {
			in_mem_undirected_vertex<edge_data_type> v(id, has_edge_data);
			it = add_edges(it, edges->end(), id, v);
			subg->add_vertex(v);
		}
		write_thread.add_vertices(subg);
	}
};

template<class edge_data_type>
serial_graph::ptr undirected_edge_graph<edge_data_type>::serialize_graph(
		large_io_creator::ptr creator) const
{
	struct timeval start, end;
	gettimeofday(&start, NULL);
	BOOST_LOG_TRIVIAL(info) << "start to serialize an undirected graph";
	serial_graph::ptr g = create_serial_graph(creator);
	std::vector<edge_stream_t> its;
	for (size_t i = 0; i < edge_lists.size(); i++)
		its.push_back(edge_lists[i]->get_stream());
	vertex_id_t max_id = get_max_vertex_id();

	std::vector<task_thread *> threads(num_threads);
	for (int i = 0; i < num_threads; i++) {
		task_thread *t = new task_thread(std::string(
					"graph-task-thread") + itoa(i), -1);
		t->start();
		threads[i] = t;
	}
	write_graph_thread *write_thread = new write_graph_thread(*g, max_id);
	write_thread->start();

	BOOST_LOG_TRIVIAL(info) << boost::format(
			"start to construct the graph. max id: %1%") % max_id;

	int thread_no = 0;
	for (vertex_id_t id = 0; id <= max_id; ) {
		std::shared_ptr<edge_list_t> v_edges
			= std::shared_ptr<edge_list_t>(new edge_list_t());
		vertex_id_t end_id = std::min(id + VERTEX_TASK_SIZE, max_id + 1);
		for (size_t i = 0; i < edge_lists.size(); i++)
			read_edges(its[i], end_id, *v_edges);

		construct_undirected_vertex_task<edge_data_type> *task
			= new construct_undirected_vertex_task<edge_data_type>(*write_thread,
					edge_graph::has_edge_data(), id, end_id, v_edges);
		threads[thread_no % num_threads]->add_task(task);
		thread_no++;
		id = end_id;
	}

	for (int i = 0; i < num_threads; i++) {
		threads[i]->wait4complete();
		threads[i]->stop();
		threads[i]->join();
		delete threads[i];
	}
	write_thread->join();
	delete write_thread;
	BOOST_LOG_TRIVIAL(info) << boost::format(
			"serial graph has %1% edges, edge graph has %2% edges")
			% g->get_num_edges() % get_num_edges();
	assert(g->get_num_edges() == get_num_edges());
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info) << boost::format(
			"It takes %1% to serialize an undirected graph") % time_diff(start, end);
	return g;
}

template<class edge_data_type>
serial_graph::ptr directed_edge_graph<edge_data_type>::serialize_graph(
		large_io_creator::ptr creator) const
{
	struct timeval start, end;
	gettimeofday(&start, NULL);
	BOOST_LOG_TRIVIAL(info) << "start to serialize a directed graph";
	serial_graph::ptr g = create_serial_graph(creator);
	assert(in_edge_lists.size() == out_edge_lists.size());
	for (size_t i = 0; i < in_edge_lists.size(); i++)
		assert(in_edge_lists[i]->size() == out_edge_lists[i]->size());

	std::vector<edge_stream_t> out_its;
	std::vector<edge_stream_t> in_its;
	for (size_t i = 0; i < out_edge_lists.size(); i++) {
		out_its.push_back(out_edge_lists[i]->get_stream());
		in_its.push_back(in_edge_lists[i]->get_stream());
	}
	vertex_id_t max_id = get_max_vertex_id();

	std::vector<task_thread *> threads(num_threads);
	for (int i = 0; i < num_threads; i++) {
		task_thread *t = new task_thread(std::string(
					"graph-task-thread") + itoa(i), -1);
		t->start();
		threads[i] = t;
	}
	write_graph_thread *write_thread = new write_graph_thread(*g, max_id);
	write_thread->start();

	BOOST_LOG_TRIVIAL(info) << boost::format(
				"start to construct the graph. max id: %1%") % max_id;

	int thread_no = 0;
	for (vertex_id_t id = 0; id <= max_id; ) {
		std::shared_ptr<edge_list_t> v_in_edges
			= std::shared_ptr<edge_list_t>(new edge_list_t());
		std::shared_ptr<edge_list_t> v_out_edges
			= std::shared_ptr<edge_list_t>(new edge_list_t());
		vertex_id_t end_id = std::min(id + VERTEX_TASK_SIZE, max_id + 1);
		for (size_t i = 0; i < in_edge_lists.size(); i++) {
			read_in_edges(in_its[i], end_id, *v_in_edges);
			read_out_edges(out_its[i], end_id, *v_out_edges);
		}

		construct_directed_vertex_task<edge_data_type> *task
			= new construct_directed_vertex_task<edge_data_type>(*write_thread,
					edge_graph::has_edge_data(), id, end_id, v_in_edges,
					v_out_edges);
		threads[thread_no % num_threads]->add_task(task);
		thread_no++;
		id = end_id;
	}

	for (int i = 0; i < num_threads; i++) {
		threads[i]->wait4complete();
		threads[i]->stop();
		threads[i]->join();
		delete threads[i];
	}
	write_thread->join();
	delete write_thread;
	assert(g->get_num_edges() == get_num_edges());
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info) << boost::format(
			"It takes %1% to serialize a directed graph") % time_diff(start, end);
	return g;
}

/**
 * This function loads edge lists from a tex file, parses them in parallel,
 * and convert the graph into the form of adjacency lists.
 */
template<class edge_data_type>
edge_graph::ptr par_load_edge_list_text(const std::vector<std::string> &files,
		bool has_edge_data, bool directed, bool in_mem)
{
	struct timeval start, end;
	gettimeofday(&start, NULL);
	BOOST_LOG_TRIVIAL(info) << "start to construct edge list";
	std::vector<task_thread *> threads(num_threads);
	for (int i = 0; i < num_threads; i++) {
		task_thread *t = new task_thread(std::string(
					"graph-task-thread") + itoa(i), -1);
		if (in_mem)
			t->set_user_data(new std_edge_vector<edge_data_type>());
		else {
#ifdef USE_STXXL
			t->set_user_data(new stxxl_edge_vector<edge_data_type>());
#else
			BOOST_LOG_TRIVIAL(error)
				<< "It doesn't support using disks to store intermediate results";
			BOOST_LOG_TRIVIAL(error)
				<< "stxxl is required to store intermediate results on disks";
			exit(1);
#endif
		}
		t->start();
		threads[i] = t;
	}
	int thread_no = 0;
	if (files.size() == 1) {
		const std::string file = files[0];
		BOOST_LOG_TRIVIAL(info) << (std::string(
					"start to read the edge list from ") + file.c_str());
		graph_file_io::ptr io;
		if (is_compressed(file)) {
#ifdef USE_GZIP
			io = graph_file_io::ptr(new gz_graph_file_io(file));
#else
			BOOST_LOG_TRIVIAL(error) << "Doesn't support reading gz file";
			BOOST_LOG_TRIVIAL(error)
				<< "zlib is required to support reading gz file";
			exit(1);
#endif
		}
		else
			io = graph_file_io::ptr(new text_graph_file_io(file));
		while (!io->eof()) {
			size_t size = 0;
			thread_task *task = new text_edge_task<edge_data_type>(
					io->read_edge_list_text(EDGE_LIST_BLOCK_SIZE, size),
					size, directed);
			threads[thread_no % num_threads]->add_task(task);
			thread_no++;
		}
	}
	else {
		for (size_t i = 0; i < files.size(); i++) {
			BOOST_LOG_TRIVIAL(info) << (std::string("read file " ) + files[i]);
			thread_task *task = new text_edge_file_task<edge_data_type>(files[i]);
			threads[thread_no % num_threads]->add_task(task);
			thread_no++;
		}
	}
	for (int i = 0; i < num_threads; i++)
		threads[i]->wait4complete();
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info) << boost::format(
			"It takes %1% seconds to construct edge list") % time_diff(start, end);
	start = end;
	BOOST_LOG_TRIVIAL(info) << "start to construct an edge graph";

	size_t num_edges = 0;
	std::vector<std::shared_ptr<edge_vector<edge_data_type> > > edge_lists(
			num_threads);
	for (int i = 0; i < num_threads; i++) {
		edge_vector<edge_data_type> *local_edges
			= (edge_vector<edge_data_type> *) threads[i]->get_user_data();
		num_edges += local_edges->size();
		edge_lists[i] = std::shared_ptr<edge_vector<edge_data_type> >(
				local_edges);
	}
	BOOST_LOG_TRIVIAL(info) << boost::format("There are %1% edges") % num_edges;

	edge_graph::ptr edge_g;
	if (directed)
		edge_g = edge_graph::ptr(new directed_edge_graph<edge_data_type>(
					edge_lists, has_edge_data));
	else
		edge_g = edge_graph::ptr(new undirected_edge_graph<edge_data_type>(
					edge_lists, has_edge_data));
	gettimeofday(&end, NULL);

	BOOST_LOG_TRIVIAL(info) << boost::format(
			"It takes %1% seconds to construct an edge graph") % time_diff(start, end);
	BOOST_LOG_TRIVIAL(info) << boost::format(
			"There are %1% edges in the edge graph") % edge_g->get_num_edges();

	for (int i = 0; i < num_threads; i++) {
		threads[i]->stop();
		threads[i]->join();
		delete threads[i];
	}

	return edge_g;
}

edge_graph::ptr parse_edge_lists(const std::vector<std::string> &edge_list_files,
		int edge_attr_type, bool directed, int nthreads, bool in_mem)
{
	BOOST_LOG_TRIVIAL(info) << "beofre load edge list";
	num_threads = nthreads;
	edge_graph::ptr g;
	switch(edge_attr_type) {
		case EDGE_COUNT:
			g = par_load_edge_list_text<edge_count>(edge_list_files, true,
					directed, in_mem);
			break;
		case EDGE_TIMESTAMP:
			g = par_load_edge_list_text<ts_edge_data>(edge_list_files, true,
					directed, in_mem);
			break;
		default:
			g = par_load_edge_list_text<empty_data>(edge_list_files, false,
					directed, in_mem);
	}
	return g;
}

serial_graph::ptr construct_graph(edge_graph::ptr edge_g,
		large_io_creator::ptr creator, int nthreads)
{
	BOOST_LOG_TRIVIAL(info) << "before sorting edges";
	num_threads = nthreads;
	struct timeval start, end;
	gettimeofday(&start, NULL);
	edge_g->sort_edges();
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info) << boost::format(
			"It takes %1% seconds to sort edge list") % time_diff(start, end);

	return edge_g->serialize_graph(creator);
}

std::pair<in_mem_graph::ptr, vertex_index::ptr> construct_mem_graph(
		const std::vector<std::string> &edge_list_files,
		const std::string &graph_name, int edge_attr_type, bool directed,
		int nthreads)
{
	edge_graph::ptr edge_g = parse_edge_lists(edge_list_files, edge_attr_type,
			directed, nthreads, true);
	serial_graph::ptr g = construct_graph(edge_g, large_io_creator::ptr(), nthreads);
	return std::pair<in_mem_graph::ptr, vertex_index::ptr>(
			((mem_serial_graph &) *g).dump_graph(graph_name), g->dump_index(true));
}

std::pair<in_mem_graph::ptr, vertex_index::ptr> construct_mem_graph(
		const std::vector<vertex_id_t> from, const std::vector<vertex_id_t> to,
		const std::string &graph_name, int edge_attr_type, bool directed,
		int num_threads)
{
	if (from.size() != to.size()) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"from vector (%1%) and to vector (%2%) have different length")
			% from.size() % to.size();
		return std::pair<in_mem_graph::ptr, vertex_index::ptr>();
	}

	size_t num_edges = from.size();
	std::vector<std::shared_ptr<edge_vector<empty_data> > > edge_lists(1);
	edge_lists[0] = std::shared_ptr<edge_vector<empty_data> >(
			new std_edge_vector<empty_data>());

	edge_graph::ptr edge_g;
	if (directed) {
		for (size_t i = 0; i < num_edges; i++)
			edge_lists[0]->push_back(edge<empty_data>(from[i], to[i]));
		edge_g = edge_graph::ptr(new directed_edge_graph<empty_data>(
					edge_lists, false));
	}
	else {
		for (size_t i = 0; i < num_edges; i++) {
			// Undirected edge graph assumes each edge has been added twice,
			// for both directions.
			edge_lists[0]->push_back(edge<empty_data>(from[i], to[i]));
			edge_lists[0]->push_back(edge<empty_data>(to[i], from[i]));
		}
		edge_g = edge_graph::ptr(new undirected_edge_graph<empty_data>(
					edge_lists, false));
	}
	serial_graph::ptr g = construct_graph(edge_g, large_io_creator::ptr(), num_threads);;
	return std::pair<in_mem_graph::ptr, vertex_index::ptr>(
			((mem_serial_graph &) *g).dump_graph(graph_name), g->dump_index(true));
}

}
