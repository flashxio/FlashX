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
#include <zlib.h>

#include "graph.h"

#include <memory>
#include <algorithm>

#include <boost/foreach.hpp>
#include <stxxl.h>

#include "thread.h"
#include "native_file.h"

#include "edge_type.h"
#include "vertex.h"

template<class edge_data_type>
class stxxl_edge_vector: public stxxl::VECTOR_GENERATOR<edge<edge_data_type> >::result
{
public:
	typedef typename stxxl::VECTOR_GENERATOR<edge<edge_data_type> >::result::const_iterator const_iterator;
	void append(stxxl_edge_vector<edge_data_type>::const_iterator it,
			stxxl_edge_vector<edge_data_type>::const_iterator end) {
		for (; it != end; it++)
			this->push_back(*it);
	}

	void append(typename std::vector<edge<edge_data_type> >::const_iterator it,
			typename std::vector<edge<edge_data_type> >::const_iterator end) {
		for (; it != end; it++)
			this->push_back(*it);
	}
};

int num_threads = 1;
const int EDGE_LIST_BLOCK_SIZE = 16 * 1024 * 1024;
const size_t SORT_BUF_SIZE = 1024 * 1024 * 1024;
static const vsize_t VERTEX_TASK_SIZE = 1024 * 128;

bool decompress = false;
bool check_graph = false;
std::string work_dir = ".";

struct timeval start_time;

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

template<class edge_data_type = empty_data>
class edge_graph
{
	bool has_data;
public:
	edge_graph(bool has_data) {
		this->has_data = has_data;
	}

	virtual ~edge_graph() {
	}

	virtual void sort_edges() = 0;
#if 0
	virtual void check_vertices(
			const std::vector<ext_mem_undirected_vertex *> &vertices,
			bool in_part) const = 0;
#endif
	virtual void construct_graph(graph *g) const = 0;
	virtual graph *create_disk_graph() const = 0;
	virtual size_t get_num_edges() const = 0;

	bool has_edge_data() const {
		return has_data;
	}
};

class subgraph;

/**
 * This is a disk-backed directed graph.
 */
template<class edge_data_type>
class disk_graph: public graph
{
	size_t num_edges;
	size_t num_vertices;
	size_t num_non_empty;
	const edge_graph<edge_data_type> *g;
	in_mem_vertex_index *index;
public:
	disk_graph(const edge_graph<edge_data_type> *g, in_mem_vertex_index *index) {
		num_edges = 0;
		num_vertices = 0;
		num_non_empty = 0;
		this->g = g;
		this->index = index;
	}

	virtual ~disk_graph() {
		delete g;
		delete index;
	}

	virtual const in_mem_vertex &get_vertex(vertex_id_t id) const {
		ABORT_MSG("get_vertex isn't supported");
	}

	virtual void add_vertex(const in_mem_vertex &v) {
		num_vertices++;
		// To get the total number of edges, I only accumulate on in-edges
		// or out-edges.
		num_edges += v.get_num_edges(edge_type::IN_EDGE);
		if (v.get_num_edges(edge_type::BOTH_EDGES) > 0)
			num_non_empty++;
		index->add_vertex(v);
	}

	const edge_graph<edge_data_type> *get_edge_graph() const {
		return g;
	}

	virtual size_t get_num_edges() const {
		return num_edges;
	}

	virtual size_t get_num_vertices() const {
		return num_vertices;
	}

	virtual bool has_edge_data() const {
		return g->has_edge_data();
	}

	virtual size_t get_num_non_empty_vertices() const {
		return num_non_empty;
	}

	virtual void dump(const std::string &index_file,
			const std::string &graph_file);

	virtual vertex_index::ptr create_vertex_index() const {
		ABORT_MSG("don't support creating vertex index");
	}

	void get_all_vertices(std::vector<vertex_id_t> &ids) const {
		ABORT_MSG("don't support getting all vertices");
	}

	virtual void print() const {
		ABORT_MSG("don't support printing");
	}

#if 0
	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const = 0;
#endif
	virtual graph_type get_graph_type() const = 0;
	virtual void finalize_graph_file(const std::string &adj_file) = 0;
	virtual void add_vertices(const subgraph &subg) = 0;
};

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
};

class subgraph
{
public:
	virtual size_t get_num_vertices() const = 0;
	virtual vertex_id_t get_start_id() const = 0;
	virtual vertex_id_t get_end_id() const = 0;
};

class directed_subgraph: public subgraph
{
	size_t in_buf_cap;
	size_t out_buf_cap;
	size_t in_buf_bytes;
	size_t out_buf_bytes;
	char *in_buf;
	char *out_buf;
	std::vector<directed_vertex_info> vertices;

	void expand_in_buf(size_t least_size) {
		while (in_buf_cap < least_size)
			in_buf_cap *= 2;
		char *tmp = new char[in_buf_cap];
		memcpy(tmp, in_buf, in_buf_bytes);
		delete [] in_buf;
		in_buf = tmp;
	}

	void expand_out_buf(size_t least_size) {
		while (out_buf_cap < least_size)
			out_buf_cap *= 2;
		char *tmp = new char[out_buf_cap];
		memcpy(tmp, out_buf, out_buf_bytes);
		delete [] out_buf;
		out_buf = tmp;
	}
public:
	directed_subgraph() {
		this->in_buf_cap = 1024 * 1024;
		this->out_buf_cap = 1024 * 1024;
		in_buf_bytes = 0;
		out_buf_bytes = 0;
		in_buf = new char[in_buf_cap];
		out_buf = new char[out_buf_cap];
	}

	~directed_subgraph() {
		delete [] in_buf;
		delete [] out_buf;
	}

	template<class edge_data_type>
	void add_vertex(const in_mem_directed_vertex<edge_data_type> &v) {
		if (!vertices.empty())
			assert(vertices.back().get_id() + 1 == v.get_id());
		vertices.push_back(directed_vertex_info(v));

		int size = v.get_serialize_size(IN_EDGE);
		if (in_buf_bytes + size > in_buf_cap)
			expand_in_buf(in_buf_bytes + size);
		assert(in_buf_bytes + size <= in_buf_cap);
		ext_mem_undirected_vertex::serialize(v, in_buf + in_buf_bytes,
				size, IN_EDGE);
		in_buf_bytes += size;

		size = v.get_serialize_size(OUT_EDGE);
		if (out_buf_bytes + size > out_buf_cap)
			expand_out_buf(out_buf_bytes + size);
		assert(out_buf_bytes + size <= out_buf_cap);
		ext_mem_undirected_vertex::serialize(v, out_buf + out_buf_bytes,
				size, OUT_EDGE);
		out_buf_bytes += size;
	}

	const directed_vertex_info &get_vertex_info(off_t idx) const {
		return vertices[idx];
	}

	size_t get_num_vertices() const {
		return vertices.size();
	}

	char *get_in_buf() const {
		return in_buf;
	}

	size_t get_in_size() const {
		return in_buf_bytes;
	}

	char *get_out_buf() const {
		return out_buf;
	}

	size_t get_out_size() const {
		return out_buf_bytes;
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

class undirected_subgraph: public subgraph
{
	size_t buf_cap;
	size_t buf_bytes;
	char *buf;
	std::vector<undirected_vertex_info> vertices;

	void expand_buf(size_t least_size) {
		while (buf_cap < least_size)
			buf_cap *= 2;
		char *tmp = new char[buf_cap];
		memcpy(tmp, buf, buf_bytes);
		delete [] buf;
		buf = tmp;
	}
public:
	undirected_subgraph() {
		this->buf_cap = 1024 * 1024;
		buf_bytes = 0;
		buf = new char[buf_cap];
	}

	~undirected_subgraph() {
		delete [] buf;
	}

	template<class edge_data_type>
	void add_vertex(const in_mem_undirected_vertex<edge_data_type> &v) {
		if (!vertices.empty())
			assert(vertices.back().get_id() + 1 == v.get_id());
		vertices.push_back(undirected_vertex_info(v));

		int size = v.get_serialize_size(OUT_EDGE);
		if (buf_bytes + size > buf_cap)
			expand_buf(buf_bytes + size);
		assert(buf_bytes + size <= buf_cap);
		ext_mem_undirected_vertex::serialize(v, buf + buf_bytes,
				size, OUT_EDGE);
		buf_bytes += size;
	}

	const undirected_vertex_info &get_vertex_info(off_t idx) const {
		return vertices[idx];
	}

	size_t get_num_vertices() const {
		return vertices.size();
	}

	char *get_buf() const {
		return buf;
	}

	size_t get_size() const {
		return buf_bytes;
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

template<class edge_data_type>
class disk_directed_graph: public disk_graph<edge_data_type>
{
	FILE *in_f;
	FILE *out_f;
	embedded_array<char> buf;
	std::string tmp_in_graph_file;
	std::string tmp_out_graph_file;
public:
	disk_directed_graph(const edge_graph<edge_data_type> *g): disk_graph<edge_data_type>(
			g, new directed_in_mem_vertex_index()) {
		tmp_in_graph_file = tempnam(work_dir.c_str(), "in-directed");
		in_f = fopen(tmp_in_graph_file.c_str(), "w");
		BOOST_VERIFY(fseek(in_f, sizeof(graph_header), SEEK_SET) == 0);
		tmp_out_graph_file = tempnam(work_dir.c_str(), "out-directed");
		out_f = fopen(tmp_out_graph_file.c_str(), "w");
	}

	~disk_directed_graph() {
		if (in_f) {
			fclose(in_f);
			in_f = NULL;
			unlink(tmp_in_graph_file.c_str());
		}
		if (out_f) {
			fclose(out_f);
			out_f = NULL;
			unlink(tmp_out_graph_file.c_str());
		}
	}

	virtual bool is_directed() const {
		return true;
	}

	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const {
	}

	virtual void add_vertex(const in_mem_vertex &v) {
		disk_graph<edge_data_type>::add_vertex(v);

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
	}

	void add_vertices(const subgraph &subg) {
		const directed_subgraph &d_subg = (const directed_subgraph &) subg;
		for (size_t i = 0; i < d_subg.get_num_vertices(); i++)
			disk_graph<edge_data_type>::add_vertex(d_subg.get_vertex_info(i));
		BOOST_VERIFY(fwrite(d_subg.get_in_buf(), d_subg.get_in_size(), 1,
					in_f) == 1);
		BOOST_VERIFY(fwrite(d_subg.get_out_buf(), d_subg.get_out_size(), 1,
					out_f) == 1);

		struct timeval curr;
		gettimeofday(&curr, NULL);
		printf("%f: write %ld and %ld bytes for v[%d, %d)\n",
				time_diff(start_time, curr), d_subg.get_in_size(),
				d_subg.get_out_size(), d_subg.get_start_id(), d_subg.get_end_id());
	}

	void copy_file(FILE *from, size_t from_size, FILE *to) {
		const size_t BUF_SIZE = 128 * 1024 * 1024;
		std::unique_ptr<char[]> buf = std::unique_ptr<char[]>(new char[BUF_SIZE]);
		size_t remain_size = from_size;
		size_t read_size = std::min(remain_size, BUF_SIZE);
		while (read_size > 0) {
			size_t ret = fread(buf.get(), read_size, 1, from);
			BOOST_VERIFY(ret == 1);
			ret = fwrite(buf.get(), read_size, 1, to);
			BOOST_VERIFY(ret == 1);
			remain_size -= read_size;
			read_size = std::min(remain_size, BUF_SIZE);
		}
	}

	virtual void finalize_graph_file(const std::string &adj_file) {
		long out_size = ftell(out_f);
		assert(out_size > 0);
		fclose(out_f);

		out_f = fopen(tmp_out_graph_file.c_str(), "r");
		assert(out_f);
		copy_file(out_f, out_size, in_f);
		fclose(out_f);
		out_f = NULL;
		unlink(tmp_out_graph_file.c_str());

		// Write the real graph header.
		graph_header header(get_graph_type(), this->get_num_vertices(),
				this->get_num_edges(),
				this->get_edge_graph()->has_edge_data() ? sizeof(edge_data_type) : 0);
		BOOST_VERIFY(fseek(in_f, 0, SEEK_SET) == 0);
		BOOST_VERIFY(fwrite(&header, sizeof(header), 1, in_f) == 1);
		fclose(in_f);
		in_f = NULL;
		BOOST_VERIFY(rename(tmp_in_graph_file.c_str(), adj_file.c_str()) == 0);
	}

	virtual graph_type get_graph_type() const {
		return graph_type::DIRECTED;
	}
};

template<class edge_data_type>
class disk_undirected_graph: public disk_graph<edge_data_type>
{
	FILE *f;
	embedded_array<char> buf;
	std::string tmp_graph_file;
public:
	disk_undirected_graph(const edge_graph<edge_data_type> *g): disk_graph<edge_data_type>(
			g, new undirected_in_mem_vertex_index()) {
		tmp_graph_file = tempnam(work_dir.c_str(), "undirected");
		f = fopen(tmp_graph_file.c_str(), "w");
		BOOST_VERIFY(fseek(f, sizeof(graph_header), SEEK_SET) == 0);
	}

	~disk_undirected_graph() {
		if (f) {
			fclose(f);
			f = NULL;
			unlink(tmp_graph_file.c_str());
		}
	}

	virtual bool is_directed() const {
		return false;
	}

	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const {
	}

	virtual size_t get_num_edges() const {
		return disk_graph<edge_data_type>::get_num_edges() / 2;
	}

	virtual void add_vertex(const in_mem_vertex &v) {
		disk_graph<edge_data_type>::add_vertex(v);
		assert(f);
		int size = v.get_serialize_size(IN_EDGE);
		buf.resize(size);
		ext_mem_undirected_vertex::serialize(v, buf.data(), size, IN_EDGE);
		BOOST_VERIFY(fwrite(buf.data(), size, 1, f) == 1);
	}

	virtual void finalize_graph_file(const std::string &adj_file) {
		// Write the real graph header.
		graph_header header(get_graph_type(), this->get_num_vertices(),
				this->get_num_edges(),
				this->get_edge_graph()->has_edge_data() ? sizeof(edge_data_type) : 0);
		BOOST_VERIFY(fseek(f, 0, SEEK_SET) == 0);
		BOOST_VERIFY(fwrite(&header, sizeof(header), 1, f) == 1);
		fclose(f);
		f = NULL;
		BOOST_VERIFY(rename(tmp_graph_file.c_str(), adj_file.c_str()) == 0);
	}

	void add_vertices(const subgraph &subg) {
		const undirected_subgraph &u_subg = (const undirected_subgraph &) subg;
		for (size_t i = 0; i < u_subg.get_num_vertices(); i++)
			disk_graph<edge_data_type>::add_vertex(u_subg.get_vertex_info(i));
		BOOST_VERIFY(fwrite(u_subg.get_buf(), u_subg.get_size(), 1, f) == 1);

		struct timeval curr;
		gettimeofday(&curr, NULL);
		printf("%f: write %ld bytes for v[%d, %d)\n",
				time_diff(start_time, curr), u_subg.get_size(),
				u_subg.get_start_id(), u_subg.get_end_id());
	}

	virtual graph_type get_graph_type() const {
		return graph_type::UNDIRECTED;
	}
};

template<class edge_data_type = empty_data>
class undirected_edge_graph: public edge_graph<edge_data_type>
{
	typedef std::vector<edge<edge_data_type> > edge_list_t;
	typedef typename stxxl_edge_vector<edge_data_type>::const_iterator edge_const_iterator;
	typedef stxxl::stream::vector_iterator2stream<edge_const_iterator> edge_stream_t;

	std::vector<std::shared_ptr<stxxl_edge_vector<edge_data_type> > > edge_lists;

	off_t add_edges(const stxxl_edge_vector<edge_data_type> &edges, off_t idx,
			vertex_id_t id, std::vector<edge<edge_data_type> > &v_edges) const;

	void read_edges(edge_stream_t &, vertex_id_t until_id, edge_list_t &v_edges) const;

	vertex_id_t get_max_vertex_id() const {
		vertex_id_t max_id = 0;
		for (size_t i = 0; i < edge_lists.size(); i++)
			max_id = std::max(edge_lists[i]->back().get_from(), max_id);
		return max_id;
	}
public:
	/**
	 * num_edges tells the edge graph that there will be num_edges
	 * edges added to the graph.
	 */
	undirected_edge_graph(
			std::vector<std::shared_ptr<stxxl_edge_vector<edge_data_type> > > &edge_lists,
			bool has_data): edge_graph<edge_data_type>(has_data) {
		this->edge_lists = edge_lists;
	}

	void sort_edges() {
		comp_edge<edge_data_type> edge_comparator;
		for (size_t i = 0; i < edge_lists.size(); i++)
			stxxl::sort(edge_lists[i]->begin(), edge_lists[i]->end(),
					edge_comparator, SORT_BUF_SIZE);
	}

	graph *create_disk_graph() const {
		return new disk_undirected_graph<edge_data_type>(this);
	}

	size_t get_num_edges() const {
		size_t num_edges = 0;
		for (size_t i = 0; i < edge_lists.size(); i++)
			num_edges += edge_lists[i]->size();
		return num_edges / 2;
	}

#if 0
	void check_vertices(
			const std::vector<ext_mem_undirected_vertex *> &vertices,
			bool in_part) const;
#endif
	void construct_graph(graph *g) const;
};

/**
 * This represents a directed graph in the form of edge list.
 * It maintains a sorted list of out-edges (sorted on the from vertices)
 * and a sorted list of in-edges (sorted on the to vertices).
 */
template<class edge_data_type = empty_data>
class directed_edge_graph: public edge_graph<edge_data_type>
{
	typedef std::vector<edge<edge_data_type> > edge_list_t;
	typedef typename stxxl_edge_vector<edge_data_type>::const_iterator edge_const_iterator;
	typedef stxxl::stream::vector_iterator2stream<edge_const_iterator> edge_stream_t;

	std::vector<std::shared_ptr<stxxl_edge_vector<edge_data_type> > > in_edge_lists;
	std::vector<std::shared_ptr<stxxl_edge_vector<edge_data_type> > > out_edge_lists;

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
public:
	/**
	 * num_edges tells the edge graph that there will be num_edges
	 * edges added to the graph.
	 */
	directed_edge_graph(
			std::vector<std::shared_ptr<stxxl_edge_vector<edge_data_type> > > &edge_lists,
			bool has_data): edge_graph<edge_data_type>(has_data) {
		this->in_edge_lists = edge_lists;
		this->out_edge_lists.resize(edge_lists.size());
		for (size_t i = 0; i < edge_lists.size(); i++)
			this->out_edge_lists[i]
				= std::shared_ptr<stxxl_edge_vector<edge_data_type> >(
						new stxxl_edge_vector<edge_data_type>(*edge_lists[i]));
	}

	void sort_edges() {
		comp_edge<edge_data_type> edge_comparator;
		comp_in_edge<edge_data_type> in_edge_comparator;
		for (size_t i = 0; i < in_edge_lists.size(); i++) {
			stxxl::sort(out_edge_lists[i]->begin(), out_edge_lists[i]->end(),
					edge_comparator, SORT_BUF_SIZE);
			stxxl::sort(in_edge_lists[i]->begin(), in_edge_lists[i]->end(),
					in_edge_comparator, SORT_BUF_SIZE);
			printf("sort edge list %ld\n", i);
		}
	}

#if 0
	void check_vertices(
			const std::vector<ext_mem_undirected_vertex *> &vertices,
			bool in_part) const;
#endif
	void construct_graph(graph *g) const;

	graph *create_disk_graph() const {
		return new disk_directed_graph<edge_data_type>(this);
	}

	size_t get_num_edges() const {
		size_t num_edges = 0;
		for (size_t i = 0; i < in_edge_lists.size(); i++)
			num_edges += in_edge_lists[i]->size();
		return num_edges;
	}
};

template<class edge_data_type>
off_t undirected_edge_graph<edge_data_type>::add_edges(
		const stxxl_edge_vector<edge_data_type> &edges, off_t idx,
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

#if 0
template<class edge_data_type>
void undirected_edge_graph<edge_data_type>::construct_graph(graph *g) const
{
	std::vector<off_t> idxs(edge_lists.size());
	vertex_id_t max_id = get_max_vertex_id();
	std::vector<edge<edge_data_type> > v_edges;
	comp_edge<edge_data_type> edge_comparator;
	for (vertex_id_t id = 0; id <= max_id; id++) {
		v_edges.clear();
		for (size_t i = 0; i < edge_lists.size(); i++) {
			idxs[i] = add_edges(*edge_lists[i], idxs[i], id, v_edges);
		}
		assert(std::is_sorted(v_edges.begin(), v_edges.end(), edge_comparator));
#if 0
		std::sort(v_edges.begin(), v_edges.end(), edge_comparator);
#endif
		in_mem_undirected_vertex<edge_data_type> v(id,
				edge_graph<edge_data_type>::has_edge_data());
		edge<edge_data_type> prev_e;
		BOOST_FOREACH(edge<edge_data_type> e, v_edges) {
			if (prev_e.get_from() == e.get_from()
					&& prev_e.get_to() == e.get_to()) {
				fprintf(stderr, "duplicate edge %u-%u\n", e.get_from(),
						e.get_to());
				continue;
			}

			v.add_edge(e);
			prev_e = e;
		}
		g->add_vertex(v);
	}
}
#endif

#if 0
template<class edge_data_type>
void undirected_edge_graph<edge_data_type>::check_vertices(
		const std::vector<ext_mem_undirected_vertex *> &vertices, bool) const
{
	printf("There are %ld edges in the graph\n", get_num_edges());
	assert(!vertices.empty());
	typename stxxl_edge_vector<edge_data_type>::const_iterator it
		= std::lower_bound(edges.begin(), edges.end(),
				edge<edge_data_type>(0, vertices[0]->get_id()),
				comp_edge<edge_data_type>());

	for (size_t i = 0; i < vertices.size(); i++) {
		ext_mem_undirected_vertex *v = vertices[i];
		for (size_t j = 0; j < v->get_num_edges(); j++, it++) {
			assert(it != edges.end());
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
	typename stxxl_edge_vector<edge_data_type>::const_iterator in_it;
	typename stxxl_edge_vector<edge_data_type>::const_iterator out_it;

	if (in_part)
		in_it = std::lower_bound(in_edges.begin(), in_edges.end(),
				edge<edge_data_type>(0, vertices[0]->get_id()),
				comp_in_edge<edge_data_type>());
	else
		out_it = std::lower_bound(out_edges.begin(), out_edges.end(),
				edge<edge_data_type>(vertices[0]->get_id(), 0),
				comp_edge<edge_data_type>());

	for (size_t i = 0; i < vertices.size(); i++) {
		// Check in-edges
		if (in_part) {
			ext_mem_undirected_vertex *v = vertices[i];
			for (size_t j = 0; j < v->get_num_edges(); j++, in_it++) {
				assert(in_it != in_edges.end());
				edge<edge_data_type> e = *in_it;
				assert(v->get_neighbor(j) == e.get_from());
				assert(v->get_id() == e.get_to());
				if (v->has_edge_data())
					assert(v->get_edge_data<edge_data_type>(j) == e.get_data());
			}
		}
		else {
			// Check out-edges
			ext_mem_undirected_vertex *v = vertices[i];
			for (size_t j = 0; j < v->get_num_edges(); j++, out_it++) {
				assert(out_it != out_edges.end());
				edge<edge_data_type> e = *out_it;
				assert(v->get_id() == e.get_from());
				assert(v->get_neighbor(j) == e.get_to());
				if (v->has_edge_data())
					assert(v->get_edge_data<edge_data_type>(j) == e.get_data());
			}
		}
	}
}
#endif

vsize_t BUF_SIZE = 1024 * 1024 * 1024;

size_t cal_vertex_size(const std::vector<ext_mem_vertex_info> &infos)
{
	assert(!infos.empty());
	return infos.back().get_off() + infos.back().get_size()
		- infos.front().get_off();
}

std::unique_ptr<char[]> read_vertices(FILE *f,
		const std::vector<ext_mem_vertex_info> &infos,
		std::vector<ext_mem_undirected_vertex *> &vertices)
{
	size_t size = cal_vertex_size(infos);
	std::unique_ptr<char[]> buf = std::unique_ptr<char[]>(new char[size]);
	off_t off_begin = infos.front().get_off();
	BOOST_VERIFY(fseek(f, off_begin, SEEK_SET) == 0);
	BOOST_VERIFY(fread(buf.get(), size, 1, f) == 1);
	BOOST_FOREACH(ext_mem_vertex_info info, infos) {
		off_t rel_off = info.get_off() - off_begin;
		vertices.push_back((ext_mem_undirected_vertex *) (buf.get() + rel_off));
	}
	return buf;
}

#if 0
template<class VertexIndexType, class GetInfoFunc, class edge_data_type>
size_t check_all_vertices(FILE *f, const VertexIndexType &idx, GetInfoFunc func,
		const edge_graph<edge_data_type> &edge_g, bool in_part)
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
		std::unique_ptr<char[]> buf = read_vertices(f, infos, vertices);
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

template<class edge_data_type>
void disk_undirected_graph<edge_data_type>::check_ext_graph(
		const std::string &index_file, const std::string &adj_file) const
{
	printf("check the graph in the external memory\n");
	default_vertex_index::ptr idx = default_vertex_index::cast(
			vertex_index::load(index_file));
	FILE *f = fopen(adj_file.c_str(), "r");
	assert(f);
	struct get_undirected_info_func {
		ext_mem_vertex_info operator()(const default_vertex_index &idx,
				vertex_id_t id) {
			return idx.get_vertex_info(id);
		}
	};
	size_t num_vertices = check_all_vertices(f, *idx, get_undirected_info_func(),
			*disk_graph<edge_data_type>::get_edge_graph(), true);
	fclose(f);
	printf("%ld vertices are checked\n", num_vertices);
}

template<class edge_data_type>
void disk_directed_graph<edge_data_type>::check_ext_graph(
		const std::string &index_file, const std::string &adj_file) const
{
	printf("check the graph in the external memory\n");
	directed_vertex_index::ptr idx = directed_vertex_index::cast(
			vertex_index::load(index_file));
	FILE *f = fopen(adj_file.c_str(), "r");
	assert(f);

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
	size_t num_vertices = check_all_vertices(f, *idx, get_in_part_info_func(),
			*disk_graph<edge_data_type>::get_edge_graph(), true);
	size_t num_vertices1 = check_all_vertices(f, *idx, get_out_part_info_func(),
			*disk_graph<edge_data_type>::get_edge_graph(), false);
	assert(num_vertices == num_vertices1);

	fclose(f);
	printf("%ld vertices are checked\n", num_vertices);
}
#endif

template<class edge_data_type>
void disk_graph<edge_data_type>::dump(const std::string &index_file,
			const std::string &graph_file)
{
	struct timeval start, end;
	gettimeofday(&start, NULL);
	assert(g);
	assert(!file_exist(graph_file));

	// Write the adjacency lists to the graph file.
	g->construct_graph(this);
	assert(g->get_num_edges() == get_num_edges());
	finalize_graph_file(graph_file);
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to dump the graph\n",
			time_diff(start, end));

	start = end;
	assert(this->get_num_edges() == g->get_num_edges());
	graph_header header(get_graph_type(), this->get_num_vertices(),
			this->get_num_edges(),
			this->get_edge_graph()->has_edge_data() ? sizeof(edge_data_type) : 0);
	index->dump(index_file, header);
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to dump the index\n",
			time_diff(start, end));
}

template<class edge_data_type = empty_data>
edge_graph<edge_data_type> *par_load_edge_list_text(
		const std::vector<std::string> &files, bool has_edge_data,
		bool directed);

class graph_file_io
{
	FILE *f;
	ssize_t file_size;
public:
	graph_file_io(const std::string file) {
		f = fopen(file.c_str(), "r");
		if (f == NULL)
			ABORT_MSG(boost::format("fail to open %1%: %2%")
					% file % strerror(errno));
		native_file local_f(file);
		file_size = local_f.get_size();
	}

	~graph_file_io() {
		if (f)
			fclose(f);
	}

	std::unique_ptr<char[]> read_edge_list_text(const size_t wanted_bytes,
			size_t &read_bytes);

	size_t get_num_remaining_bytes() const {
		off_t curr_off = ftell(f);
		return file_size - curr_off;
	}
};

/**
 * It read a text of an edge list roughly the size of the wanted bytes.
 * The returned text may be a little more than the wanted bytes, but
 * it's guaranteed that all lines are complete.
 * The returned string ends with '\0'.
 */
std::unique_ptr<char[]> graph_file_io::read_edge_list_text(
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

static std::unique_ptr<char[]> read_file(const std::string &file_name,
		size_t &size)
{
	native_file local_f(file_name);
	size = local_f.get_size();
	FILE *f = fopen(file_name.c_str(), "r");
	assert(f);
	char *buf = new char[size];
	BOOST_VERIFY(fread(buf, size, 1, f) == 1);
	return std::unique_ptr<char[]>(buf);
}

static std::unique_ptr<char[]> read_gz_file(const std::string &file_name,
		size_t &size)
{
	printf("read gz file: %s\n", file_name.c_str());
	const size_t BUF_SIZE = 1024 * 1024 * 16;
	std::vector<std::shared_ptr<char> > bufs;
	gzFile f = gzopen(file_name.c_str(), "rb");
	size_t out_size = 0;
	while (!gzeof(f)) {
		char *buf = new char[BUF_SIZE];
		bufs.push_back(std::shared_ptr<char>(buf));
		int ret = gzread(f, buf, BUF_SIZE);
		assert(ret > 0);
		out_size += ret;
	}
	printf("get %ld bytes from %s\n", out_size, file_name.c_str());

	size = out_size;
	char *out_buf = new char[out_size];
	std::unique_ptr<char[]> ret_buf(out_buf);
	for (size_t i = 0; i < bufs.size(); i++) {
		char *buf = bufs[i].get();
		assert(out_size > 0);
		size_t buf_size = std::min(BUF_SIZE, out_size);
		memcpy(out_buf, buf, buf_size);
		out_buf += buf_size;
		out_size -= buf_size;
	}
	assert(out_size == 0);
	gzclose(f);
	return ret_buf;
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
		*line_end = 0;
		if (*(line_end - 1) == '\r')
			*(line_end - 1) = 0;
		edge<edge_data_type> e;
		int num = parse_edge_list_line(line, e);
		if (num > 0)
			edges.push_back(e);
		num_edges += num;
		line = line_end + 1;
		assert(line - line_end <= (ssize_t) size);
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
		stxxl_edge_vector<edge_data_type> *local_edge_buf
			= (stxxl_edge_vector<edge_data_type> *) thread::get_curr_thread()->get_user_data();
		local_edge_buf->append(edges.cbegin(), edges.cend());

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

	void run() {
		size_t size =  0;
		std::unique_ptr<char[]> data;
		if (decompress)
			data = read_gz_file(file_name, size);
		else
			data = read_file(file_name, size);

		std::vector<edge<edge_data_type> > edges;
		parse_edge_list_text(data.get(), size, edges);
		stxxl_edge_vector<edge_data_type> *local_edge_buf
			= (stxxl_edge_vector<edge_data_type> *) thread::get_curr_thread()->get_user_data();
		local_edge_buf->append(edges.cbegin(), edges.cend());
		printf("There are %lld edges in thread %d\n", local_edge_buf->size(),
				thread::get_curr_thread()->get_id());
	}
};

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

template<class edge_data_type>
class write_graph_thread: public thread
{
	typedef std::shared_ptr<subgraph> subgraph_ptr;
	struct subgraph_comp {
		bool operator()(const subgraph_ptr &g1, const subgraph_ptr &g2) {
			return g1->get_start_id() > g2->get_start_id();
		}
	};

	std::vector<subgraph_ptr> added_subgraphs;
	std::priority_queue<subgraph_ptr, std::vector<subgraph_ptr>, subgraph_comp> subgraphs;
	pthread_spinlock_t lock;
	disk_graph<edge_data_type> &g;
	vertex_id_t curr_id;
	vertex_id_t max_id;
public:
	write_graph_thread(disk_graph<edge_data_type> &_g,
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

template<class edge_data_type>
void write_graph_thread<edge_data_type>::run()
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
	printf("write %d vertices\n", curr_id);
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
	write_graph_thread<edge_data_type> &write_thread;
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
	construct_directed_vertex_task(
			write_graph_thread<edge_data_type> &_write_thread,
			bool has_edge_data, vertex_id_t start_id, vertex_id_t end_id,
			std::shared_ptr<edge_list_t> in_edges,
			std::shared_ptr<edge_list_t> out_edges): write_thread(_write_thread) {
		this->in_edges = in_edges;
		this->out_edges = out_edges;
		this->start_id = start_id;
		this->end_id = end_id;
		this->has_edge_data = has_edge_data;
		struct timeval curr;
		gettimeofday(&curr, NULL);
		printf("%f: create a task for [%d, %d), %ld in-edges and %ld out-edges\n",
				time_diff(start_time, curr), start_id, end_id, in_edges->size(),
				out_edges->size());
	}

	~construct_directed_vertex_task() {
		struct timeval curr;
		gettimeofday(&curr, NULL);
		printf("%f: task completes for [%d, %d)\n", time_diff(start_time, curr),
				start_id, end_id);
	}

	void run() {
		comp_edge<edge_data_type> edge_comparator;
		comp_in_edge<edge_data_type> in_edge_comparator;
		std::sort(in_edges->begin(), in_edges->end(), in_edge_comparator);
		std::sort(out_edges->begin(), out_edges->end(), edge_comparator);

		std::shared_ptr<directed_subgraph> subg
			= std::shared_ptr<directed_subgraph>(new directed_subgraph());
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
	write_graph_thread<edge_data_type> &write_thread;
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
	construct_undirected_vertex_task(
			write_graph_thread<edge_data_type> &_write_thread,
			bool has_edge_data, vertex_id_t start_id, vertex_id_t end_id,
			std::shared_ptr<edge_list_t> edges): write_thread(_write_thread) {
		this->edges = edges;
		this->start_id = start_id;
		this->end_id = end_id;
		this->has_edge_data = has_edge_data;
		struct timeval curr;
		gettimeofday(&curr, NULL);
		printf("%f: create a task for [%d, %d), %ld edges\n",
				time_diff(start_time, curr), start_id, end_id, edges->size());
	}

	~construct_undirected_vertex_task() {
		struct timeval curr;
		gettimeofday(&curr, NULL);
		printf("%f: task completes for [%d, %d)\n", time_diff(start_time, curr),
				start_id, end_id);
	}

	void run() {
		comp_edge<edge_data_type> edge_comparator;
		std::sort(edges->begin(), edges->end(), edge_comparator);

		std::shared_ptr<undirected_subgraph> subg
			= std::shared_ptr<undirected_subgraph>(new undirected_subgraph());
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
void undirected_edge_graph<edge_data_type>::construct_graph(graph *g) const
{
	std::vector<edge_stream_t> its;
	for (size_t i = 0; i < edge_lists.size(); i++)
		its.push_back(edge_stream_t(edge_lists[i]->cbegin(),
				edge_lists[i]->cend()));
	vertex_id_t max_id = get_max_vertex_id();

	gettimeofday(&start_time, NULL);
	std::vector<task_thread *> threads(num_threads);
	for (int i = 0; i < num_threads; i++) {
		task_thread *t = new task_thread(std::string(
					"graph-task-thread") + itoa(i), -1);
		t->start();
		threads[i] = t;
	}
	write_graph_thread<edge_data_type> *write_thread
		= new write_graph_thread<edge_data_type>(
				(disk_undirected_graph<edge_data_type> &) *g, max_id);
	write_thread->start();

	printf("start to construct the graph. max id: %d\n", max_id);

	int thread_no = 0;
	for (vertex_id_t id = 0; id <= max_id; ) {
		std::shared_ptr<edge_list_t> v_edges
			= std::shared_ptr<edge_list_t>(new edge_list_t());
		vertex_id_t end_id = std::min(id + VERTEX_TASK_SIZE, max_id + 1);
		for (size_t i = 0; i < edge_lists.size(); i++)
			read_edges(its[i], end_id, *v_edges);

		construct_undirected_vertex_task<edge_data_type> *task
			= new construct_undirected_vertex_task<edge_data_type>(*write_thread,
					edge_graph<edge_data_type>::has_edge_data(),
					id, end_id, v_edges);
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
}

template<class edge_data_type>
void directed_edge_graph<edge_data_type>::construct_graph(graph *g) const
{
	assert(in_edge_lists.size() == out_edge_lists.size());
	for (size_t i = 0; i < in_edge_lists.size(); i++)
		assert(in_edge_lists[i]->size() == out_edge_lists[i]->size());

	std::vector<edge_stream_t> out_its;
	std::vector<edge_stream_t> in_its;
	for (size_t i = 0; i < out_edge_lists.size(); i++) {
		out_its.push_back(edge_stream_t(out_edge_lists[i]->cbegin(),
				out_edge_lists[i]->cend()));
		in_its.push_back(edge_stream_t(in_edge_lists[i]->cbegin(),
				in_edge_lists[i]->cend()));
	}
	vertex_id_t max_id = get_max_vertex_id();

	gettimeofday(&start_time, NULL);
	std::vector<task_thread *> threads(num_threads);
	for (int i = 0; i < num_threads; i++) {
		task_thread *t = new task_thread(std::string(
					"graph-task-thread") + itoa(i), -1);
		t->start();
		threads[i] = t;
	}
	write_graph_thread<edge_data_type> *write_thread
		= new write_graph_thread<edge_data_type>(
				(disk_directed_graph<edge_data_type> &) *g, max_id);
	write_thread->start();

	printf("start to construct the graph. max id: %d\n", max_id);

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
					edge_graph<edge_data_type>::has_edge_data(),
					id, end_id, v_in_edges, v_out_edges);
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
}

/**
 * This function loads edge lists from a tex file, parses them in parallel,
 * and convert the graph into the form of adjacency lists.
 */
template<class edge_data_type>
edge_graph<edge_data_type> *par_load_edge_list_text(
		const std::vector<std::string> &files, bool has_edge_data, bool directed)
{
	struct timeval start, end;
	gettimeofday(&start, NULL);
	std::vector<task_thread *> threads(num_threads);
	for (int i = 0; i < num_threads; i++) {
		task_thread *t = new task_thread(std::string(
					"graph-task-thread") + itoa(i), -1);
		t->set_user_data(new stxxl_edge_vector<edge_data_type>());
		t->start();
		threads[i] = t;
	}
	int thread_no = 0;
	if (files.size() == 1) {
		const std::string file = files[0];
		printf("start to read the edge list from %s\n", file.c_str());
		graph_file_io io(file);
		while (io.get_num_remaining_bytes() > 0) {
			size_t size = 0;
			thread_task *task = new text_edge_task<edge_data_type>(
					io.read_edge_list_text(EDGE_LIST_BLOCK_SIZE, size),
					size, directed);
			threads[thread_no % num_threads]->add_task(task);
			thread_no++;
		}
	}
	else {
		for (size_t i = 0; i < files.size(); i++) {
			printf("read file %s\n", files[i].c_str());
			thread_task *task = new text_edge_file_task<edge_data_type>(files[i]);
			threads[thread_no % num_threads]->add_task(task);
			thread_no++;
		}
	}
	for (int i = 0; i < num_threads; i++)
		threads[i]->wait4complete();
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to construct edge list\n",
			time_diff(start, end));

	size_t mem_size = 0;
	size_t num_edges = 0;
	std::vector<std::shared_ptr<stxxl_edge_vector<edge_data_type> > > edge_lists(
			num_threads);
	for (int i = 0; i < num_threads; i++) {
		stxxl_edge_vector<edge_data_type> *local_edges
			= (stxxl_edge_vector<edge_data_type> *) threads[i]->get_user_data();
		num_edges += local_edges->size();
		mem_size += local_edges->capacity() * sizeof(edge<edge_data_type>);
		edge_lists[i] = std::shared_ptr<stxxl_edge_vector<edge_data_type> >(
				local_edges);
	}
	printf("There are %ld edges and use %ld bytes\n", num_edges, mem_size);

	edge_graph<edge_data_type> *edge_g;
	if (directed)
		edge_g = new directed_edge_graph<edge_data_type>(edge_lists,
				has_edge_data);
	else
		edge_g = new undirected_edge_graph<edge_data_type>(edge_lists,
				has_edge_data);

	start = end;
	printf("There are %ld edges in the edge graph\n", edge_g->get_num_edges());

	for (int i = 0; i < num_threads; i++) {
		threads[i]->stop();
		threads[i]->join();
		delete threads[i];
	}

	return edge_g;
}

template<class edge_data_type = empty_data>
graph *construct_graph(
		const std::vector<std::string> &edge_list_files,
		bool has_edge_data, bool directed)
{
	printf("beofre load edge list\n");
	struct timeval start, end;
	edge_graph<edge_data_type> *edge_g
		= par_load_edge_list_text<edge_data_type>(edge_list_files,
				has_edge_data, directed);

	printf("before sorting edges\n");
	gettimeofday(&start, NULL);
	edge_g->sort_edges();
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to sort edge list\n", time_diff(start, end));

	return edge_g->create_disk_graph();
}

void print_usage()
{
	fprintf(stderr, "convert an edge list to adjacency lists\n");
	fprintf(stderr,
			"el2al [options] adj_list_file index_file edge_list_files (or directories)\n");
	fprintf(stderr, "-u: undirected graph\n");
	fprintf(stderr, "-v: verify the created adjacency list\n");
	fprintf(stderr, "-t type: the type of edge data. Supported type: ");
	for (int i = 0; i < type_map_size; i++) {
		fprintf(stderr, "%s, ", edge_type_map[i].str.c_str());
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "-m: merge multiple edge lists into a single graph. \n");
	fprintf(stderr, "-w: write the graph to a file\n");
	fprintf(stderr, "-T: the number of threads to process in parallel\n");
	fprintf(stderr, "-W dir: the working directory\n");
	fprintf(stderr, "-D: decompress data\n");
}

graph *construct_graph(const std::vector<std::string> &edge_list_files,
		int edge_attr_type, bool directed)
{
	graph *g = NULL;
	switch(edge_attr_type) {
		case EDGE_COUNT:
			g = construct_graph<edge_count>(
					edge_list_files, true, directed);
			break;
		case EDGE_TIMESTAMP:
			g = construct_graph<ts_edge_data>(
					edge_list_files, true, directed);
			break;
		default:
			g = construct_graph<>(edge_list_files, false,
					directed);
	}
	return g;
}

int main(int argc, char *argv[])
{
	int opt;
	bool directed = true;
	int num_opts = 0;
	char *type_str = NULL;
	bool merge_graph = false;
	bool write_graph = false;
	while ((opt = getopt(argc, argv, "uvt:mwT:W:D")) != -1) {
		num_opts++;
		switch (opt) {
			case 'u':
				directed = false;
				break;
			case 'v':
				check_graph = true;
				break;
			case 't':
				type_str = optarg;
				num_opts++;
				break;
			case 'm':
				merge_graph = true;
				break;
			case 'w':
				write_graph = true;
				break;
			case 'T':
				num_threads = atoi(optarg);
				num_opts++;
				break;
			case 'W':
				work_dir = optarg;
				num_opts++;
				break;
			case 'D':
				decompress = true;
				break;
			default:
				print_usage();
		}
	}
	argv += 1 + num_opts;
	argc -= 1 + num_opts;
	if (argc < 3) {
		print_usage();
		exit(-1);
	}

	int edge_attr_type = DEFAULT_TYPE;
	if (type_str) {
		edge_attr_type = conv_edge_type_str2int(type_str);
	}
	printf("work dir: %s\n", work_dir.c_str());

	std::string adjacency_list_file = argv[0];
	adjacency_list_file += std::string("-v") + itoa(CURR_VERSION);
	std::string index_file = argv[1];
	index_file += std::string("-v") + itoa(CURR_VERSION);
	std::vector<std::string> edge_list_files;
	for (int i = 2; i < argc; i++) {
		native_dir dir(argv[i]);
		if (dir.is_dir()) {
			std::vector<std::string> files;
			std::string dir_name = argv[i];
			dir.read_all_files(files);
			for (size_t i = 0; i < files.size(); i++)
				edge_list_files.push_back(dir_name + "/" + files[i]);
		}
		else
			edge_list_files.push_back(argv[i]);
	}

	for (size_t i = 0; i < edge_list_files.size(); i++)
		printf("edge list file: %s\n", edge_list_files[i].c_str());

	if (merge_graph) {
		graph *g = construct_graph(edge_list_files, edge_attr_type, directed);
		// Write the constructed individual graph to a file.
		if (write_graph)
			g->dump(index_file, adjacency_list_file);
		printf("There are %ld vertices, %ld non-empty vertices and %ld edges\n",
				g->get_num_vertices(), g->get_num_non_empty_vertices(),
				g->get_num_edges());
#if 0
		if (check_graph)
			g->check_ext_graph(index_file, adjacency_list_file);
#endif
		delete g;
	}
	else {
		std::vector<std::string> graph_files;
		std::vector<std::string> index_files;
		if (edge_list_files.size() > 1) {
			for (size_t i = 0; i < edge_list_files.size(); i++) {
				graph_files.push_back(adjacency_list_file + "-" + itoa(i));
				index_files.push_back(index_file + "-" + itoa(i));
			}
		}
		else {
			graph_files.push_back(adjacency_list_file);
			index_files.push_back(index_file);
		}

		for (size_t i = 0; i < edge_list_files.size(); i++) {
			// construct individual graphs.
			std::vector<std::string> files(1);
			files[0] = edge_list_files[i];
			graph *g = construct_graph(files, edge_attr_type, directed);
			// Write the constructed individual graph to a file.
			if (write_graph)
				g->dump(index_files[i], graph_files[i]);
			printf("There are %ld vertices, %ld non-empty vertices and %ld edges\n",
					g->get_num_vertices(), g->get_num_non_empty_vertices(),
					g->get_num_edges());
#if 0
			if (check_graph)
				g->check_ext_graph(index_files[i], graph_files[i]);
#endif
			delete g;
		}
	}
}
