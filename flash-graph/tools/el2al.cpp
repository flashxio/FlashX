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

#include "graph.h"

#include <memory>
#include <algorithm>

#include <boost/foreach.hpp>
#include <stxxl.h>

#include "thread.h"
#include "native_file.h"

#include "edge_type.h"

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
const int EDGE_LIST_BLOCK_SIZE = 1 * 1024 * 1024;
const char *delimiter = "\t";

bool compress = false;
bool simplify = false;
bool print_graph = false;
bool check_graph = false;
std::string work_dir = ".";

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
	virtual edge_graph<edge_count> *compress_edges() const = 0;
	virtual edge_graph<edge_data_type> *simplify_edges() const = 0;
	virtual void construct_graph(graph *g) const = 0;
	virtual graph *create_disk_graph() const = 0;
	virtual void add_edges(stxxl_edge_vector<edge_data_type> &edges) = 0;
	virtual size_t get_num_edges() const = 0;
	virtual void check_vertices(
			const std::vector<ext_mem_undirected_vertex *> &vertices,
			bool in_part) const = 0;

	bool has_edge_data() const {
		return has_data;
	}
};

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
		assert(0);
	}

	void get_all_vertices(std::vector<vertex_id_t> &ids) const {
		assert(0);
	}

	virtual void print() const {
		assert(0);
	}

	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const = 0;
	virtual graph_type get_graph_type() const = 0;
	virtual void finalize_graph_file(const std::string &adj_file) = 0;
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
		int ret = fseek(in_f, sizeof(graph_header), SEEK_SET);
		assert(ret == 0);
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

	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const;

	virtual void add_vertex(const in_mem_vertex &v) {
		disk_graph<edge_data_type>::add_vertex(v);

		assert(in_f);
		assert(out_f);
		int size = v.get_serialize_size(IN_EDGE);
		buf.resize(size);
		ext_mem_undirected_vertex::serialize(v, buf.data(), size, IN_EDGE);
		ssize_t ret = fwrite(buf.data(), size, 1, in_f);
		assert(ret == 1);

		size = v.get_serialize_size(OUT_EDGE);
		buf.resize(size);
		ext_mem_undirected_vertex::serialize(v, buf.data(), size, OUT_EDGE);
		ret = fwrite(buf.data(), size, 1, out_f);
		assert(ret == 1);
	}

	void copy_file(FILE *from, size_t from_size, FILE *to) {
		const size_t BUF_SIZE = 128 * 1024 * 1024;
		std::unique_ptr<char[]> buf = std::unique_ptr<char[]>(new char[BUF_SIZE]);
		size_t remain_size = from_size;
		size_t read_size = std::min(remain_size, BUF_SIZE);
		while (read_size > 0) {
			size_t ret = fread(buf.get(), read_size, 1, from);
			assert(ret == 1);
			ret = fwrite(buf.get(), read_size, 1, to);
			assert(ret == 1);
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
		int seek_ret = fseek(in_f, 0, SEEK_SET);
		assert(seek_ret == 0);
		ssize_t ret = fwrite(&header, sizeof(header), 1, in_f);
		assert(ret == 1);
		fclose(in_f);
		in_f = NULL;
		assert(rename(tmp_in_graph_file.c_str(), adj_file.c_str()) == 0);
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
		int ret = fseek(f, sizeof(graph_header), SEEK_SET);
		assert(ret == 0);
	}

	~disk_undirected_graph() {
		if (f) {
			fclose(f);
			f = NULL;
			unlink(tmp_graph_file.c_str());
		}
	}

	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const;

	virtual size_t get_num_edges() const {
		return disk_graph<edge_data_type>::get_num_edges() / 2;
	}

	virtual void add_vertex(const in_mem_vertex &v) {
		disk_graph<edge_data_type>::add_vertex(v);
		assert(f);
		int size = v.get_serialize_size(IN_EDGE);
		buf.resize(size);
		ext_mem_undirected_vertex::serialize(v, buf.data(), size, IN_EDGE);
		ssize_t ret = fwrite(buf.data(), size, 1, f);
		assert(ret == 1);
	}

	virtual void finalize_graph_file(const std::string &adj_file) {
		// Write the real graph header.
		graph_header header(get_graph_type(), this->get_num_vertices(),
				this->get_num_edges(),
				this->get_edge_graph()->has_edge_data() ? sizeof(edge_data_type) : 0);
		int seek_ret = fseek(f, 0, SEEK_SET);
		assert(seek_ret == 0);
		ssize_t ret = fwrite(&header, sizeof(header), 1, f);
		assert(ret == 1);
		fclose(f);
		f = NULL;
		assert(rename(tmp_graph_file.c_str(), adj_file.c_str()) == 0);
	}

	virtual graph_type get_graph_type() const {
		return graph_type::UNDIRECTED;
	}
};

template<class edge_data_type = empty_data>
class undirected_edge_graph: public edge_graph<edge_data_type>
{
	stxxl_edge_vector<edge_data_type> edges;
public:
	undirected_edge_graph(bool has_data): edge_graph<edge_data_type>(has_data) {
	}

	/**
	 * num_edges tells the edge graph that there will be num_edges
	 * edges added to the graph.
	 */
	undirected_edge_graph(size_t num_edges,
			bool has_data): edge_graph<edge_data_type>(has_data) {
		edges.reserve(num_edges);
	}

	void sort_edges() {
		comp_edge<edge_data_type> edge_comparator;
//#ifdef MEMCHECK
		stxxl::sort(edges.begin(), edges.end(), edge_comparator, 1024 * 1024 * 10);
//#else
//		__gnu_parallel::sort(edges.begin(), edges.end(), edge_comparator);
//#endif
	}

	graph *create_disk_graph() const {
		return new disk_undirected_graph<edge_data_type>(this);
	}

	/**
	 * This is used by compress_edges, so the edge graph has all the edges.
	 */
	void add_edge(const edge<edge_data_type> &e) {
		edges.push_back(e);
	}

	/**
	 * The input edges are read from the edge list file.
	 * Each edge may appear only once, so we need to reverse them for
	 * the other end of the edges.
	 */
	void add_edges(stxxl_edge_vector<edge_data_type> &edges) {
		this->edges.append(edges.begin(), edges.end());
		BOOST_FOREACH(edge<edge_data_type> &e, edges) {
			if (e.has_edge_data())
				this->edges.push_back(edge<edge_data_type>(e.get_to(),
							e.get_from(), e.get_data()));
			else
				this->edges.push_back(edge<edge_data_type>(e.get_to(),
							e.get_from()));
		}
	}

	size_t get_num_edges() const {
		return edges.size() / 2;
	}

	edge_graph<edge_count> *compress_edges() const;
	edge_graph<edge_data_type> *simplify_edges() const;
	void construct_graph(graph *g) const;
	void check_vertices(
			const std::vector<ext_mem_undirected_vertex *> &vertices,
			bool in_part) const;
};

/**
 * This represents a directed graph in the form of edge list.
 * It maintains a sorted list of out-edges (sorted on the from vertices)
 * and a sorted list of in-edges (sorted on the to vertices).
 */
template<class edge_data_type = empty_data>
class directed_edge_graph: public edge_graph<edge_data_type>
{
	stxxl_edge_vector<edge_data_type> in_edges;
	stxxl_edge_vector<edge_data_type> out_edges;
public:
	directed_edge_graph(bool has_data): edge_graph<edge_data_type>(has_data) {
	}

	/**
	 * num_edges tells the edge graph that there will be num_edges
	 * edges added to the graph.
	 */
	directed_edge_graph(size_t num_edges,
			bool has_data): edge_graph<edge_data_type>(has_data) {
		in_edges.reserve(num_edges);
		out_edges.reserve(num_edges);
	}

	void sort_edges() {
		comp_edge<edge_data_type> edge_comparator;
		comp_in_edge<edge_data_type> in_edge_comparator;
//#ifdef MEMCHECK
		stxxl::sort(out_edges.begin(), out_edges.end(), edge_comparator, 1024 * 1024 * 10);
		stxxl::sort(in_edges.begin(), in_edges.end(), in_edge_comparator, 1024 * 1024 * 10);
//#else
//		__gnu_parallel::sort(out_edges.begin(), out_edges.end(), edge_comparator);
//		__gnu_parallel::sort(in_edges.begin(), in_edges.end(), in_edge_comparator);
//#endif
	}

	edge_graph<edge_count> *compress_edges() const;
	edge_graph<edge_data_type> *simplify_edges() const;
	void construct_graph(graph *g) const;

	graph *create_disk_graph() const {
		return new disk_directed_graph<edge_data_type>(this);
	}

	void add_edge(const edge<edge_data_type> &e) {
		in_edges.push_back(e);
		out_edges.push_back(e);
	}

	void add_edges(stxxl_edge_vector<edge_data_type> &edges) {
		in_edges.append(edges.begin(), edges.end());
		out_edges.append(edges.begin(), edges.end());
	}

	size_t get_num_edges() const {
		assert(in_edges.size() == out_edges.size());
		return in_edges.size();
	}

	void check_vertices(
			const std::vector<ext_mem_undirected_vertex *> &vertices,
			bool in_part) const;
};

template<class edge_data_type>
edge_graph<edge_count> *
undirected_edge_graph<edge_data_type>::compress_edges() const
{
	undirected_edge_graph<edge_count> *new_graph
		= new undirected_edge_graph<edge_count>(true);
	printf("before: %ld edges\n", get_num_edges());
	if (!edges.empty()) {
		vertex_id_t from = edges[0].get_from();
		vertex_id_t to = edges[0].get_to();
		int num_duplicates = 1;
		for (size_t i = 1; i < edges.size(); i++) {
			if (edges[i].get_from() == from && edges[i].get_to() == to) {
				num_duplicates++;
			}
			else {
				edge_count c(num_duplicates);
				edge<edge_count> e(from, to, num_duplicates);
				new_graph->add_edge(e);

				num_duplicates = 1;
				from = edges[i].get_from();
				to = edges[i].get_to();
			}
		}
		edge_count c(num_duplicates);
		edge<edge_count> e(from, to, num_duplicates);
		new_graph->add_edge(e);
	}

	printf("after: %ld edges\n", new_graph->get_num_edges());
	return new_graph;
}

template<class edge_data_type>
edge_graph<edge_data_type> *
undirected_edge_graph<edge_data_type>::simplify_edges() const
{
	undirected_edge_graph<edge_data_type> *new_graph
		= new undirected_edge_graph<edge_data_type>(false);
	std::cout << "before: " << edges.size() << " edges\n";
	if (!edges.empty()) {
		new_graph->edges.push_back(edges[0]);
		for (size_t i = 1; i < edges.size(); i++) {
			if (edges[i].get_from() != new_graph->edges.back().get_from()
					|| edges[i].get_to() != new_graph->edges.back().get_to())
				new_graph->edges.push_back(edges[i]);
		}
	}
	std::cout << "after: " << new_graph->edges.size() << " edges\n";
	return new_graph;
}

template<class edge_data_type>
void undirected_edge_graph<edge_data_type>::construct_graph(graph *g) const
{
	vertex_id_t curr = 0;
	in_mem_undirected_vertex<edge_data_type> v(curr,
			edge_graph<edge_data_type>::has_edge_data());
	size_t idx = 0;
	size_t num_edges = edges.size();

	// Add remaining out-edges.
	while (idx < num_edges) {
		while (idx < num_edges
				&& edges[idx].get_from() == curr) {
			v.add_edge(edges[idx++]);
		}
		g->add_vertex(v);
		vertex_id_t prev = curr + 1;
		if (idx < num_edges)
			curr = edges[idx].get_from();
		else
			break;
		// The vertices without edges won't show up in the edge list,
		// but we need to fill the gap in the vertex Id space with empty
		// vertices.
		while (prev < curr) {
			v = in_mem_undirected_vertex<edge_data_type>(prev,
					edge_graph<edge_data_type>::has_edge_data());
			prev++;
			g->add_vertex(v);
		}
		v = in_mem_undirected_vertex<edge_data_type>(curr,
				edge_graph<edge_data_type>::has_edge_data());
	}
}

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
	int ret = fseek(f, off_begin, SEEK_SET);
	assert(ret == 0);
	ret = fread(buf.get(), size, 1, f);
	assert(ret == 1);
	BOOST_FOREACH(ext_mem_vertex_info info, infos) {
		off_t rel_off = info.get_off() - off_begin;
		vertices.push_back((ext_mem_undirected_vertex *) (buf.get() + rel_off));
	}
	return buf;
}

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
		if (f == NULL) {
			perror("fopen");
			assert(0);
		}
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
	int ret = fread(line_buf, read_bytes, 1, f);
	assert(ret == 1);
	line_buf[read_bytes] = 0;

	return std::unique_ptr<char[]>(line_buf);
}

size_t parse_edge_list_line(char *line, edge<ts_edge_data> &e)
{
	int len = strlen(line);
	/*
	 * The format of a line should be
	 * from_vertex to_vertex "time" weight
	 * Fields are separated by tabs.
	 */
	if (line[len - 1] == '\n')
		line[len - 1] = 0;
	if (line[len - 2] == '\r')
		line[len - 2] = 0;
	if (line[0] == '#')
		return -1;
	char *second = strstr(line, delimiter);
	assert(second);
	*second = 0;
	second += strlen(delimiter);

	char *third = strstr(second, delimiter);
	assert(third);
	*third = 0;
	third += strlen(delimiter);
	if (*third == '"')
		third++;

	if (!isnumeric(line) || !isnumeric(second)) {
		printf("%s\t%s\t%s\n", line, second, third);
		return -1;
	}
	long lfrom = atol(line);
	long lto = atol(second);
	assert(lfrom >= 0 && lfrom < MAX_VERTEX_ID);
	assert(lto >= 0 && lto < MAX_VERTEX_ID);
	vertex_id_t from = lfrom;
	vertex_id_t to = lto;
	time_t timestamp = atol(third);
	ts_edge_data data(timestamp);
	e = edge<ts_edge_data>(from, to, data);
	return 1;
}

int parse_edge_list_line(char *line, edge<edge_count> &e)
{
	if (line[0] == '#')
		return 0;
	char *second = strstr(line, delimiter);
	if (second == NULL) {
		fprintf(stderr, "wrong format 1: %s\n", line);
		return -1;
	}
	*second = 0;
	second += strlen(delimiter);
	char *third = strstr(second, delimiter);
	if (third == NULL) {
		fprintf(stderr, "wrong format 2: %s\n", second);
		return -1;
	}
	*third = 0;
	third += strlen(delimiter);
	if (!isnumeric(line) || !isnumeric(second) || !isnumeric(third)) {
		fprintf(stderr, "wrong format 3: %s\t%s\t%s\n", line, second, third);
		return -1;
	}
	long lfrom = atol(line);
	long lto = atol(second);
	assert(lfrom >= 0 && lfrom < MAX_VERTEX_ID);
	assert(lto >= 0 && lto < MAX_VERTEX_ID);
	vertex_id_t from = lfrom;
	vertex_id_t to = lto;
	edge_count c(atol(third));
	e = edge<edge_count>(from, to, c);

	return 1;
}

int parse_edge_list_line(char *line, edge<> &e)
{
	if (line[0] == '#')
		return 0;
	char *second = strstr(line, delimiter);
	if (second == NULL) {
		fprintf(stderr, "wrong format 1: %s\n", line);
		return -1;
	}
	*second = 0;
	second += strlen(delimiter);
	if (!isnumeric(line) || !isnumeric(second)) {
		fprintf(stderr, "wrong format 2: %s\t%s\n", line, second);
		return -1;
	}
	long lfrom = atol(line);
	long lto = atol(second);
	assert(lfrom >= 0 && lfrom < MAX_VERTEX_ID);
	assert(lto >= 0 && lto < MAX_VERTEX_ID);
	vertex_id_t from = lfrom;
	vertex_id_t to = lto;
	e = edge<>(from, to);

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
public:
	text_edge_task(std::unique_ptr<char[]> line_buf, size_t size) {
		this->line_buf = std::move(line_buf);
		this->size = size;
	}

	void run() {
		std::vector<edge<edge_data_type> > edges;
		parse_edge_list_text(line_buf.get(), size, edges);
		stxxl_edge_vector<edge_data_type> *local_edge_buf
			= (stxxl_edge_vector<edge_data_type> *) thread::get_curr_thread()->get_user_data();
		local_edge_buf->append(edges.cbegin(), edges.cend());
	}
};

template<class edge_data_type>
edge_graph<edge_count> *
directed_edge_graph<edge_data_type>::compress_edges() const
{
	directed_edge_graph<edge_count> *new_graph
		= new directed_edge_graph<edge_count>(true);
	printf("before: %ld edges\n", get_num_edges());
	if (!in_edges.empty()) {
		vertex_id_t from = in_edges[0].get_from();
		vertex_id_t to = in_edges[0].get_to();
		int num_duplicates = 1;
		for (size_t i = 1; i < in_edges.size(); i++) {
			if (in_edges[i].get_from() == from && in_edges[i].get_to() == to) {
				num_duplicates++;
			}
			else {
				edge_count c(num_duplicates);
				edge<edge_count> e(from, to, num_duplicates);
				new_graph->add_edge(e);

				num_duplicates = 1;
				from = in_edges[i].get_from();
				to = in_edges[i].get_to();
			}
		}
		edge_count c(num_duplicates);
		edge<edge_count> e(from, to, num_duplicates);
		new_graph->add_edge(e);
	}
	new_graph->sort_edges();

	printf("after: %ld edges\n", new_graph->get_num_edges());
	return new_graph;
}

template<class edge_data_type>
edge_graph<edge_data_type> *
directed_edge_graph<edge_data_type>::simplify_edges() const
{
	directed_edge_graph<edge_data_type> *new_graph
		= new directed_edge_graph<edge_data_type>(false);
	std::cout << "before: " << in_edges.size() << " in-edges and "
		<< out_edges.size() << " out-edges\n";
	if (!in_edges.empty()) {
		new_graph->in_edges.push_back(in_edges[0]);
		for (size_t i = 1; i < in_edges.size(); i++) {
			if (in_edges[i].get_from() != new_graph->in_edges.back().get_from()
					|| in_edges[i].get_to() != new_graph->in_edges.back().get_to())
				new_graph->in_edges.push_back(in_edges[i]);
		}
	}

	if (!out_edges.empty()) {
		new_graph->out_edges.push_back(out_edges[0]);
		for (size_t i = 1; i < out_edges.size(); i++) {
			if (out_edges[i].get_from() != new_graph->out_edges.back().get_from()
					|| out_edges[i].get_to() != new_graph->out_edges.back().get_to())
				new_graph->out_edges.push_back(out_edges[i]);
		}
	}
	std::cout << "after: " << new_graph->in_edges.size() << " in-edges and "
		<< new_graph->out_edges.size() << " out-edges\n";
	return new_graph;
}

template<class edge_data_type>
void directed_edge_graph<edge_data_type>::construct_graph(graph *g) const
{
	vertex_id_t curr = 0;
	in_mem_directed_vertex<edge_data_type> v(curr,
			edge_graph<edge_data_type>::has_edge_data());
	size_t out_idx = 0;
	size_t in_idx = 0;
	size_t num_edges = in_edges.size();
	assert(in_edges.size() == out_edges.size());
	while (out_idx < num_edges && in_idx < num_edges) {
		while (out_idx < num_edges
				&& out_edges[out_idx].get_from() == curr) {
			v.add_out_edge(out_edges[out_idx++]);
		}
		while (in_idx < num_edges
				&& in_edges[in_idx].get_to() == curr) {
			v.add_in_edge(in_edges[in_idx++]);
		}
		g->add_vertex(v);
		vertex_id_t prev = curr + 1;
		if (out_idx < num_edges && in_idx < num_edges)
			curr = min(out_edges[out_idx].get_from(),
					in_edges[in_idx].get_to());
		else if (out_idx < num_edges)
			curr = out_edges[out_idx].get_from();
		else if (in_idx < num_edges)
			curr = in_edges[in_idx].get_to();
		else
			break;
		// The vertices without edges won't show up in the edge list,
		// but we need to fill the gap in the vertex Id space with empty
		// vertices.
		while (prev < curr) {
			v = in_mem_directed_vertex<edge_data_type>(prev,
					edge_graph<edge_data_type>::has_edge_data());
			prev++;
			g->add_vertex(v);
		}
		v = in_mem_directed_vertex<edge_data_type>(curr,
				edge_graph<edge_data_type>::has_edge_data());
	}

	// Add remaining out-edges.
	while (out_idx < num_edges) {
		while (out_idx < num_edges
				&& out_edges[out_idx].get_from() == curr) {
			v.add_out_edge(out_edges[out_idx++]);
		}
		g->add_vertex(v);
		vertex_id_t prev = curr + 1;
		if (out_idx < num_edges)
			curr = out_edges[out_idx].get_from();
		else
			break;
		// The vertices without edges won't show up in the edge list,
		// but we need to fill the gap in the vertex Id space with empty
		// vertices.
		while (prev < curr) {
			v = in_mem_directed_vertex<edge_data_type>(prev,
					edge_graph<edge_data_type>::has_edge_data());
			prev++;
			g->add_vertex(v);
		}
		v = in_mem_directed_vertex<edge_data_type>(curr,
				edge_graph<edge_data_type>::has_edge_data());
	}

	// Add remaining in-edges
	while (in_idx < num_edges) {
		while (in_idx < num_edges
				&& in_edges[in_idx].get_to() == curr) {
			v.add_in_edge(in_edges[in_idx++]);
		}
		g->add_vertex(v);
		vertex_id_t prev = curr + 1;
		if (in_idx < num_edges)
			curr = in_edges[in_idx].get_to();
		else
			break;
		// The vertices without edges won't show up in the edge list,
		// but we need to fill the gap in the vertex Id space with empty
		// vertices.
		while (prev < curr) {
			v = in_mem_directed_vertex<edge_data_type>(prev,
					edge_graph<edge_data_type>::has_edge_data());
			prev++;
			g->add_vertex(v);
		}
		v = in_mem_directed_vertex<edge_data_type>(curr,
				edge_graph<edge_data_type>::has_edge_data());
	}
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
	printf("start to read the edge list\n");
	for (size_t i = 0; i < files.size(); i++) {
		printf("read file %s\n", files[i].c_str());
		graph_file_io io(files[i]);
		while (io.get_num_remaining_bytes() > 0) {
			size_t size = 0;
			thread_task *task = new text_edge_task<edge_data_type>(
					io.read_edge_list_text(EDGE_LIST_BLOCK_SIZE, size),
					size);
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
	for (int i = 0; i < num_threads; i++) {
		stxxl_edge_vector<edge_data_type> *local_edges
			= (stxxl_edge_vector<edge_data_type> *) threads[i]->get_user_data();
		num_edges += local_edges->size();
		mem_size += local_edges->capacity() * sizeof(edge<edge_data_type>);
	}
	printf("There are %ld edges and use %ld bytes\n", num_edges, mem_size);

	edge_graph<edge_data_type> *edge_g;
	if (directed)
		edge_g = new directed_edge_graph<edge_data_type>(num_edges,
				has_edge_data);
	else
		edge_g = new undirected_edge_graph<edge_data_type>(num_edges,
				has_edge_data);
	start = end;
	for (int i = 0; i < num_threads; i++) {
		stxxl_edge_vector<edge_data_type> *local_edges
			= (stxxl_edge_vector<edge_data_type> *) threads[i]->get_user_data();
		edge_g->add_edges(*local_edges);
		delete local_edges;
	}
	gettimeofday(&end, NULL);
	printf("There are %ld edges in the edge graph\n", edge_g->get_num_edges());
	printf("It takes %f seconds to combine edge list\n", time_diff(start, end));

	for (int i = 0; i < num_threads; i++) {
		threads[i]->stop();
		threads[i]->join();
		delete threads[i];
	}

	return edge_g;
}

template<class edge_data_type = empty_data>
graph *construct_graph_compressed(
		const std::vector<std::string> &edge_list_files, bool directed)
{
	struct timeval start, end;
	edge_graph<edge_data_type> *edge_g
		= par_load_edge_list_text<edge_data_type>(edge_list_files,
				true, directed);

	gettimeofday(&start, NULL);
	edge_g->sort_edges();
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to sort edge list\n", time_diff(start, end));

	start = end;
	size_t orig_num_edges = edge_g->get_num_edges();
	edge_graph<edge_count> *new_edge_g = edge_g->compress_edges();
	delete edge_g;
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to compress edge list from %ld to %ld\n",
			time_diff(start, end), orig_num_edges,
			new_edge_g->get_num_edges());

	return new_edge_g->create_disk_graph();
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

	if (simplify) {
		start = end;
		size_t orig_num_edges = edge_g->get_num_edges();
		edge_graph<edge_data_type> *new_edge_g = edge_g->simplify_edges();
		delete edge_g;
		edge_g = new_edge_g;
		gettimeofday(&end, NULL);
		printf("It takes %f seconds to remove duplicated edges from %ld to %ld\n",
				time_diff(start, end), orig_num_edges, edge_g->get_num_edges());
	}

	return edge_g->create_disk_graph();
}

void print_usage()
{
	fprintf(stderr, "convert an edge list to adjacency lists\n");
	fprintf(stderr,
			"el2al [options] adj_list_file index_file edge_list_files (or directories)\n");
	fprintf(stderr, "-u: undirected graph\n");
	fprintf(stderr, "-d delimiter: the delimiter to seperate the input edge list\n");
	fprintf(stderr, "-c: compress the graph (remove duplicated edges)\n");
	fprintf(stderr, "-p: print adjacency list\n");
	fprintf(stderr, "-v: verify the created adjacency list\n");
	fprintf(stderr, "-t type: the type of edge data. Supported type: ");
	for (int i = 0; i < type_map_size; i++) {
		fprintf(stderr, "%s, ", edge_type_map[i].str.c_str());
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "-m: merge multiple edge lists into a single graph. \n");
	fprintf(stderr, "-w: write the graph to a file\n");
	fprintf(stderr, "-s: simplify a graph (remove duplicated edges)\n");
	fprintf(stderr, "-T: the number of threads to process in parallel\n");
	fprintf(stderr, "-W dir: the working directory\n");
}

graph *construct_graph(const std::vector<std::string> &edge_list_files,
		int edge_attr_type, bool directed)
{
	graph *g = NULL;
	switch(edge_attr_type) {
		case EDGE_COUNT:
			if (compress)
				g = construct_graph_compressed<edge_count>(
						edge_list_files, directed);
			else
				g = construct_graph<edge_count>(
						edge_list_files, true, directed);
			break;
		case EDGE_TIMESTAMP:
			if (compress)
				g = construct_graph_compressed<ts_edge_data>(
						edge_list_files, directed);
			else
				g = construct_graph<ts_edge_data>(
						edge_list_files, true, directed);
			break;
		default:
			if (compress)
				g = construct_graph_compressed<>(edge_list_files,
						directed);
			else
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
	while ((opt = getopt(argc, argv, "ud:cpvt:mwsT:W:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'u':
				directed = false;
				break;
			case 'd':
				delimiter = optarg;
				num_opts++;
				break;
			case 'c':
				compress = true;
				break;
			case 'p':
				print_graph = true;
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
			case 's':
				simplify = true;
				break;
			case 'T':
				num_threads = atoi(optarg);
				num_opts++;
				break;
			case 'W':
				work_dir = optarg;
				num_opts++;
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
		if (print_graph)
			g->print();
		if (check_graph)
			g->check_ext_graph(index_file, adjacency_list_file);
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
			if (print_graph)
				g->print();
			if (check_graph)
				g->check_ext_graph(index_files[i], graph_files[i]);
			delete g;
		}
	}
}
