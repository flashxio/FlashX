#ifndef __GRAPH_H__
#define __GRAPH_H__

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

#include <string>
#include <set>

#include "native_file.h"

#include "vertex_index.h"
#include "vertex.h"

class graph
{
public:
	typedef std::shared_ptr<graph> ptr;

	virtual ~graph() {
	}

	virtual bool is_directed() const = 0;
	virtual const in_mem_vertex &get_vertex(vertex_id_t id) const = 0;
	virtual void add_vertex(const in_mem_vertex &v) = 0;
	virtual void get_all_vertices(std::vector<vertex_id_t> &ids) const = 0;
	virtual void dump(const std::string &index_file,
			const std::string &graph_file) = 0;
	virtual void dump_as_edge_list(const std::string &file) const {
		// TODO
		ABORT_MSG("dump_as_edge_list isn't implemented");
	}
	virtual size_t get_num_edges() const = 0;
	virtual size_t get_num_vertices() const = 0;
	virtual bool has_edge_data() const = 0;
	virtual size_t get_num_non_empty_vertices() const = 0;
	virtual void print() const = 0;
	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const = 0;
	// Merge the graph to this graph.
	virtual void merge(graph::ptr g) {
		ABORT_MSG("merge isn't implemented");
	}
};

class in_mem_subgraph: public graph
{
protected:
	size_t num_non_empty;
	bool has_data;
	size_t num_edges;

	in_mem_subgraph(bool has_data) {
		num_edges = 0;
		num_non_empty = 0;
		this->has_data = has_data;
	}
public:
	static graph::ptr create(graph_type type, bool has_data);

	virtual void print() const {
		ABORT_MSG("print isn't implemented");
	}
	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const {
		ABORT_MSG("check_ext_graph isn't implemented");
	}
	virtual void dump(const std::string &index_file,
			const std::string &graph_file) {
		ABORT_MSG("dump isn't implemented");
	}
	virtual void dump_as_edge_list(const std::string &file) const {
		// TODO
		ABORT_MSG("dump_as_edge_list isn't implemented");
	}
	virtual size_t get_num_edges() const {
		return num_edges;
	}
	virtual bool has_edge_data() const {
		return has_data;
	}
	virtual size_t get_num_non_empty_vertices() const {
		return num_non_empty;
	}
	/**
	 * This compresses the subgraph and generates a graph whose vertex IDs
	 * are adjacent to each other.
	 */
	virtual graph::ptr compress() const {
		return graph::ptr();
	}

	// Merge the graph to this graph.
	virtual void merge(graph::ptr g) {
		std::vector<vertex_id_t> vertices;
		g->get_all_vertices(vertices);
		for (size_t i = 0; i < vertices.size(); i++) {
			vertex_id_t id = vertices[i];
			add_vertex(g->get_vertex(id));
		}
	}
};

template<class edge_data_type = empty_data>
class in_mem_undirected_subgraph: public in_mem_subgraph
{
	typedef std::map<vertex_id_t, in_mem_undirected_vertex<edge_data_type> > vmap_t;
	vmap_t vertices;

	in_mem_undirected_subgraph(bool has_data): in_mem_subgraph(has_data) {
	}
public:
	static graph::ptr create(bool has_data) {
		return graph::ptr(new in_mem_undirected_subgraph<edge_data_type>(
					has_data));
	}

	virtual bool is_directed() const {
		return false;
	}

	virtual void add_vertex(const in_mem_vertex &v) {
		const in_mem_undirected_vertex<edge_data_type> &un_v
			= (const in_mem_undirected_vertex<edge_data_type> &) v;
		vertices.insert(typename vmap_t::value_type(v.get_id(), un_v));
		num_edges += un_v.get_num_edges();
		if (un_v.get_num_edges() > 0)
			num_non_empty++;
	}

	virtual const in_mem_vertex &get_vertex(vertex_id_t id) const {
		typename vmap_t::const_iterator it = vertices.find(id);
		assert(it != vertices.end());
		return it->second;
	}

	virtual void get_all_vertices(std::vector<vertex_id_t> &ids) const {
		for (typename vmap_t::const_iterator it = vertices.begin();
				it != vertices.end(); it++)
			ids.push_back(it->first);
	}
	virtual size_t get_num_vertices() const {
		return vertices.size();
	}
};

template<class edge_data_type = empty_data>
class in_mem_directed_subgraph: public in_mem_subgraph
{
	typedef std::map<vertex_id_t, in_mem_directed_vertex<edge_data_type> > vmap_t;
	vmap_t vertices;

	in_mem_directed_subgraph(bool has_data): in_mem_subgraph(has_data) {
	}
public:
	static graph::ptr create(bool has_data) {
		return graph::ptr(new in_mem_directed_subgraph<edge_data_type>(has_data));
	}

	virtual bool is_directed() const {
		return true;
	}

	virtual void add_vertex(const in_mem_vertex &v) {
		const in_mem_directed_vertex<edge_data_type> &dv
			= (const in_mem_directed_vertex<edge_data_type> &) v;
		vertices.insert(typename vmap_t::value_type(v.get_id(), dv));
		num_edges += dv.get_num_edges(edge_type::IN_EDGE);
		if (dv.get_num_edges(edge_type::BOTH_EDGES) > 0)
			num_non_empty++;
	}

	virtual const in_mem_vertex &get_vertex(vertex_id_t id) const {
		typename vmap_t::const_iterator it = vertices.find(id);
		assert(it != vertices.end());
		return it->second;
	}

	virtual void get_all_vertices(std::vector<vertex_id_t> &ids) const {
		for (typename vmap_t::const_iterator it = vertices.begin();
				it != vertices.end(); it++)
			ids.push_back(it->first);
	}
	virtual size_t get_num_vertices() const {
		return vertices.size();
	}
};

inline graph::ptr in_mem_subgraph::create(graph_type type, bool has_data)
{
	if (type == graph_type::DIRECTED)
		return in_mem_directed_subgraph<>::create(has_data);
	else if (type == graph_type::UNDIRECTED)
		return in_mem_undirected_subgraph<>::create(has_data);
	else
		ABORT_MSG("wrong graph type");
}

size_t read_edge_list_text(const std::string &file,
		std::vector<edge<> > &edges);

#if 0
template<class edge_data_type = empty_data>
class undirected_graph: public graph
{
	std::vector<in_mem_undirected_vertex<edge_data_type> > vertices;
	undirected_in_mem_vertex_index  index;

	undirected_graph() {
	}
public:
#if 0
	static undirected_graph *load_edge_list_text(const std::string &file) {
		std::vector<edge<edge_data_type> > edges;
		read_edge_list_text(file, edges);
		return create(edges.data(), edges.size());
	}

	static undirected_graph *load_adjacency_list(const std::string &file) {
		assert(0);
	}
	static ptr create(edge<edge_data_type> edges[], size_t num_edges);
#endif
	static ptr create(bool has_data) {
		return ptr(new undirected_graph<edge_data_type>());
	}

	virtual bool is_directed() const {
		return false;
	}

	bool has_edge_data() const {
		return false;
	}

	void add_vertex(const in_mem_vertex &v1) {
		const in_mem_undirected_vertex<edge_data_type> &v
			= (const in_mem_undirected_vertex<edge_data_type> &) v1;
		vertices.push_back(v);
		index.add_vertex(v);
	}

	virtual const in_mem_vertex &get_vertex(vertex_id_t id) const {
		return vertices[id];
	}

	void get_all_vertices(std::vector<vertex_id_t> &ids) const {
		for (size_t i = 0; i < vertices.size(); i++)
			ids.push_back(vertices[i].get_id());
	}

	void dump(const std::string &index_file,
			const std::string &graph_file) {
		ABORT_MSG("dump isn't implemented");
	}

	size_t get_num_edges() const {
		size_t num_edges = 0;
		for (size_t i = 0; i < vertices.size(); i++)
			num_edges += vertices[i].get_num_edges();
		return num_edges;
	}

	size_t get_num_vertices() const {
		return vertices.size();
	}

	size_t get_num_non_empty_vertices() const {
		size_t num_vertices = 0;
		for (size_t i = 0; i < vertices.size(); i++)
			if (vertices[i].get_num_edges() > 0)
				num_vertices++;
		return num_vertices;
	}

	virtual void print() const {
		ABORT_MSG("print isn't implemented");
	}

	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const {
		ABORT_MSG("check_ext_graph isn't implemented");
	}
};
#endif

template<class edge_data_type>
void check_vertex(const in_mem_directed_vertex<edge_data_type> &in_v,
		ext_mem_undirected_vertex *ext_in_v, ext_mem_undirected_vertex *ext_out_v)
{
	TEST(ext_in_v->get_id() == ext_out_v->get_id());
	TEST(ext_in_v->get_id() == in_v.get_id());
	TEST(ext_in_v->get_num_edges() == in_v.get_num_in_edges());
	TEST(ext_out_v->get_num_edges() == in_v.get_num_out_edges());
	edge_const_iterator<edge_data_type> in_it2
		= in_v.get_in_edge_begin();
	edge_const_iterator<edge_data_type> in_end2
		= in_v.get_in_edge_end();
	for (size_t i = 0; i < ext_in_v->get_num_edges(); i++, ++in_it2) {
		edge<edge_data_type> e2 = *in_it2;
		TEST(ext_in_v->get_neighbor(i) == e2.get_from());
		TEST(ext_in_v->get_id() == e2.get_to());
		if (ext_in_v->has_edge_data())
			TEST(ext_in_v->get_edge_data(i) == e2.get_data());
	}
	TEST(in_it2 == in_end2);

	edge_const_iterator<edge_data_type> out_it2 = in_v.get_out_edge_begin();
	edge_const_iterator<edge_data_type> out_end2 = in_v.get_out_edge_end();
	for (size_t i = 0; i < ext_out_v->get_num_edges(); i++, ++out_it2) {
		edge<edge_data_type> e2 = *out_it2;
		TEST(ext_out_v->get_id() == e2.get_from());
		TEST(ext_out_v->get_neighbor(i) == e2.get_to());
		if (ext_out_v->has_edge_data())
			TEST(ext_out_v->get_edge_data(i) == e2.get_data());
	}
	TEST(out_it2 == out_end2);
}

#if 0
template<class edge_data_type = empty_data>
class directed_graph: public graph
{
	bool has_data;
	size_t num_in_edges;
	size_t num_out_edges;
	size_t num_non_empty_vertices;

	typedef std::pair<vertex_id_t, in_mem_directed_vertex<edge_data_type> > v_pair_t;
	typedef std::map<vertex_id_t, in_mem_directed_vertex<edge_data_type> > v_map_t;
	v_map_t vertices;
	directed_in_mem_vertex_index index;

	directed_graph(bool has_data) {
		this->has_data = has_data;
		num_in_edges = 0;
		num_out_edges = 0;
		num_non_empty_vertices = 0;
	}

	bool exist_vertex(vertex_id_t id) const {
		typename v_map_t::const_iterator it = vertices.find(id);
		return it != vertices.end();
	}
public:
#if 0
	struct delete_graph {
		void operator()(directed_graph<edge_data_type> *g) {
			directed_graph<edge_data_type>::destroy(g);
		}
	};

	typedef std::unique_ptr<directed_graph<edge_data_type>, delete_graph> unique_ptr;

	static unique_ptr load(const std::string &index_file,
			const std::string &graph_file) {
		directed_vertex_index *index = directed_vertex_index::load(index_file);
		
		native_file file(graph_file);
		size_t adj_file_size = file.get_size();
		FILE *adj_f = fopen(graph_file.c_str(), "r");
		assert(adj_f);
		char *adj_buf = new char[adj_file_size];
		size_t ret = fread(adj_buf, adj_file_size, 1, adj_f);
		assert(ret == 1);
		fclose(adj_f);

		graph_header *header = (graph_header *) adj_buf;
		header->verify();
		directed_graph<edge_data_type>::unique_ptr g
			= unique_ptr(new directed_graph<edge_data_type>(
						header->has_edge_data()));
		for (vertex_id_t id = 0; id < index->get_num_vertices(); id++) {
			size_t size = index->get_vertex_size(id);
			off_t off = index->get_vertex_off(id);
			ext_mem_directed_vertex *v = (ext_mem_directed_vertex *) (adj_buf + off);
			assert(v->get_size() == size);
			g->vertices.emplace_back(v);
			assert(g->has_data == v->has_edge_data());
		}
		vertex_index::destroy(index);
		delete [] adj_buf;
		return g;
	}

	static void destroy(directed_graph<edge_data_type> *g) {
		delete g;
	}
#endif
	static ptr create(bool has_data) {
		return ptr(new directed_graph<edge_data_type>(has_data));
	}

	virtual bool is_directed() const {
		return true;
	}

	bool has_edge_data() const {
		return has_data;
	}

	void add_vertex(const in_mem_vertex &v1) {
		const in_mem_directed_vertex<edge_data_type> &v
			= (const in_mem_directed_vertex<edge_data_type> &) v1;
		assert(v.has_edge_data() == has_data);
		std::pair<typename v_map_t::iterator, bool> ret = vertices.insert(v_pair_t(v.get_id(), v));
		BOOST_VERIFY(ret.second);
		num_in_edges += v.get_num_in_edges();
		num_out_edges += v.get_num_out_edges();
		if (v.get_num_edges(edge_type::IN_EDGE)
				+ v.get_num_edges(edge_type::OUT_EDGE) > 0)
			num_non_empty_vertices++;
		index.add_vertex(v);
	}

	virtual const in_mem_vertex &get_vertex(vertex_id_t id) const {
		return vertices[id];
	}

	void get_all_vertices(std::vector<vertex_id_t> &ids) const {
		for (typename v_map_t::const_iterator it = vertices.begin();
				it != vertices.end(); it++)
			ids.push_back(it->second.get_id());
	}

	void dump(const std::string &index_file,
			const std::string &graph_file) {
		assert(!file_exist(index_file));
		assert(!file_exist(graph_file));
		std::string in_graph_file = graph_file + ".in";
		std::string out_graph_file = graph_file + ".out";
		FILE *in_f = fopen(in_graph_file.c_str(), "w");
		FILE *out_f = fopen(out_graph_file.c_str(), "w");
		if (in_f == NULL || out_f == NULL) {
			ABORT_MSG(std::string("fail to open: ") + strerror(errno));
		}

		graph_header header(graph_type::DIRECTED, vertices.size(),
				get_num_edges() / 2, has_data ? sizeof(edge_data_type) : 0);
		BOOST_VERIFY(fwrite(&header, sizeof(header), 1, in_f) == 1);

		std::vector<char> buf;
		size_t in_size = 0;
		size_t out_size = 0;
		for (typename v_map_t::const_iterator it = vertices.begin();
				it != vertices.end(); it++) {
			const in_mem_directed_vertex<edge_data_type> &v = it->second;
			int mem_size = v.get_serialize_size(IN_EDGE);
			buf.resize(mem_size);
			ext_mem_undirected_vertex::serialize(v, buf.data(), mem_size,
					IN_EDGE);
			ssize_t ret = fwrite(buf.data(), mem_size, 1, in_f);
			in_size += mem_size;
			assert(ret == 1);

			mem_size = v.get_serialize_size(OUT_EDGE);
			buf.resize(mem_size);
			ext_mem_undirected_vertex::serialize(v, buf.data(), mem_size,
					OUT_EDGE);
			ret = fwrite(buf.data(), mem_size, 1, out_f);
			out_size += mem_size;
			assert(ret == 1);
		}
		assert(in_size == out_size);
		fclose(out_f);

		// TODO need to append the out-parts to the in-parts
		// and rename the file for in-parts.
		assert(0);
		fclose(in_f);

		index.dump(index_file, header);
	}

	size_t get_num_in_edges() const {
		return num_in_edges;
	}

	size_t get_num_out_edges() const {
		return num_out_edges;
	}

	size_t get_num_edges() const {
		assert(get_num_in_edges() == get_num_out_edges());
		// The total number of edges should be the total number of in-edges
		// or out-edges.
		return get_num_in_edges();
	}

	size_t get_num_vertices() const {
		return vertices.size();
	}

	size_t get_num_non_empty_vertices() const {
		return num_non_empty_vertices;
	}

	virtual void print() const {
		for (typename v_map_t::const_iterator it = vertices.begin();
				it != vertices.end(); it++) {
			const in_mem_directed_vertex<edge_data_type> &v = it->second;
			if (v.get_num_in_edges() + v.get_num_out_edges() > 0)
				v.print();
		}
	}

	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const {
		printf("check the graph in the external memory\n");
		directed_vertex_index::ptr index = directed_vertex_index::load(index_file);
		
		native_file file(adj_file);
		size_t adj_file_size = file.get_size();
		FILE *adj_f = fopen(adj_file.c_str(), "r");
		assert(adj_f);
		char *adj_buf = new char[adj_file_size];
		BOOST_VERIFY(fread(adj_buf, adj_file_size, 1, adj_f) == 1);
		fclose(adj_f);

		graph_header *header = (graph_header *) adj_buf;
		header->verify();
		for (vertex_id_t id = 0; id < index->get_num_vertices(); id++) {
			ext_mem_vertex_info info = index->get_vertex_info_in(id);
			ext_mem_undirected_vertex *in_v
				= (ext_mem_undirected_vertex *) (adj_buf + info.get_off());
			TEST(in_v->get_size() == info.get_size());

			info = index->get_vertex_info_out(id);
			ext_mem_undirected_vertex *out_v
				= (ext_mem_undirected_vertex *) (adj_buf + info.get_off());
			TEST(out_v->get_size() == info.get_size());

			typename v_map_t::const_iterator it = vertices.find(id);
			TEST(it != vertices.end());
			check_vertex(it->second, in_v, out_v);
		}
		delete [] adj_buf;
	}

	virtual void merge(graph::ptr g) {
		directed_graph<edge_data_type> &other
			= (directed_graph<edge_data_type> &) *g;
		assert(this->has_data == other.has_data);
		for (typename v_map_t::const_iterator it = other.vertices.begin();
				it != other.vertices.end(); it++) {
			std::pair<typename v_map_t::iterator, bool> ret = this->vertices.insert(*it);
			BOOST_VERIFY(ret.second);
		}
		this->num_in_edges += other.num_in_edges;
		this->num_out_edges += other.num_out_edges;
		this->num_non_empty_vertices += other.num_non_empty_vertices;
	}

	virtual void dump_as_edge_list(const std::string &file) const {
		assert(!file_exist(file));
		FILE *f = fopen(file.c_str(), "w");
		if (f == NULL)
			ABORT_MSG(boost::format("fail to open %1%: %2%")
					% file % strerror(errno));

		for (typename v_map_t::const_iterator vit = vertices.begin();
				vit != vertices.end(); vit++) {
			edge_const_iterator<edge_data_type> eit
				= vit->second.get_out_edge_begin();
			edge_const_iterator<edge_data_type> end
				= vit->second.get_out_edge_end();
			for (; eit != end; ++eit) {
				edge<edge_data_type> e = *eit;
				assert(e.get_from() == vit->first);
				// We only print the edges inside the graph or subgraph.
				if (exist_vertex(e.get_to()))
					fprintf(f, "%u\t%u\n", e.get_from(), e.get_to());
			}
		}

		fclose(f);
	}
};
#endif

static inline void unique_merge(const std::vector<vertex_id_t> &v1,
		const std::vector<vertex_id_t> &v2, std::vector<vertex_id_t> &v)
{
	std::vector<vertex_id_t>::const_iterator it1 = v1.begin();
	std::vector<vertex_id_t>::const_iterator it2 = v2.begin();
	while (it1 != v1.end() && it2 != v2.end()) {
		if (*it1 > *it2) {
			v.push_back(*it2);
			it2++;
		}
		else if (*it1 < *it2) {
			v.push_back(*it1);
			it1++;
		}
		else {
			v.push_back(*it1);
			it1++;
			it2++;
		}
	}

	while (it1 != v1.end()) {
		v.push_back(*it1);
		it1++;
	}

	while (it2 != v2.end()) {
		v.push_back(*it2);
		it2++;
	}
}

#endif
