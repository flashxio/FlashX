#ifndef __GRAPH_H__
#define __GRAPH_H__

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

#include <string>
#include <set>

#include "native_file.h"

#include "vertex_index.h"
#include "vertex.h"

namespace fg
{

class graph
{
public:
	typedef std::shared_ptr<graph> ptr;

	virtual ~graph() {
	}

	virtual bool is_directed() const = 0;
	virtual const in_mem_vertex &get_vertex(vertex_id_t id) const = 0;
	virtual in_mem_vertex &get_vertex(vertex_id_t id) = 0;
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

class in_mem_graph;

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
	typedef std::shared_ptr<in_mem_subgraph> ptr;

	static ptr create(graph_type type, bool has_data);

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
	/*
	 * This compresses the subgraph and generates a graph whose vertex IDs
	 * are adjacent to each other.
	 */
	std::pair<std::shared_ptr<in_mem_graph>, std::shared_ptr<vertex_index> > serialize(
			const std::string &name, bool compress) const;
	void compress();

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
	static in_mem_subgraph::ptr create(bool has_data) {
		return in_mem_subgraph::ptr(new in_mem_undirected_subgraph<edge_data_type>(
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

	virtual in_mem_vertex &get_vertex(vertex_id_t id) {
		typename vmap_t::iterator it = vertices.find(id);
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
	static in_mem_subgraph::ptr create(bool has_data) {
		return in_mem_subgraph::ptr(
				new in_mem_directed_subgraph<edge_data_type>(has_data));
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

	virtual in_mem_vertex &get_vertex(vertex_id_t id) {
		typename vmap_t::iterator it = vertices.find(id);
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

}

#endif
