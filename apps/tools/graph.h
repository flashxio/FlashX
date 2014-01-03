#ifndef __GRAPH_H__
#define __GRAPH_H__

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

#include <string>

#include "vertex_index.h"
#include "vertex.h"

size_t read_edge_list_text(const std::string &file,
		std::vector<edge<> > &edges);

template<class edge_data_type = empty_data>
class undirected_graph
{
	std::vector<in_mem_undirected_vertex<edge_data_type> > vertices;

	static undirected_graph *create(edge<edge_data_type> edges[], size_t num_edges);
public:
	static undirected_graph *load_edge_list_text(const std::string &file) {
		std::vector<edge<edge_data_type> > edges;
		read_edge_list_text(file, edges);
		return create(edges.data(), edges.size());
	}

	static undirected_graph *load_adjacency_list(const std::string &file) {
		assert(0);
	}

	static void destroy(undirected_graph *g) {
		delete g;
	}

	void add_vertex(const in_mem_undirected_vertex<edge_data_type> &v) {
		vertices.push_back(v);
	}

	vertex_index *create_vertex_index() const {
		return vertex_index::create<in_mem_undirected_vertex<edge_data_type> >(
				vertices);
	}

	void dump(const std::string &file) const;

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
};

template<class edge_data_type = empty_data>
class directed_graph
{
	std::vector<in_mem_directed_vertex<edge_data_type> > vertices;
public:
	static void destroy(directed_graph *g) {
		delete g;
	}

	void add_vertex(const in_mem_directed_vertex<edge_data_type> &v) {
		assert(vertices.size() == v.get_id());
		vertices.push_back(v);
	}

	vertex_index *create_vertex_index() const {
		return vertex_index::create<in_mem_directed_vertex<edge_data_type> >(
				vertices);
	}

	void dump(const std::string &file) const;

	size_t get_num_in_edges() const {
		size_t num_in_edges = 0;
		for (size_t i = 0; i < vertices.size(); i++)
			num_in_edges += vertices[i].get_num_in_edges();
		return num_in_edges;
	}

	size_t get_num_out_edges() const {
		size_t num_out_edges = 0;
		for (size_t i = 0; i < vertices.size(); i++)
			num_out_edges += vertices[i].get_num_out_edges();
		return num_out_edges;
	}

	size_t get_num_vertices() const {
		return vertices.size();
	}

	size_t get_num_non_empty_vertices() const {
		size_t num_vertices = 0;
		for (size_t i = 0; i < vertices.size(); i++)
			if (vertices[i].get_num_in_edges() > 0
					|| vertices[i].get_num_out_edges() > 0)
				num_vertices++;
		return num_vertices;
	}
};

#endif
