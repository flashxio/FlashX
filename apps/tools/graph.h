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
#include <set>

#include "vertex_index.h"
#include "vertex.h"

size_t read_edge_list_text(const std::string &file,
		std::vector<edge<> > &edges);

class in_mem_graph
{
public:
	virtual void get_all_vertices(std::vector<vertex_id_t> &ids) const = 0;
	virtual vertex_index *create_vertex_index() const = 0;
	virtual void dump(const std::string &file) const = 0;
	virtual size_t get_num_edges() const = 0;
	virtual size_t get_num_vertices() const = 0;
	virtual size_t get_num_non_empty_vertices() const = 0;
	virtual void print() const = 0;
};

template<class edge_data_type = empty_data>
class undirected_graph: public in_mem_graph
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

	void get_all_vertices(std::vector<vertex_id_t> &ids) const {
		for (size_t i = 0; i < vertices.size(); i++)
			ids.push_back(vertices[i].get_id());
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

	virtual void print() const {
		assert(0);
	}
};

template<class edge_data_type = empty_data>
class directed_graph: public in_mem_graph
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

	const in_mem_directed_vertex<edge_data_type> *get_vertex(
			vertex_id_t id) const {
		for (size_t i = 0; i < vertices.size(); i++)
			if (vertices[i].get_id() == id)
				return &vertices[i];
		return NULL;
	}

	void get_all_vertices(std::vector<vertex_id_t> &ids) const {
		for (size_t i = 0; i < vertices.size(); i++)
			ids.push_back(vertices[i].get_id());
	}

	vertex_index *create_vertex_index() const {
		return vertex_index::create<in_mem_directed_vertex<edge_data_type> >(
				vertices);
	}

	void dump(const std::string &file) const {
		FILE *f = fopen(file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			assert(0);
		}

		for (size_t i = 0; i < vertices.size(); i++) {
			int mem_size = vertices[i].get_serialize_size();
			char *buf = new char[mem_size];
			ext_mem_directed_vertex::serialize<edge_data_type>(vertices[i],
					buf, mem_size);
			ssize_t ret = fwrite(buf, mem_size, 1, f);
			delete [] buf;
			assert(ret == 1);
		}

		fclose(f);
	}

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

	size_t get_num_edges() const {
		return get_num_in_edges() + get_num_out_edges();
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

	virtual void print() const {
		assert(0);
	}
};

template<class edge_data_type = empty_data>
class ts_directed_graph: public in_mem_graph
{
	std::vector<ts_in_mem_directed_vertex<edge_data_type> > vertices;
public:
	static ts_directed_graph<edge_data_type> *merge_graphs(
			const std::vector<directed_graph<edge_data_type> *> &graphs) {
		// Get all vertex Ids.
		std::set<vertex_id_t> vertex_ids;
		for (unsigned i = 0; i < graphs.size(); i++) {
			std::vector<vertex_id_t> ids;
			graphs[i]->get_all_vertices(ids);
			vertex_ids.insert(ids.begin(), ids.end());
		}

		ts_directed_graph<edge_data_type> *g
			= new ts_directed_graph<edge_data_type>();
		// Construct one time-series vertex at a time.
		for (std::set<vertex_id_t>::const_iterator it = vertex_ids.begin();
				it != vertex_ids.end(); it++) {
			vertex_id_t id = *it;
			ts_in_mem_directed_vertex<edge_data_type> ts_v(id);
			for (unsigned i = 0; i < graphs.size(); i++) {
				const in_mem_directed_vertex<edge_data_type> *v = graphs[i]->get_vertex(id);
				if (v && v->get_num_in_edges() + v->get_num_out_edges() > 0)
					ts_v.add_timestamp(i, *v);
			}
			g->vertices.push_back(ts_v);
		}
		return g;
	}

	void get_all_vertices(std::vector<vertex_id_t> &ids) const {
		for (size_t i = 0; i < vertices.size(); i++)
			ids.push_back(vertices[i].get_id());
	}

	virtual vertex_index *create_vertex_index() const {
		return vertex_index::create<ts_in_mem_directed_vertex<edge_data_type> >(
				vertices);
	}

	virtual void dump(const std::string &file) const {
		FILE *f = fopen(file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			assert(0);
		}

		for (size_t i = 0; i < vertices.size(); i++) {
			int mem_size = vertices[i].get_serialize_size();
			char *buf = new char[mem_size];
			ts_ext_mem_directed_vertex::serialize<edge_data_type>(vertices[i],
					buf, mem_size);
			ssize_t ret = fwrite(buf, mem_size, 1, f);
			delete [] buf;
			assert(ret == 1);
		}

		fclose(f);
	}

	size_t get_num_in_edges() const {
		size_t num = 0;
		for (size_t i = 0; i < vertices.size(); i++)
			num += vertices[i].get_num_in_edges();
		return num;
	}

	size_t get_num_out_edges() const {
		size_t num = 0;
		for (size_t i = 0; i < vertices.size(); i++)
			num += vertices[i].get_num_out_edges();
		return num;
	}

	virtual size_t get_num_edges() const {
		return get_num_in_edges() + get_num_out_edges();
	}

	virtual size_t get_num_vertices() const {
		return vertices.size();
	}

	virtual size_t get_num_non_empty_vertices() const {
		size_t num = 0;
		for (size_t i = 0; i < vertices.size(); i++) {
			if (vertices[i].get_num_edges() > 0)
				num++;
		}
		return num;
	}

	virtual void print() const {
		for (size_t i = 0; i < vertices.size(); i++) {
			if (vertices[i].get_num_edges() > 0)
				vertices[i].print();
		}
	}
};

#endif
