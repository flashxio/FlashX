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

#include "native_file.h"

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
	virtual bool has_edge_data() const = 0;
	virtual size_t get_num_non_empty_vertices() const = 0;
	virtual void print() const = 0;
	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const = 0;
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

	bool has_edge_data() const {
		return false;
	}

	void add_vertex(const in_mem_undirected_vertex<edge_data_type> &v) {
		vertices.push_back(v);
	}

	void get_all_vertices(std::vector<vertex_id_t> &ids) const {
		for (size_t i = 0; i < vertices.size(); i++)
			ids.push_back(vertices[i].get_id());
	}

	vertex_index *create_vertex_index() const {
		graph_header header(graph_type::UNDIRECTED, vertices.size(),
				get_num_edges(), false);
		return vertex_index::create<in_mem_undirected_vertex<edge_data_type> >(
				header, vertices);
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

	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const {
		assert(0);
	}
};

template<class edge_data_type>
void check_vertex(in_mem_directed_vertex<edge_data_type> in_v,
		ext_mem_directed_vertex *ext_v)
{
	assert(ext_v->get_id() == in_v.get_id());
	assert(ext_v->get_num_in_edges() == in_v.get_num_in_edges());
	assert(ext_v->get_num_out_edges() == in_v.get_num_out_edges());
	edge_const_iterator<edge_data_type> in_it1
		= ext_v->get_in_edge_begin<edge_data_type>();
	edge_const_iterator<edge_data_type> in_end1
		= ext_v->get_in_edge_end<edge_data_type>();
	edge_const_iterator<edge_data_type> in_it2
		= in_v.get_in_edge_begin();
	edge_const_iterator<edge_data_type> in_end2
		= in_v.get_in_edge_end();
	while (in_it1 != in_end1 && in_it2 != in_end2) {
		edge<edge_data_type> e1 = *in_it1;
		edge<edge_data_type> e2 = *in_it2;
		assert(e1.get_from() == e2.get_from());
		assert(e1.get_to() == e2.get_to());
		++in_it1;
		++in_it2;
	}
	assert(in_it1 == in_end1 && in_it2 == in_end2);

	edge_const_iterator<edge_data_type> out_it1
		= ext_v->get_out_edge_begin<edge_data_type>();
	edge_const_iterator<edge_data_type> out_end1
		= ext_v->get_out_edge_end<edge_data_type>();
	edge_const_iterator<edge_data_type> out_it2 = in_v.get_out_edge_begin();
	edge_const_iterator<edge_data_type> out_end2 = in_v.get_out_edge_end();
	while (out_it1 != out_end1 && out_it2 != out_end2) {
		edge<edge_data_type> e1 = *out_it1;
		edge<edge_data_type> e2 = *out_it2;
		assert(e1.get_from() == e2.get_from());
		assert(e1.get_to() == e2.get_to());
		++out_it1;
		++out_it2;
	}
	assert(out_it1 == out_end1 && out_it2 == out_end2);
}

template<class edge_data_type = empty_data>
class directed_graph: public in_mem_graph
{
	bool has_data;
	std::vector<in_mem_directed_vertex<edge_data_type> > vertices;
public:
	static void destroy(directed_graph *g) {
		delete g;
	}

	directed_graph(bool has_data) {
		this->has_data = has_data;
	}

	bool has_edge_data() const {
		return has_data;
	}

	void add_vertex(const in_mem_directed_vertex<edge_data_type> &v) {
		assert(vertices.size() == v.get_id());
		assert(v.has_edge_data() == has_data);
		vertices.push_back(v);
	}

	typename std::vector<in_mem_directed_vertex<edge_data_type> >::const_iterator begin() const {
		return vertices.begin();
	}

	typename std::vector<in_mem_directed_vertex<edge_data_type> >::const_iterator end() const {
		return vertices.end();
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
		graph_header header(graph_type::DIRECTED, vertices.size(),
				get_num_edges(), has_data);
		return vertex_index::create<in_mem_directed_vertex<edge_data_type> >(
				header, vertices);
	}

	void dump(const std::string &file) const {
		FILE *f = fopen(file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			assert(0);
		}

		graph_header header(graph_type::DIRECTED, vertices.size(),
				get_num_edges(), has_data);
		ssize_t ret = fwrite(&header, sizeof(header), 1, f);
		assert(ret == 1);

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
		for (size_t i = 0; i < vertices.size(); i++) {
			if (vertices[i].get_num_in_edges()
					+ vertices[i].get_num_out_edges() > 0)
				vertices[i].print();
		}
	}

	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const {
		printf("check the graph in the external memory\n");
		vertex_index *index = vertex_index::load(index_file);
		
		native_file file(adj_file);
		size_t adj_file_size = file.get_size();
		FILE *adj_f = fopen(adj_file.c_str(), "r");
		assert(adj_f);
		char *adj_buf = new char[adj_file_size];
		size_t ret = fread(adj_buf, adj_file_size, 1, adj_f);
		assert(ret == 1);
		fclose(adj_f);

		graph_header *header = (graph_header *) adj_buf;
		header->verify();
		for (vertex_id_t id = 0; id < index->get_num_vertices(); id++) {
			int size = index->get_vertex_size(id);
			off_t off = index->get_vertex_off(id);
			ext_mem_directed_vertex *v = (ext_mem_directed_vertex *) (adj_buf + off);
			assert(v->get_size() == (size_t) size);
			check_vertex(vertices[id], v);
		}
		vertex_index::destroy(index);
		delete [] adj_buf;
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

template<class edge_data_type>
void check_vertex(ts_in_mem_directed_vertex<edge_data_type> in_v,
		ts_ext_mem_directed_vertex *ext_v)
{
	assert(ext_v->get_id() == in_v.get_id());
	assert(ext_v->get_num_edges() == in_v.get_num_edges());
	assert(ext_v->get_num_timestamps() == in_v.get_num_timestamps());
	std::vector<int> all_timestamps;
	in_v.get_all_timestamps(all_timestamps);
	assert(all_timestamps.size() == (size_t) ext_v->get_num_timestamps());
	for (std::vector<int>::const_iterator it = all_timestamps.begin();
			it != all_timestamps.end(); it++) {
		assert(ext_v->get_num_in_edges(*it)
				== in_v.get_num_in_edges(*it));
		assert(ext_v->get_num_out_edges(*it)
				== in_v.get_num_out_edges(*it));
		edge_const_iterator<edge_data_type> in_it1
			= ext_v->get_in_edge_begin<edge_data_type>(*it);
		edge_const_iterator<edge_data_type> in_end1
			= ext_v->get_in_edge_end<edge_data_type>(*it);
		edge_const_iterator<edge_data_type> in_it2 = in_v.get_in_edge_begin(*it);
		edge_const_iterator<edge_data_type> in_end2 = in_v.get_in_edge_end(*it);
		while (in_it1 != in_end1 && in_it2 != in_end2) {
			edge<edge_data_type> e1 = *in_it1;
			edge<edge_data_type> e2 = *in_it2;
			assert(e1.get_from() == e2.get_from());
			assert(e1.get_to() == e2.get_to());
			++in_it1;
			++in_it2;
		}
		assert(in_it1 == in_end1 && in_it2 == in_end2);

		edge_const_iterator<edge_data_type> out_it1
			= ext_v->get_out_edge_begin<edge_data_type>(*it);
		edge_const_iterator<edge_data_type> out_end1
			= ext_v->get_out_edge_end<edge_data_type>(*it);
		edge_const_iterator<edge_data_type> out_it2
			= in_v.get_out_edge_begin(*it);
		edge_const_iterator<edge_data_type> out_end2
			= in_v.get_out_edge_end(*it);
		while (out_it1 != out_end1 && out_it2 != out_end2) {
			edge<edge_data_type> e1 = *out_it1;
			edge<edge_data_type> e2 = *out_it2;
			assert(e1.get_from() == e2.get_from());
			assert(e1.get_to() == e2.get_to());
			++out_it1;
			++out_it2;
		}
		assert(out_it1 == out_end1 && out_it2 == out_end2);
	}
}

template<class edge_data_type = empty_data>
class ts_directed_graph: public in_mem_graph
{
	int max_num_timestamps;
	bool has_data;
	std::vector<ts_in_mem_directed_vertex<edge_data_type> > vertices;

	ts_directed_graph() {
		has_data = false;
		max_num_timestamps = 0;
	}
public:
	static ts_directed_graph<edge_data_type> *merge_graphs(
			const std::vector<directed_graph<edge_data_type> *> &graphs) {
		if (graphs.empty())
			return NULL;

		// Get all vertex Ids.
		bool has_edge_data = graphs[0]->has_edge_data();
		std::vector<vertex_id_t> vertex_ids;
		for (unsigned i = 0; i < graphs.size(); i++) {
			std::vector<vertex_id_t> ids;
			// The vertices in the graph should have been sorted.
			graphs[i]->get_all_vertices(ids);
			std::vector<vertex_id_t> tmp;
			unique_merge(vertex_ids, ids, tmp);
			vertex_ids = tmp;
			assert(has_edge_data == graphs[i]->has_edge_data());
		}

		ts_directed_graph<edge_data_type> *g
			= new ts_directed_graph<edge_data_type>();
		g->has_data = has_edge_data;
		g->max_num_timestamps = graphs.size();
		std::vector<typename std::vector<in_mem_directed_vertex<edge_data_type> >::const_iterator> its;
		for (unsigned i = 0; i < graphs.size(); i++) {
			its.push_back(graphs[i]->begin());
		}
		// Construct one time-series vertex at a time.
		for (std::vector<vertex_id_t>::const_iterator it = vertex_ids.begin();
				it != vertex_ids.end(); it++) {
			vertex_id_t id = *it;
			ts_in_mem_directed_vertex<edge_data_type> ts_v(id, has_edge_data);
			for (unsigned i = 0; i < its.size(); i++) {
				if (its[i] == graphs[i]->end())
					continue;
				assert(its[i]->get_id() >= id);
				if (its[i]->get_id() == id
						&& its[i]->get_num_in_edges() + its[i]->get_num_out_edges() > 0)
					ts_v.add_timestamp(i, *its[i]);
				if (its[i]->get_id() == id)
					its[i]++;
			}
			g->add_vertex(ts_v);
		}
		return g;
	}

	void add_vertex(const ts_in_mem_directed_vertex<edge_data_type> &v) {
		vertices.push_back(v);
	}

	bool has_edge_data() const {
		return has_data;
	}

	void get_all_vertices(std::vector<vertex_id_t> &ids) const {
		for (size_t i = 0; i < vertices.size(); i++)
			ids.push_back(vertices[i].get_id());
	}

	virtual vertex_index *create_vertex_index() const {
		graph_header header(TS_DIRECTED, vertices.size(), get_num_edges(),
				has_data, max_num_timestamps);
		return vertex_index::create<ts_in_mem_directed_vertex<edge_data_type> >(
				header, vertices);
	}

	virtual void dump(const std::string &file) const {
		FILE *f = fopen(file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			assert(0);
		}

		graph_header header(graph_type::TS_DIRECTED, vertices.size(),
				get_num_edges(), has_data, max_num_timestamps);
		ssize_t ret = fwrite(&header, sizeof(header), 1, f);
		assert(ret == 1);

		for (size_t i = 0; i < vertices.size(); i++) {
			int mem_size = vertices[i].get_serialize_size();
			char *buf = new char[mem_size];
			ts_ext_mem_directed_vertex::serialize<edge_data_type>(vertices[i],
					buf, mem_size);

			// Test the correctness of ts_ext_mem_directed_vertex.
			ts_ext_mem_directed_vertex *ext_v = (ts_ext_mem_directed_vertex *) buf;
			check_vertex(vertices[i], ext_v);

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

	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const {
		printf("check the graph in the external memory\n");
		vertex_index *index = vertex_index::load(index_file);
		
		native_file file(adj_file);
		size_t adj_file_size = file.get_size();
		FILE *adj_f = fopen(adj_file.c_str(), "r");
		assert(adj_f);
		char *adj_buf = new char[adj_file_size];
		size_t ret = fread(adj_buf, adj_file_size, 1, adj_f);
		assert(ret == 1);
		fclose(adj_f);

		graph_header *header = (graph_header *) adj_buf;
		header->verify();
		for (vertex_id_t id = 0; id < index->get_num_vertices(); id++) {
			int size = index->get_vertex_size(id);
			off_t off = index->get_vertex_off(id);
			ts_ext_mem_directed_vertex *v = (ts_ext_mem_directed_vertex *) (adj_buf + off);
			assert(v->get_size() == (size_t) size);
			check_vertex(vertices[id], v);
		}
		vertex_index::destroy(index);
		delete [] adj_buf;
	}
};

#endif
