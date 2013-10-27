#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <string>

#include "vertex_index.h"
#include "vertex.h"

class undirected_graph
{
	std::vector<in_mem_undirected_vertex> vertices;

	static undirected_graph *create(edge edges[], size_t num_edges);
public:
	static undirected_graph *load_edge_list(const std::string &file);
	static undirected_graph *load_edge_list_text(const std::string &file);
	static undirected_graph *load_adjacency_list(const std::string &file);
	static void destroy(undirected_graph *g) {
		delete g;
	}

	void add_vertex(const in_mem_undirected_vertex &v) {
		vertices.push_back(v);
	}

	vertex_index *create_vertex_index() const;
	void dump(const std::string &file) const;

	size_t get_num_edges() const;
	size_t get_num_vertices() const {
		return vertices.size();
	}
	size_t get_num_non_empty_vertices() const;
};

class directed_graph
{
	std::vector<in_mem_directed_vertex> vertices;

	static directed_graph *create(edge edges[], size_t num_edges);
public:
	static directed_graph *load_edge_list(const std::string &file);
	static directed_graph *load_edge_list_text(const std::string &file);
	static directed_graph *load_adjacency_list(const std::string &file);
	static void destroy(directed_graph *g) {
		delete g;
	}

	void add_vertex(const in_mem_directed_vertex &v) {
		vertices.push_back(v);
	}

	vertex_index *create_vertex_index() const;
	void dump(const std::string &file) const;

	size_t get_num_in_edges() const;
	size_t get_num_out_edges() const;
	size_t get_num_vertices() const {
		return vertices.size();
	}
	size_t get_num_non_empty_vertices() const;
};

#endif
