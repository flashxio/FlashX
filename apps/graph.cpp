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

#include <algorithm>

#include "graph.h"
#include "common.h"

size_t read_edge_list(const std::string &file, std::vector<edge> &edge_vector)
{
	ssize_t file_size = get_file_size(file.c_str());
	assert(file_size > (ssize_t) sizeof(edge));
	assert(file_size % sizeof(edge) == 0);
	int num_edges = file_size / sizeof(edge);
	edge *edges = new edge[num_edges];
	FILE *f = fopen(file.c_str(), "r");
	if (f == NULL) {
		perror("fopen");
		assert(0);
	}
	size_t ret = fread(edges, file_size, 1, f);
	assert(ret == 1);
	fclose(f);
	edge_vector.assign(edges, edges + num_edges);
	delete [] edges;
	return edge_vector.size();
}

size_t read_edge_list_text(const std::string &file, std::vector<edge> &edges)
{
	FILE *f = fopen(file.c_str(), "r");
	if (f == NULL) {
		perror("fopen");
		assert(0);
	}

	ssize_t read;
	size_t len = 0;
	char *line = NULL;
	while ((read = getline(&line, &len, f)) != -1) {
		if (line[read - 1] == '\n')
			line[read - 1] = 0;
		if (line[read - 2] == '\r')
			line[read - 2] = 0;
		if (line[0] == '#')
			continue;
		char *second = strchr(line, '\t');
		assert(second);
		*second = 0;
		second++;
		if (!isnumeric(line) || !isnumeric(second)) {
			printf("%s\t%s\n", line, second);
			continue;
		}
		vertex_id_t from = atol(line);
		vertex_id_t to = atol(second);
		edges.push_back(edge(from, to));
	}
	fclose(f);
	return edges.size();
}

static struct comp_edge {
	bool operator() (const edge &e1, const edge &e2) {
		if (e1.get_from() == e2.get_from())
			return e1.get_to() < e2.get_to();
		else
			return e1.get_from() < e2.get_from();
	}
} edge_comparator;

undirected_graph *undirected_graph::create(edge edges[], size_t num_edges)
{
	undirected_graph *g = new undirected_graph();
	// Each edge appears twice and in different directions.
	// When we sort the edges with the first vertex id, we only need
	// a single scan to construct the graph in the form of
	// the adjacency list.
	edge *tmp = new edge[num_edges * 2];
	for (size_t i = 0; i < num_edges; i++) {
		tmp[2 * i] = edges[i];
		tmp[2 * i + 1] = edge(edges[i].get_to(), edges[i].get_from());
	}
	edges = tmp;
	num_edges *= 2;
	std::sort(edges, edges + num_edges, edge_comparator);

	vertex_id_t curr = edges[0].get_from();
	in_mem_undirected_vertex v(curr);
	for (size_t i = 0; i < num_edges; i++) {
		vertex_id_t id = edges[i].get_from();
		if (curr == id) {
			// We have to make sure the edge doesn't exist.
			assert(!v.has_edge(edges[i].get_to()));
			v.add_edge(edges[i].get_to());
		}
		else {
			g->add_vertex(v);
			vertex_id_t prev = curr + 1;
			curr = id;
			// The vertices without edges won't show up in the edge list,
			// but we need to fill the gap in the vertex Id space with empty
			// vertices.
			while (prev < curr) {
				v = in_mem_undirected_vertex(prev);
				prev++;
				g->add_vertex(v);
			}
			v = in_mem_undirected_vertex(curr);
			v.add_edge(edges[i].get_to());
		}
	}
	g->add_vertex(v);
	delete [] edges;
	return g;
}

undirected_graph *undirected_graph::load_edge_list(const std::string &file)
{
	std::vector<edge> edges;
	read_edge_list(file, edges);
	return create(edges.data(), edges.size());
}

undirected_graph *undirected_graph::load_edge_list_text(const std::string &file)
{
	std::vector<edge> edges;
	read_edge_list_text(file, edges);
	return create(edges.data(), edges.size());
}

undirected_graph *undirected_graph::load_adjacency_list(const std::string &file)
{
	assert(0);
}

vertex_index *undirected_graph::create_vertex_index() const
{
	return vertex_index::create<in_mem_undirected_vertex>(vertices);
}

void undirected_graph::dump(const std::string &file) const
{
	FILE *f = fopen(file.c_str(), "w");
	if (f == NULL) {
		perror("fopen");
		assert(0);
	}

	for (size_t i = 0; i < vertices.size(); i++) {
		int mem_size = vertices[i].get_serialize_size();
		char *buf = new char[mem_size];
		ext_mem_undirected_vertex::serialize(vertices[i], buf, mem_size);
		ssize_t ret = fwrite(buf, mem_size, 1, f);
		delete [] buf;
		assert(ret == 1);
	}

	fclose(f);
}

size_t undirected_graph::get_num_edges() const
{
	size_t num_edges = 0;
	for (size_t i = 0; i < vertices.size(); i++)
		num_edges += vertices[i].get_num_edges();
	return num_edges;
}

size_t undirected_graph::get_num_non_empty_vertices() const
{
	size_t num_vertices = 0;
	for (size_t i = 0; i < vertices.size(); i++)
		if (vertices[i].get_num_edges() > 0)
			num_vertices++;
	return num_vertices;
}

static struct comp_in_edge {
	bool operator() (const edge &e1, const edge &e2) {
		if (e1.get_to() == e2.get_to())
			return e1.get_from() < e2.get_from();
		else
			return e1.get_to() < e2.get_to();
	}
} in_edge_comparator;

directed_graph *directed_graph::create(edge edges[], size_t num_edges)
{
	directed_graph *g = new directed_graph();

	std::sort(edges, edges + num_edges, edge_comparator);
	edge *sorted_out_edges = edges;

	edge *copied_edges = new edge[num_edges];
	memcpy(copied_edges, edges, num_edges * sizeof(edges[0]));
	std::sort(copied_edges, copied_edges + num_edges, in_edge_comparator);
	edge *sorted_in_edges = copied_edges;

	vertex_id_t curr = min(sorted_out_edges[0].get_from(),
			sorted_in_edges[0].get_to());
	in_mem_directed_vertex v(curr);
	size_t out_idx = 0;
	size_t in_idx = 0;
	while (out_idx < num_edges && in_idx < num_edges) {
		while (sorted_out_edges[out_idx].get_from() == curr
				&& out_idx < num_edges) {
			v.add_out_edge(sorted_out_edges[out_idx++].get_to());
		}
		while (sorted_in_edges[in_idx].get_to() == curr
				&& in_idx < num_edges) {
			v.add_in_edge(sorted_in_edges[in_idx++].get_from());
		}
		g->add_vertex(v);
		vertex_id_t prev = curr + 1;
		if (out_idx < num_edges && in_idx < num_edges)
			curr = min(sorted_out_edges[out_idx].get_from(),
					sorted_in_edges[in_idx].get_to());
		else if (out_idx < num_edges)
			curr = sorted_out_edges[out_idx].get_from();
		else if (in_idx < num_edges)
			curr = sorted_in_edges[in_idx].get_to();
		else
			break;
		// The vertices without edges won't show up in the edge list,
		// but we need to fill the gap in the vertex Id space with empty
		// vertices.
		while (prev < curr) {
			v = in_mem_directed_vertex(prev);
			prev++;
			g->add_vertex(v);
		}
		v = in_mem_directed_vertex(curr);
	}

	assert(g->get_num_in_edges() == num_edges);
	assert(g->get_num_out_edges() == num_edges);
	delete [] copied_edges;
	return g;
}

directed_graph *directed_graph::load_edge_list(const std::string &file)
{
	std::vector<edge> edges;
	read_edge_list(file, edges);
	return create(edges.data(), edges.size());
}

directed_graph *directed_graph::load_edge_list_text(const std::string &file)
{
	std::vector<edge> edges;
	read_edge_list_text(file, edges);
	return create(edges.data(), edges.size());
}

directed_graph *directed_graph::load_adjacency_list(const std::string &file)
{
	assert(0);
	return NULL;
}

vertex_index *directed_graph::create_vertex_index() const
{
	return vertex_index::create<in_mem_directed_vertex>(vertices);
}

void directed_graph::dump(const std::string &file) const
{
	FILE *f = fopen(file.c_str(), "w");
	if (f == NULL) {
		perror("fopen");
		assert(0);
	}

	for (size_t i = 0; i < vertices.size(); i++) {
		int mem_size = vertices[i].get_serialize_size();
		char *buf = new char[mem_size];
		ext_mem_directed_vertex::serialize(vertices[i], buf, mem_size);
		ssize_t ret = fwrite(buf, mem_size, 1, f);
		delete [] buf;
		assert(ret == 1);
	}

	fclose(f);
}

size_t directed_graph::get_num_in_edges() const
{
	size_t num_in_edges = 0;
	for (size_t i = 0; i < vertices.size(); i++)
		num_in_edges += vertices[i].get_num_in_edges();
	return num_in_edges;
}

size_t directed_graph::get_num_out_edges() const
{
	size_t num_out_edges = 0;
	for (size_t i = 0; i < vertices.size(); i++)
		num_out_edges += vertices[i].get_num_out_edges();
	return num_out_edges;
}

size_t directed_graph::get_num_non_empty_vertices() const
{
	size_t num_vertices = 0;
	for (size_t i = 0; i < vertices.size(); i++)
		if (vertices[i].get_num_in_edges() > 0
				|| vertices[i].get_num_out_edges() > 0)
			num_vertices++;
	return num_vertices;
}
