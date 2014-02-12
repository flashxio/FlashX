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
#include "native_file.h"

size_t read_edge_list_text(const std::string &file, std::vector<edge<> > &edges)
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
		edges.push_back(edge<>(from, to));
	}
	fclose(f);
	return edges.size();
}

template<class edge_data_type = empty_data>
struct comp_edge {
	bool operator() (const edge<edge_data_type> &e1,
			const edge<edge_data_type> &e2) {
		if (e1.get_from() == e2.get_from())
			return e1.get_to() < e2.get_to();
		else
			return e1.get_from() < e2.get_from();
	}
};

template<class edge_data_type>
undirected_graph<edge_data_type> *undirected_graph<edge_data_type>::create(
		edge<edge_data_type> edges[], size_t num_edges)
{
	undirected_graph<edge_data_type> *g = new undirected_graph<edge_data_type>();
	// Each edge appears twice and in different directions.
	// When we sort the edges with the first vertex id, we only need
	// a single scan to construct the graph in the form of
	// the adjacency list.
	edge<edge_data_type> *tmp = new edge<edge_data_type>[num_edges * 2];
	for (size_t i = 0; i < num_edges; i++) {
		tmp[2 * i] = edges[i];
		tmp[2 * i + 1] = edge<edge_data_type>(edges[i].get_to(),
				edges[i].get_from());
	}
	edges = tmp;
	num_edges *= 2;
	comp_edge<edge_data_type> edge_comparator;
	std::sort(edges, edges + num_edges, edge_comparator);

	vertex_id_t curr = edges[0].get_from();
	in_mem_undirected_vertex<edge_data_type> v(curr);
	for (size_t i = 0; i < num_edges; i++) {
		vertex_id_t id = edges[i].get_from();
		if (curr == id) {
			// We have to make sure the edge doesn't exist.
			v.add_edge(edges[i]);
		}
		else {
			g->add_vertex(v);
			vertex_id_t prev = curr + 1;
			curr = id;
			// The vertices without edges won't show up in the edge list,
			// but we need to fill the gap in the vertex Id space with empty
			// vertices.
			while (prev < curr) {
				v = in_mem_undirected_vertex<edge_data_type>(prev);
				prev++;
				g->add_vertex(v);
			}
			v = in_mem_undirected_vertex<edge_data_type>(curr);
			v.add_edge(edges[i]);
		}
	}
	g->add_vertex(v);
	delete [] edges;
	return g;
}

template<class edge_data_type>
void undirected_graph<edge_data_type>::dump(const std::string &index_file,
		const std::string &graph_file)
{
	assert(!file_exist(index_file));
	assert(!file_exist(graph_file));
	FILE *f = fopen(graph_file.c_str(), "w");
	if (f == NULL) {
		perror("fopen");
		assert(0);
	}

	graph_header header(graph_type::UNDIRECTED, vertices.size(),
			get_num_edges(), false);
	ssize_t ret = fwrite(&header, sizeof(header), 1, f);
	assert(ret == 1);
	for (size_t i = 0; i < vertices.size(); i++) {
		int mem_size = vertices[i].get_serialize_size();
		char *buf = new char[mem_size];
		ext_mem_undirected_vertex::serialize<edge_data_type>(vertices[i],
				buf, mem_size);
		ssize_t ret = fwrite(buf, mem_size, 1, f);
		delete [] buf;
		assert(ret == 1);
	}

	fclose(f);

	vertex_index *index = create_vertex_index();
	index->dump(index_file);
	vertex_index::destroy(index);
}

template class undirected_graph<>;
