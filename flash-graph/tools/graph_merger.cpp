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
#include "edge_type.h"
#include "ext_mem_vertex_iterator.h"

#if 0

bool print_graph = false;
bool check_graph = false;

template<class in_vertex_type, class out_vertex_type>
void merge_dump(std::vector<ext_mem_vertex_iterator> &its, graph_type type,
		const std::string &adjacency_list_file, const std::string &index_file)
{
	FILE *adjacency_fd = fopen(adjacency_list_file.c_str(), "w");
	assert(adjacency_fd);
	// We don't know enough information about the graph, let's just
	// write some dump data to occupy the space for the graph header.
	graph_header header;
	size_t ret = fwrite(&header, sizeof(header), 1, adjacency_fd);
	assert(ret == 1);

	in_mem_vertex_index *in_mem_index;
	if (type == graph_type::DIRECTED)
		in_mem_index = new directed_in_mem_vertex_index();
	else if (type == graph_type::TS_DIRECTED)
		in_mem_index = new ts_directed_in_mem_vertex_index();
	else
		assert(0);

	bool has_next;
	size_t tot_vertices = 0;
	size_t tot_non_empty = 0;
	size_t tot_edges = 0;
	do {
		has_next = false;

		// Get all vertices of the same vertex id.
		std::vector<const in_vertex_type *> vertices;
		for (size_t i = 0; i < its.size(); i++) {
			if (its[i].has_next()) {
				in_vertex_type *v = its[i].next<in_vertex_type>();
				vertices.push_back(v);
				has_next |= its[i].has_next();
			}
			else
				vertices.push_back(NULL);
		}
		assert(vertices.size() == its.size());
		tot_vertices++;

		// Merge the vertices
		typename out_vertex_type::unique_ptr v = out_vertex_type::merge(vertices);
		int num_edges = v->get_num_in_edges() + v->get_num_out_edges();
		if (num_edges > 0) {
			tot_non_empty++;
			tot_edges += num_edges;
		}
		
		// Output the merged vertex.
		in_mem_index->add_vertex((char *) v.get());
		size_t ret = fwrite((char *) v.get(), v->get_size(), 1, adjacency_fd);
		assert(ret);
	} while (has_next);
	tot_edges /= 2;

	assert(its.size() > 0);
	bool has_edge_data = its[0].get_graph_header().has_edge_data();
	for (size_t i = 1; i < its.size(); i++)
		assert(its[i].get_graph_header().has_edge_data() == has_edge_data);

	int num_timestamps = 0;
	if (type == graph_type::TS_DIRECTED || type == graph_type::TS_UNDIRECTED)
		num_timestamps = its.size();
	header = graph_header(type, tot_vertices, tot_edges, has_edge_data,
			num_timestamps);
	int seek_ret = fseek(adjacency_fd, 0, SEEK_SET);
	assert(seek_ret == 0);
	ret = fwrite(&header, sizeof(header), 1, adjacency_fd);
	assert(ret == 1);
	fclose(adjacency_fd);

	in_mem_index->dump(index_file, header);
	delete in_mem_index;

	printf("There are %ld vertices, %ld non-empty vertices and %ld edges\n",
			tot_vertices, tot_non_empty, tot_edges);
}

#if 0
template<class edge_data_type>
void verify_ts_directed_graph(const std::vector<std::string> &in_graph_files,
		const std::vector<std::string> &in_index_files,
		const std::string &out_graph_file, const std::string &out_index_file)
{
	assert(in_graph_files.size() == in_index_files.size());
	std::vector<typename directed_graph<edge_data_type>::unique_ptr> graphs;
	for (size_t i = 0; i < in_graph_files.size(); i++) {
		graphs.emplace_back(directed_graph<edge_data_type>::load(
					in_index_files[i], in_graph_files[i]));
	}
	typename ts_directed_graph<edge_data_type>::unique_ptr ts_g
		= ts_directed_graph<edge_data_type>::merge_graphs(graphs);
	ts_g->check_ext_graph(out_index_file, out_graph_file);
}

template<class edge_data_type = empty_data>
void verify_graph(const std::vector<std::string> &in_graph_files,
		const std::vector<std::string> &in_index_files,
		const std::string &out_graph_file, const std::string &out_index_file,
		graph_type out_type)
{
	switch(out_type) {
		case graph_type::TS_DIRECTED:
			verify_ts_directed_graph<edge_count>(in_graph_files,
					in_index_files, out_graph_file, out_index_file);
			break;
		default:
			assert(0);
	}
}
#endif

void print_usage()
{
	fprintf(stderr, "merge multiple graphs to a single graph\n");
	fprintf(stderr,
			"graph_merger [options] adj_list_file index_file graph_file_list index_file_list\n");
	fprintf(stderr, "-p: print adjacency list\n");
	fprintf(stderr, "-v: verify the created adjacency list\n");
	fprintf(stderr, "-t: merge graphs into a time-series graph. \n");
	fprintf(stderr, "-T type: the type of edge data. Supported type: ");
	for (int i = 0; i < type_map_size; i++) {
		fprintf(stderr, "%s, ", edge_type_map[i].str.c_str());
	}
	fprintf(stderr, "\n");
}

size_t read_file_list(const std::string &file, std::vector<std::string> &file_list)
{
	FILE *f = fopen(file.c_str(), "r");
	assert(f);

	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	while ((read = getline(&line, &len, f)) > 0) {
		if (line[read] == '\n')
			line[read] = 0;
		if (read > 0 && line[read - 1] == '\n')
			line[read - 1] = 0;
		file_list.push_back(line);
	}

	fclose(f);
	return file_list.size();
}

int main(int argc, char *argv[])
{
	int opt;
	int num_opts = 0;
	bool ts_merge = false;
	char *type_str = NULL;
	while ((opt = getopt(argc, argv, "pvtT:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'p':
				print_graph = true;
				break;
			case 'v':
				check_graph = true;
				break;
			case 't':
				ts_merge = true;
				break;
			case 'T':
				type_str = optarg;
				num_opts++;
				break;
			default:
				print_usage();
		}
	}
	argv += 1 + num_opts;
	argc -= 1 + num_opts;
	if (argc < 4) {
		print_usage();
		exit(-1);
	}

	std::string out_graph_file = argv[0];
	out_graph_file += std::string("-v") + itoa(CURR_VERSION);
	std::string out_index_file = argv[1];
	out_index_file += std::string("-v") + itoa(CURR_VERSION);
	std::string graph_list_file = argv[2];
	std::string index_list_file = argv[3];

	std::vector<std::string> in_graph_files;
	std::vector<std::string> in_index_files;
	read_file_list(graph_list_file, in_graph_files);
	read_file_list(index_list_file, in_index_files);
	assert(in_graph_files.size() == in_index_files.size());

	// Initialize vertex iterators.
	std::vector<ext_mem_vertex_iterator> its(in_index_files.size());
	for (size_t i = 0; i < in_index_files.size(); i++)
		its[i].init(in_index_files[i], in_graph_files[i]);

	// Check the input graphs.
	bool directed = its[0].get_graph_header().is_directed_graph();
	graph_type in_type = its[0].get_graph_header().get_graph_type();
	// We don't want to merge two time-series graphs.
	assert(in_type == graph_type::DIRECTED
			|| in_type == graph_type::UNDIRECTED);
	for (size_t i = 0; i < its.size(); i++) {
		assert(its[i].get_graph_header().is_directed_graph() == directed);
		in_type = its[i].get_graph_header().get_graph_type();
		assert(in_type == graph_type::DIRECTED
				|| in_type == graph_type::UNDIRECTED);
	}

	// If we get more than one graph, we need to merge them and will have
	// a time-series graph.
	if (its.size() > 1) {
		struct timeval start, end;
		gettimeofday(&start, NULL);

		graph_type out_type;
		if (ts_merge && directed) {
			out_type = graph_type::TS_DIRECTED;
			merge_dump<ext_mem_directed_vertex, ts_ext_mem_directed_vertex>(its,
					out_type, out_graph_file, out_index_file);
		}
		else if (ts_merge && !directed) {
			assert(0);
		}
		else if (directed) {
			out_type = graph_type::DIRECTED;
			merge_dump<ext_mem_directed_vertex, ext_mem_directed_vertex>(its,
					out_type, out_graph_file, out_index_file);
		}
		else {
			assert(0);
		}

		gettimeofday(&end, NULL);
		printf("It takes %f seconds to merge graphs\n",
				time_diff(start, end));

#if 0
		if (check_graph) {
			int edge_type = DEFAULT_TYPE;
			if (type_str)
				edge_type = conv_edge_type_str2int(type_str);
			switch(edge_type) {
				case DEFAULT_TYPE:
					verify_graph<>(in_graph_files, in_index_files,
							out_graph_file, out_index_file, out_type);
					break;
				case EDGE_COUNT:
					verify_graph<edge_count>(in_graph_files, in_index_files,
							out_graph_file, out_index_file, out_type);
					break;
				default:
					assert(0);
			}
		}
#endif
	}
}
#endif
