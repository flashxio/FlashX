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

#include <libgen.h>

#include "common.h"
#include "native_file.h"

#include "utils.h"

struct str2int_pair {
	std::string str;
	int number;
};
static struct str2int_pair edge_type_map[] = {
	{"count", EDGE_COUNT},
	{"timestamp", EDGE_TIMESTAMP},
};
static int type_map_size = sizeof(edge_type_map) / sizeof(edge_type_map[0]);

static inline int conv_edge_type_str2int(const std::string &type_str)
{
	for (int i = 0; i < type_map_size; i++) {
		if (edge_type_map[i].str == type_str) {
			return edge_type_map[i].number;
		}
	}
	return DEFAULT_TYPE;
}

static bool check_graph = false;
static std::string work_dir = ".";

void print_usage()
{
	fprintf(stderr, "convert an edge list to adjacency lists\n");
	fprintf(stderr,
			"el2al [options] adj_list_file index_file edge_list_files (or directories)\n");
	fprintf(stderr, "-u: undirected graph\n");
	fprintf(stderr, "-v: verify the created adjacency list\n");
	fprintf(stderr, "-t type: the type of edge data. Supported type: ");
	for (int i = 0; i < type_map_size; i++) {
		fprintf(stderr, "%s, ", edge_type_map[i].str.c_str());
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "-m: merge multiple edge lists into a single graph. \n");
	fprintf(stderr, "-w: write the graph to a file\n");
	fprintf(stderr, "-T: the number of threads to process in parallel\n");
}

int main(int argc, char *argv[])
{
	int opt;
	int num_threads = 1;
	bool directed = true;
	int num_opts = 0;
	char *type_str = NULL;
	bool merge_graph = false;
	bool write_graph = false;
	while ((opt = getopt(argc, argv, "uvt:mwT:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'u':
				directed = false;
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
			case 'T':
				num_threads = atoi(optarg);
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

	std::string adjacency_list_file = argv[0];
	adjacency_list_file += std::string("-v") + itoa(CURR_VERSION);
	work_dir = dirname(argv[0]);
	printf("work dir: %s\n", work_dir.c_str());

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
		edge_graph::ptr edge_g = parse_edge_lists(edge_list_files, edge_attr_type,
				directed, num_threads);
		disk_serial_graph::ptr g
			= std::static_pointer_cast<disk_serial_graph, serial_graph>(
					construct_graph(edge_g, work_dir, num_threads));
		// Write the constructed individual graph to a file.
		if (write_graph) {
			assert(!file_exist(adjacency_list_file));
			assert(!file_exist(index_file));
			g->dump(index_file, adjacency_list_file, true);
		}
		printf("There are %ld vertices, %ld non-empty vertices and %ld edges\n",
				g->get_num_vertices(), g->get_num_non_empty_vertices(),
				g->get_num_edges());
		if (check_graph)
			g->check_ext_graph(*edge_g, index_file, adjacency_list_file);
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

			edge_graph::ptr edge_g = parse_edge_lists(files, edge_attr_type,
					directed, num_threads);
			disk_serial_graph::ptr g
				= std::static_pointer_cast<disk_serial_graph, serial_graph>(
						construct_graph(edge_g, work_dir, num_threads));
			// Write the constructed individual graph to a file.
			if (write_graph) {
				assert(!file_exist(graph_files[i]));
				assert(!file_exist(index_files[i]));
				g->dump(index_files[i], graph_files[i], true);
			}
			printf("There are %ld vertices, %ld non-empty vertices and %ld edges\n",
					g->get_num_vertices(), g->get_num_non_empty_vertices(),
					g->get_num_edges());
			if (check_graph)
				g->check_ext_graph(*edge_g, index_files[i], graph_files[i]);
		}
	}
}
