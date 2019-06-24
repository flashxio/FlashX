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

#include "el2fg.h"

namespace {
void print_usage()
{
	fprintf(stderr, "convert an edge list to adjacency lists\n");
	fprintf(stderr, "el2fg [options] adj_list_file index_file "
            "edge_list_files (or directories)\n");
	fprintf(stderr, "-u: undirected graph\n");
	fprintf(stderr, "-v: verify the created adjacency list\n");
	fprintf(stderr, "-t type: the type of edge data. Supported type: ");
	for (int i = 0; i < fg::utils::type_map_size; i++) {
		fprintf(stderr, "%s, ", fg::utils::edge_type_map[i].str.c_str());
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "-m: merge multiple edge lists into a single graph. \n");
	fprintf(stderr, "-w: write the graph to a file\n");
	fprintf(stderr, "-T: the number of threads to process in parallel\n");
	fprintf(stderr, "-d: store intermediate data on disks\n");
	fprintf(stderr, "-c: the SAFS configuration file\n");
	fprintf(stderr, "-W: the working directory\n");
	fprintf(stderr, "-b: the size of the buffer for sorting edge lists\n");
	fprintf(stderr, "-B: the size of the buffer for writing the graph\n");
}
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
	bool on_disk = false;
	size_t sort_buf_size = 0;
	size_t write_buf_size = 0;
	std::string conf_file;
	std::string work_dir = ".";
	while ((opt = getopt(argc, argv, "uvt:mwT:dc:W:b:B:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'u':
				directed = false;
				break;
			case 'v':
				fg::utils::check_graph = true;
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
			case 'd':
				on_disk = true;
				break;
			case 'c':
				conf_file = optarg;
				num_opts++;
				break;
			case 'W':
				work_dir = optarg;
				num_opts++;
				break;
			case 'b':
				sort_buf_size = str2size(optarg);
				num_opts++;
				break;
			case 'B':
				write_buf_size = str2size(optarg);
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

	std::string adjacency_list_file = argv[0];
	std::string index_file = argv[1];

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

    fg::utils::el2fg(edge_list_files, adjacency_list_file, index_file, directed,
            num_threads, type_str, merge_graph, write_graph, on_disk,
            sort_buf_size, write_buf_size, conf_file, work_dir);
}
