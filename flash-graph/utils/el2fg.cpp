/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
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

#include <string.h>
#include <assert.h>

#include <vector>
#include <string>

#include "common.h"

#include "native_file.h"
#include "log.h"
#include "FG_basic_types.h"
#include "vertex.h"
#include "graph_file_header.h"
#include "vertex_index.h"
#include "in_mem_storage.h"
#include "FGlib.h"
#include "utils.h"
#include "fg_utils.h"

#include "generic_type.h"
#include "data_io.h"
#include "data_frame.h"
#include "vector.h"
#include "vector_vector.h"
#include "matrix_config.h"
#include "sparse_matrix.h"

using namespace fm;

void print_usage()
{
	fprintf(stderr, "convert an edge list to adjacency lists\n");
	fprintf(stderr, "el2al [options] conf_file edge_file graph_name\n");
	fprintf(stderr, "-u: undirected graph\n");
	fprintf(stderr, "-U: unqiue edges\n");
	fprintf(stderr, "-e: use external memory\n");
	fprintf(stderr, "-s size: sort buffer size\n");
	fprintf(stderr, "-g size: groupby buffer size\n");
	fprintf(stderr, "-t type: the edge attribute type\n");
	fprintf(stderr, "-d delim: specified the string as delimiter\n");
}

int main(int argc, char *argv[])
{
	bool directed = true;
	bool in_mem = true;
	bool uniq_edge = false;
	size_t sort_buf_size = 1UL * 1024 * 1024 * 1024;
	size_t groupby_buf_size = 1UL * 1024 * 1024 * 1024;
	int opt;
	int num_opts = 0;
	std::string edge_attr_type;
	std::string delim = "auto";
	while ((opt = getopt(argc, argv, "uUes:g:t:d:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'u':
				directed = false;
				break;
			case 'U':
				uniq_edge = true;
				break;
			case 'e':
				in_mem = false;
				break;
			case 's':
				sort_buf_size = str2size(optarg);
				num_opts++;
				break;
			case 'g':
				groupby_buf_size = str2size(optarg);
				num_opts++;
				break;
			case 't':
				edge_attr_type = optarg;
				num_opts++;
				break;
			case 'd':
				delim = optarg;
				num_opts++;
				break;
			default:
				print_usage();
				exit(1);
		}
	}

	argv += 1 + num_opts;
	argc -= 1 + num_opts;
	if (argc < 2) {
		print_usage();
		exit(1);
	}

	std::string conf_file = argv[0];
	std::string file_name = argv[1];
	std::string graph_name = argv[2];
	std::string adj_file = graph_name + ".adj";
	std::string index_file = graph_name + ".index";

	std::vector<std::string> files;
	safs::native_file f(file_name);
	if (f.exist() && !f.is_dir())
		files.push_back(file_name);
	else if (f.exist() && f.is_dir()) {
		safs::native_dir d(file_name);
		d.read_all_files(files);
		for (size_t i = 0; i < files.size(); i++)
			files[i] = file_name + "/" + files[i];
	}
	else {
		fprintf(stderr, "The input file %s doesn't exist\n", file_name.c_str());
		return -1;
	}
	printf("read edges from %ld files\n", files.size());

	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);
	matrix_conf.set_sort_buf_size(sort_buf_size);
	matrix_conf.set_groupby_buf_size(groupby_buf_size);
	printf("sort buf size: %ld, groupby buf size: %ld\n",
			matrix_conf.get_sort_buf_size(), matrix_conf.get_groupby_buf_size());
	fg::set_deduplicate(uniq_edge);

	{
		struct timeval start, end;
		/*
		 * We only need to indicate here whether we use external memory or not.
		 * If the columns in the data frame are in external memory, it'll always
		 * use external data containers for the remaining of processing.
		 */
		printf("start to read and parse edge list\n");
		gettimeofday(&start, NULL);
		data_frame::ptr df = fg::utils::read_edge_list(files, in_mem, delim,
				edge_attr_type, directed);
		assert(df);
		gettimeofday(&end, NULL);
		printf("It takes %.3f seconds to parse the edge lists\n",
				time_diff(start, end));
		printf("There are %ld edges\n", df->get_num_entries());
		assert(df->get_num_entries() > 0);

		fg::edge_list::ptr el = fg::edge_list::create(df, directed);
		printf("start to construct FlashGraph graph\n");
		fg::FG_graph::ptr graph = create_fg_graph(graph_name, el);

		if (graph && graph->get_index_data())
			graph->get_index_data()->dump(index_file);
		if (graph && graph->get_graph_data())
			graph->get_graph_data()->dump(adj_file);
	}
	destroy_flash_matrix();

	return 0;
}
