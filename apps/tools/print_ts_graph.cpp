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

/**
 * This is a test file to print a small time-series graph.
 */

#include <stdio.h>

#include <string>

#include "native_file.h"
#include "vertex_index.h"

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "print_ts_graph adj_list_file index_file\n");
		return -1;
	}

	const std::string adj_file_name = argv[1];
	const std::string index_file_name = argv[2];

	native_file adj_file(adj_file_name);
	ssize_t adj_file_size = adj_file.get_size();
	char *adj_list = new char[adj_file_size];
	FILE *f = fopen(adj_file_name.c_str(), "r");
	if (f == NULL) {
		perror("fopen");
		return -1;
	}
	size_t ret = fread(adj_list, adj_file_size, 1, f);
	assert(ret == 1);

	default_vertex_index *index = default_vertex_index::load(index_file_name);
	int num_edges = 0;
	int num_vertices = 0;
	for (size_t i = 0; i < index->get_num_vertices(); i++) {
		off_t off = index->get_vertex_off(i);
		int size = index->get_vertex_size(i);
		ts_ext_mem_directed_vertex *v
			= (ts_ext_mem_directed_vertex *) (adj_list + off);
		if (v->get_num_edges() > 0) {
			v->print();
			num_vertices++;
			num_edges += v->get_num_edges();
		}
		assert(size == (int) v->get_size());
		assert(i == v->get_id());
	}
	printf("There are %d vertices and %d edges\n", num_vertices, num_edges);
}
