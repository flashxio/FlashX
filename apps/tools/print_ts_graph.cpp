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
