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
#include "fg_utils.h"

#include "vertex_index.h"

void read_data(char *data, size_t size, size_t num, off_t off, FILE *stream)
{
	int iret = fseek(stream, off, SEEK_SET);
	assert(iret == 0);
	size_t ret = fread(data, size, num, stream);
	assert(ret == num);
}

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "el2crs graph_file index_file crs_file\n");
		return -1;
	}

	std::string graph_file = std::string(argv[1]);
	std::string index_file = std::string(argv[2]);
	std::string crs_file = std::string(argv[3]);

	fg::vertex_index::ptr vindex = fg::vertex_index::load(index_file);

	FILE *f = fopen(graph_file.c_str(), "r");
	if (f == NULL) {
		fprintf(stderr, "can't open %s: %s\n", graph_file.c_str(),
				strerror(errno));
		return -1;
	}

#if 0
    // TODO
	size_t out_size = fg::get_out_size(vindex);
	off_t out_off = fg::get_out_off(vindex);
	detail::simple_raw_array out_data(out_size, -1);
	read_data(out_data.get_raw(), out_size, 1, out_off, f);
	std::vector<off_t> out_offs(vindex->get_num_vertices() + 1);
	fg::init_out_offs(vindex, out_offs);
	for (size_t i = 0; i < out_offs.size(); i++)
		out_offs[i] -= out_off;
	vector_vector::ptr out_adjs = vector_vector::create(
			out_data, out_offs, get_scalar_type<char>());

	export_crs(out_adjs, crs_file);
	out_adjs = NULL;

	if (vindex->get_graph_header().is_directed_graph()) {
		std::string t_crs_file = std::string("t_") + crs_file;

		size_t in_size = fg::get_in_size(vindex);
		off_t in_off = fg::get_in_off(vindex);
		detail::simple_raw_array in_data(in_size, -1);
		read_data(in_data.get_raw(), in_size, 1, in_off, f);
		std::vector<off_t> in_offs(vindex->get_num_vertices() + 1);
		fg::init_in_offs(vindex, in_offs);
		for (size_t i = 0; i < in_offs.size(); i++)
			in_offs[i] -= in_off;
		vector_vector::ptr in_adjs = vector_vector::create(
				in_data, in_offs, get_scalar_type<char>());

		export_crs(in_adjs, t_crs_file);
	}
	fclose(f);
#endif
}
