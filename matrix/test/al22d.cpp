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

#include "vertex_index.h"

#include "mem_vector_vector.h"
#include "fm_utils.h"
#include "fg_utils.h"

using namespace fm;

void verify_2d_matrix(const std::string &mat_file, const std::string &mat_idx_file)
{
	SpM_2d_index::ptr idx = SpM_2d_index::load(mat_idx_file);
	SpM_2d_storage::ptr mat = SpM_2d_storage::load(mat_file, idx);
	mat->verify();
}

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
		fprintf(stderr, "el2al graph_file index_file matrix_name [block_height [block_width]]\n");
		return -1;
	}

	std::string graph_file = std::string(argv[1]);
	std::string index_file = std::string(argv[2]);
	std::string mat_name = std::string(argv[3]);
	size_t block_height = block_max_num_rows;
	size_t block_width = block_height;
	if (argc >= 5)
		block_width = block_height = std::atoi(argv[4]);
	if (argc >= 6)
		block_width = std::atoi(argv[5]);

	std::string mat_file = mat_name + ".mat";
	std::string mat_idx_file = mat_name + ".mat_idx";

	block_2d_size block_size(block_height, block_width);
	fg::vertex_index::ptr vindex = fg::vertex_index::load(index_file);

	FILE *f = fopen(graph_file.c_str(), "r");
	if (f == NULL) {
		fprintf(stderr, "can't open %s: %s\n", graph_file.c_str(),
				strerror(errno));
		return -1;
	}

	size_t out_size = get_out_size(vindex);
	off_t out_off = get_out_off(vindex);
	detail::raw_data_array out_data(out_size);
	read_data(out_data.get_raw(), out_size, 1, out_off, f);
	std::vector<off_t> out_offs(vindex->get_num_vertices() + 1);
	init_out_offs(vindex, out_offs);
	mem_vector_vector::ptr out_adjs = mem_vector_vector::create(
			out_data, out_offs, get_scalar_type<char>());

	// Construct 2D partitioning of the adjacency matrix.
	export_2d_matrix(out_adjs, block_size, mat_file, mat_idx_file);
	verify_2d_matrix(mat_file, mat_idx_file);
	out_adjs = NULL;

	if (vindex->get_graph_header().is_directed_graph()) {
		std::string t_mat_file = mat_name + "_t.mat";
		std::string t_mat_idx_file = mat_name + "_t.mat_idx";

		size_t in_size = get_in_size(vindex);
		off_t in_off = get_in_off(vindex);
		detail::raw_data_array in_data(in_size);
		read_data(in_data.get_raw(), in_size, 1, in_off, f);
		std::vector<off_t> in_offs(vindex->get_num_vertices() + 1);
		init_in_offs(vindex, in_offs);
		mem_vector_vector::ptr in_adjs = mem_vector_vector::create(
				in_data, in_offs, get_scalar_type<char>());

		// Construct 2D partitioning of the adjacency matrix.
		export_2d_matrix(in_adjs, block_size, t_mat_file, t_mat_idx_file);
		verify_2d_matrix(t_mat_file, t_mat_idx_file);
	}

	fclose(f);
}
