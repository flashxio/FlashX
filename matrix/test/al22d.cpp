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

using namespace fm;

void verify_2d_matrix(const std::string &mat_file, const std::string &mat_idx_file)
{
	SpM_2d_index::ptr idx = SpM_2d_index::load(mat_idx_file);
	SpM_2d_storage::ptr mat = SpM_2d_storage::load(mat_file, idx);
	mat->verify();
}

size_t get_out_size(fg::vertex_index::ptr vindex)
{
	if (vindex->is_compressed()) {
		fg::vsize_t num_vertices = vindex->get_num_vertices();
		if (vindex->get_graph_header().is_directed_graph()) {
			fg::in_mem_cdirected_vertex_index::ptr dindex
				= fg::in_mem_cdirected_vertex_index::create(*vindex);
			fg::directed_vertex_entry dentry = dindex->get_vertex(num_vertices - 1);
			return dentry.get_out_off() + dindex->get_out_size(
					num_vertices - 1) - vindex->get_out_part_loc();
		}
		else {
			fg::in_mem_cundirected_vertex_index::ptr uindex
				= fg::in_mem_cundirected_vertex_index::create(*vindex);
			fg::vertex_offset off = uindex->get_vertex(num_vertices - 1);
			return off.get_off() + uindex->get_size(
					num_vertices - 1) - vindex->get_header_size();
		}
	}
	else {
		if (vindex->get_graph_header().is_directed_graph()) {
			fg::directed_vertex_index::ptr dindex
				= fg::directed_vertex_index::cast(vindex);
			return dindex->get_graph_size() - vindex->get_out_part_loc();
		}
		else {
			fg::undirected_vertex_index::ptr uindex
				= fg::undirected_vertex_index::cast(vindex);
			return uindex->get_graph_size() - vindex->get_header_size();
		}
	}
}

size_t get_in_size(fg::vertex_index::ptr vindex)
{
	assert(vindex->get_graph_header().is_directed_graph());
	return vindex->get_out_part_loc() - vindex->get_header_size();
}

size_t get_out_off(fg::vertex_index::ptr vindex)
{
	if (vindex->get_graph_header().is_directed_graph())
		return vindex->get_out_part_loc();
	else
		return vindex->get_header_size();
}

size_t get_in_off(fg::vertex_index::ptr vindex)
{
	assert(vindex->get_graph_header().is_directed_graph());
	return vindex->get_header_size();
}

void init_out_offs(fg::vertex_index::ptr vindex, std::vector<off_t> &out_offs)
{
	size_t num_vertices = vindex->get_num_vertices();
	assert(num_vertices + 1 == out_offs.size());
	if (vindex->is_compressed()) {
		if (vindex->get_graph_header().is_directed_graph()) {
			fg::in_mem_cdirected_vertex_index::ptr dindex
				= fg::in_mem_cdirected_vertex_index::create(*vindex);
			out_offs[0] = 0;
			for (size_t i = 1; i <= num_vertices; i++)
				out_offs[i] = out_offs[i - 1] + dindex->get_out_size(i - 1);
		}
		else {
			fg::in_mem_cundirected_vertex_index::ptr uindex
				= fg::in_mem_cundirected_vertex_index::create(*vindex);
			out_offs[0] = 0;
			for (size_t i = 1; i <= num_vertices; i++)
				out_offs[i] = out_offs[i - 1] + uindex->get_size(i - 1);
		}
		assert((size_t) out_offs[num_vertices] == get_out_size(vindex));
	}
	else {
		if (vindex->get_graph_header().is_directed_graph()) {
			off_t out_part_loc = vindex->get_out_part_loc();
			fg::directed_vertex_index::ptr dindex
				= fg::directed_vertex_index::cast(vindex);
			for (size_t i = 0; i < num_vertices; i++)
				out_offs[i] = dindex->get_vertex(i).get_out_off()
					- out_part_loc;
			out_offs[num_vertices] = get_out_size(vindex);
		}
		else {
			fg::undirected_vertex_index::ptr uindex
				= fg::undirected_vertex_index::cast(vindex);
			for (size_t i = 0; i < num_vertices; i++)
				out_offs[i] = uindex->get_vertex(i).get_off()
					- vindex->get_header_size();
			out_offs[num_vertices] = get_out_size(vindex);
		}
	}
	for (size_t i = 1; i <= num_vertices; i++)
		assert(out_offs[i] > out_offs[i - 1]);
}

void init_in_offs(fg::vertex_index::ptr vindex, std::vector<off_t> &in_offs)
{
	size_t num_vertices = vindex->get_num_vertices();
	assert(num_vertices + 1 == in_offs.size());
	assert(vindex->get_graph_header().is_directed_graph());
	if (vindex->is_compressed()) {
		fg::in_mem_cdirected_vertex_index::ptr dindex
			= fg::in_mem_cdirected_vertex_index::create(*vindex);
		in_offs[0] = 0;
		for (size_t i = 1; i <= num_vertices; i++)
			in_offs[i] = in_offs[i - 1] + dindex->get_in_size(i - 1);
		assert((size_t) in_offs[num_vertices] == get_in_size(vindex));
	}
	else {
		fg::directed_vertex_index::ptr dindex
			= fg::directed_vertex_index::cast(vindex);
		for (size_t i = 0; i < num_vertices; i++)
			in_offs[i] = dindex->get_vertex(i).get_in_off()
				- vindex->get_header_size();
		in_offs[num_vertices] = get_in_size(vindex);
	}
	for (size_t i = 1; i <= num_vertices; i++)
		assert(in_offs[i] > in_offs[i - 1]);
}

void read_data(char *data, size_t size, size_t num, off_t off, FILE *stream)
{
	int iret = fseek(stream, off, SEEK_SET);
	assert(iret == 0);
	size_t ret = fread(data, size, num, stream);
	assert(ret == num);
}

struct deleter
{
	void operator()(char *buf) {
		free(buf);
	}
};

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "el2al graph_file index_file matrix_name\n");
		return -1;
	}

	std::string graph_file = std::string(argv[1]);
	std::string index_file = std::string(argv[2]);
	std::string mat_name = std::string(argv[3]);

	std::string mat_file = mat_name + ".mat";
	std::string mat_idx_file = mat_name + ".mat_idx";

	size_t block_height
		= ((size_t) std::numeric_limits<unsigned short>::max()) + 1;
	block_2d_size block_size(block_height, block_height);
	fg::vertex_index::ptr vindex = fg::vertex_index::load(index_file);

	FILE *f = fopen(graph_file.c_str(), "r");
	if (f == NULL) {
		fprintf(stderr, "can't open %s: %s\n", graph_file.c_str(),
				strerror(errno));
		return -1;
	}

	size_t out_size = get_out_size(vindex);
	off_t out_off = get_out_off(vindex);
	std::shared_ptr<char> out_data = std::shared_ptr<char>(
			(char *) malloc(out_size), deleter());
	read_data(out_data.get(), out_size, 1, out_off, f);
	std::vector<off_t> out_offs(vindex->get_num_vertices() + 1);
	init_out_offs(vindex, out_offs);
	mem_vector_vector::ptr out_adjs = type_mem_vector_vector<char>::create(
			out_data, out_size, out_offs);

	// Construct 2D partitioning of the adjacency matrix.
	export_2d_matrix(out_adjs, block_size, mat_file, mat_idx_file);
	verify_2d_matrix(mat_file, mat_idx_file);
	out_adjs = NULL;

	if (vindex->get_graph_header().is_directed_graph()) {
		std::string t_mat_file = mat_name + "_t.mat";
		std::string t_mat_idx_file = mat_name + "_t.mat_idx";

		size_t in_size = get_in_size(vindex);
		off_t in_off = get_in_off(vindex);
		std::shared_ptr<char> in_data = std::shared_ptr<char>(
				(char *) malloc(in_size), deleter());
		read_data(in_data.get(), in_size, 1, in_off, f);
		std::vector<off_t> in_offs(vindex->get_num_vertices() + 1);
		init_in_offs(vindex, in_offs);
		mem_vector_vector::ptr in_adjs = type_mem_vector_vector<char>::create(
				in_data, in_size, in_offs);

		// Construct 2D partitioning of the adjacency matrix.
		export_2d_matrix(in_adjs, block_size, t_mat_file, t_mat_idx_file);
		verify_2d_matrix(t_mat_file, t_mat_idx_file);
	}

	fclose(f);
}
