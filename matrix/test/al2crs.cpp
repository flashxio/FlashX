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
#include "crs_header.h"
#include "local_vec_store.h"

using namespace fm;

void read_data(char *data, size_t size, size_t num, off_t off, FILE *stream)
{
	int iret = fseek(stream, off, SEEK_SET);
	assert(iret == 0);
	size_t ret = fread(data, size, num, stream);
	assert(ret == num);
}

class adj2crs_apply_operate: public arr_apply_operate
{
public:
	adj2crs_apply_operate(): arr_apply_operate(0) {
	}
	virtual void run(const local_vec_store &in,
			local_vec_store &out) const;

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<char>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<crs_idx_t>();
	}
};

void adj2crs_apply_operate::run(const local_vec_store &in,
		local_vec_store &out) const
{
	const fg::ext_mem_undirected_vertex *v
		= (const fg::ext_mem_undirected_vertex *) in.get_raw_arr();
	size_t num_edges = v->get_num_edges();
	out.resize(num_edges);
	for (size_t i = 0; i < num_edges; i++)
		out.set<crs_idx_t>(i, v->get_neighbor(i));
}

void export_crs(mem_vector_vector::ptr adjs, const std::string &output_file)
{
	mem_vector_vector::ptr col_idxs = mem_vector_vector::cast(
			adjs->apply(adj2crs_apply_operate()));
	mem_vector::ptr col_vec = mem_vector::cast(col_idxs->cat());
	std::vector<crs_idx_t> offs(adjs->get_num_vecs() + 1);
	for (size_t i = 0; i < adjs->get_num_vecs(); i++) {
		offs[i + 1] = offs[i] + col_idxs->get_length(i);
	}

	FILE *f = fopen(output_file.c_str(), "w");
	if (f == NULL) {
		fprintf(stderr, "can't open %s: %s\n", output_file.c_str(), strerror(errno));
		exit(1);
	}

	// This is an adjacency matrix of a graph, so it has the same number of
	// rows and columns.
	crs_header header(adjs->get_num_vecs(), adjs->get_num_vecs(),
			col_vec->get_length());
	size_t ret = fwrite(&header, sizeof(header), 1, f);
	if (ret == 0) {
		fprintf(stderr, "can't write header: %s\n", strerror(errno));
		exit(1);
	}

	ret = fwrite(offs.data(), offs.size() * sizeof(offs[0]), 1, f);
	if (ret == 0) {
		fprintf(stderr, "can't write the row pointers: %s\n", strerror(errno));
		exit(1);
	}
	ret = fwrite(col_vec->get_raw_arr(),
			col_vec->get_length() * col_vec->get_entry_size(), 1, f);
	if (ret == 0) {
		fprintf(stderr, "can't write col idx: %s\n", strerror(errno));
		exit(1);
	}

	fclose(f);
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

	size_t out_size = get_out_size(vindex);
	off_t out_off = get_out_off(vindex);
	detail::raw_data_array out_data(out_size);
	read_data(out_data.get_raw(), out_size, 1, out_off, f);
	std::vector<off_t> out_offs(vindex->get_num_vertices() + 1);
	init_out_offs(vindex, out_offs);
	mem_vector_vector::ptr out_adjs = mem_vector_vector::create(
			out_data, out_offs, get_scalar_type<char>());

	export_crs(out_adjs, crs_file);
	out_adjs = NULL;

	if (vindex->get_graph_header().is_directed_graph()) {
		std::string t_crs_file = std::string("t_") + crs_file;

		size_t in_size = get_in_size(vindex);
		off_t in_off = get_in_off(vindex);
		detail::raw_data_array in_data(in_size);
		read_data(in_data.get_raw(), in_size, 1, in_off, f);
		std::vector<off_t> in_offs(vindex->get_num_vertices() + 1);
		init_in_offs(vindex, in_offs);
		mem_vector_vector::ptr in_adjs = mem_vector_vector::create(
				in_data, in_offs, get_scalar_type<char>());

		export_crs(in_adjs, t_crs_file);
	}

	fclose(f);
}
