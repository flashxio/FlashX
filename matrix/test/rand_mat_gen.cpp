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

#include "vv_store.h"
#include "fm_utils.h"
#include "mem_vv_store.h"
#include "sparse_matrix.h"

using namespace fm;

class rand_vertex: public set_vv_operate
{
	size_t max_neigh;
	size_t edge_data_size;
public:
	rand_vertex(size_t max_neigh, size_t edge_data_size) {
		this->max_neigh = max_neigh;
		this->edge_data_size = edge_data_size;
	}

	virtual void set(off_t arr_idx, void *arr, size_t num_eles) const;

	virtual const scalar_type &get_type() const {
		return get_scalar_type<char>();
	}
};

void rand_vertex::set(off_t arr_idx, void *arr, size_t num_eles) const
{
	size_t num_edges = fg::ext_mem_undirected_vertex::vsize2num_edges(
			num_eles, edge_data_size);
	new (arr) fg::ext_mem_undirected_vertex(arr_idx, num_edges, edge_data_size);
	fg::ext_mem_undirected_vertex *v = (fg::ext_mem_undirected_vertex *) arr;
	assert(v->get_size() == num_eles);
	for (size_t i = 0; i < num_edges; i++)
		v->set_neighbor(i, random() % max_neigh);
}

class sum_row: public arr_apply_operate
{
public:
	virtual void run(const local_vec_store &in, local_vec_store &out) const {
		const fg::ext_mem_undirected_vertex *v
			= (const fg::ext_mem_undirected_vertex *) in.get_raw_arr();
		assert(v->get_size() == in.get_length());
		fg::vertex_id_t sum = 0;
		for (size_t i = 0; i < v->get_num_edges(); i++)
			sum += v->get_neighbor(i);
		out.set<fg::vertex_id_t>(0, sum);
	}

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<fg::vertex_id_t>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<fg::vertex_id_t>();
	}
	virtual size_t get_num_out_eles(size_t num_input) const {
		return 1;
	}
};

void print_usage()
{
	fprintf(stderr, "rand_mat_gen nrow ncol nnz\n");
	fprintf(stderr, "-t type: the type of non-zero entries\n");
	fprintf(stderr, "-b size: the block size of the sparse matrix\n");
}

int main(int argc, char *argv[])
{
	int opt;
	int num_opts = 0;
	std::string ele_type = "binary";
	size_t block_size = 32 * 1024;
	size_t max_num_edges = 1000;
	bool verify = false;
	std::string output_name;
	while ((opt = getopt(argc, argv, "t:b:m:vo:")) != -1) {
		num_opts++;
		switch (opt) {
			case 't':
				ele_type = optarg;
				num_opts++;
				break;
			case 'b':
				block_size = str2size(optarg);
				num_opts++;
				break;
			case 'm':
				max_num_edges = str2size(optarg);
				num_opts++;
				break;
			case 'v':
				verify = true;
				break;
			case 'o':
				output_name = optarg;
				num_opts++;
				break;
			default:
				print_usage();
				exit(1);
		}
	}

	argv += 1 + num_opts;
	argc -= 1 + num_opts;
	if (argc < 3) {
		print_usage();
		return -1;
	}

	size_t nrow = atol(argv[0]);
	size_t ncol = atol(argv[1]);
	size_t nnz = atol(argv[2]);
	printf("create a matrix with %ld rows and %ld cols and %ld non-zero entries\n",
			nrow, ncol, nnz);
	const scalar_type *entry_type = NULL;
	if (ele_type == "integer")
		entry_type = &get_scalar_type<int>();
	else if (ele_type == "long")
		entry_type = &get_scalar_type<long>();
	else if (ele_type == "float")
		entry_type = &get_scalar_type<float>();
	else if (ele_type == "double")
		entry_type = &get_scalar_type<double>();
	size_t edge_data_size = 0;
	if (entry_type)
		edge_data_size = entry_type->get_size();
	max_num_edges = std::min(max_num_edges, ncol);

	// To generate a random matrix, we need to generate the FlashGraph format
	// of a random graph.

	// This indicates the offset of the adjacency list of each vertex.
	// The graph header isn't included.
	std::vector<off_t> offs(nrow + 1);
	offs[0] = 0;
	size_t remain_nnz = nnz;
	for (size_t i = 1; i < offs.size(); i++) {
		size_t max_num = std::min(max_num_edges, remain_nnz);
		size_t num_edges = random() % max_num;
		remain_nnz -= num_edges;
		offs[i] = offs[i - 1] + fg::ext_mem_undirected_vertex::num_edges2vsize(
				num_edges, edge_data_size);
	}
	printf("There remains %ld non-zero entries\n", remain_nnz);
	detail::mem_vv_store::ptr store = detail::mem_vv_store::create(offs,
			get_scalar_type<char>());
	store->set_data(rand_vertex(ncol, edge_data_size));

	block_2d_size bsize(block_size, block_size);
	vector_vector::ptr adj_vv = vector_vector::create(store);
	auto res = create_2d_matrix(adj_vv, ncol, bsize, entry_type);
	res.second->verify();
	if (!output_name.empty()) {
		res.first->dump(output_name + ".mat_idx");
		res.second->dump(output_name + ".mat");
	}

	if (verify) {
		res.second->verify();

		// Perform sparse matrix multiplication with a vector of sequence numbers.
		sparse_matrix::ptr mat = sparse_matrix::create(res.first, res.second);
		printf("The sparse matrix has %ld rows and %ld cols\n", mat->get_num_rows(),
				mat->get_num_cols());
		vector::ptr vec = create_seq_vector<fg::vertex_id_t>(0, ncol - 1, 1);
		detail::smp_vec_store::ptr spmv_tmp = detail::smp_vec_store::create(nrow,
				get_scalar_type<fg::vertex_id_t>());
		mat->multiply(detail::mem_vec_store::cast(
					vec->get_raw_store()), spmv_tmp);
		vector::ptr spmv_res = vector::create(spmv_tmp);

		// For a binary matrix, this operation is equivalent to sum all of
		// the column indexes in a row.
		vector_vector::ptr agg_vv = adj_vv->apply(sum_row());
		vector::ptr agg_res = agg_vv->cat();
		detail::smp_vec_store::const_ptr agg_tmp
			= detail::smp_vec_store::cast(agg_res->get_raw_store());
		assert(agg_res->equals(*spmv_res));
	}
}
