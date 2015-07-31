/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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

#include "safs_file.h"

#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "EM_dense_matrix.h"

#include "eigensolver.h"

using namespace fm;
using namespace fm::eigen;

class eigen_Operator: public spm_function
{
	sparse_matrix::ptr mat;
public:
	eigen_Operator(sparse_matrix::ptr mat) {
		this->mat = mat;
	}

	virtual dense_matrix::ptr run(dense_matrix::ptr x) const {
		assert(x->get_type() == get_scalar_type<double>());
		const detail::mem_matrix_store &mem_in
			= static_cast<const detail::mem_matrix_store &>(x->get_data());
		detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
				mat->get_num_rows(), mem_in.get_num_cols(),
				matrix_layout_t::L_COL, mem_in.get_type(), mem_in.get_num_nodes());
		mat->multiply<double>(mem_in, *res);
		return dense_matrix::create(res);
	}

	virtual size_t get_num_cols() const {
		return mat->get_num_cols();
	}

	virtual size_t get_num_rows() const {
		return mat->get_num_rows();
	}
};

class apply1_2: public bulk_uoperate
{
public:
	virtual void runA(size_t num_eles, const void *in_arr,
			void *out_arr) const {
		const double *t_in = (const double *) in_arr;
		double *t_out = (double *) out_arr;
		for (size_t i = 0; i < num_eles; i++)
			t_out[i] = 1 / std::sqrt(t_in[i]);
	}

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<double>();
	}

	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<double>();
	}
};

/*
 * This is used to compute the eigenpairs of normalized adjacency matrix.
 */
class NA_eigen_Operator: public spm_function
{
	sparse_matrix::ptr mat;
	vector::ptr deg_vec1_2;
public:
	NA_eigen_Operator(sparse_matrix::ptr mat) {
		this->mat = mat;
		// Get V = 1
		vector::ptr vec = create_vector<double>(mat->get_num_cols(), 1);
		// Get degree of each row.
		detail::mem_vec_store::ptr deg = detail::mem_vec_store::create(
				mat->get_num_rows(), matrix_conf.get_num_nodes(),
				get_scalar_type<double>());
		printf("len1: %ld, len2: %ld\n", vec->get_length(), deg->get_length());
		mat->multiply<double>(dynamic_cast<const detail::mem_vec_store &>(
					vec->get_data()), *deg);
		vec = vector::create(deg);
		// Get D^-1/2.
		deg_vec1_2 = vec->sapply(bulk_uoperate::const_ptr(new apply1_2()));
	}

	virtual dense_matrix::ptr run(dense_matrix::ptr x) const {
		assert(x->get_type() == get_scalar_type<double>());
		dense_matrix::ptr tmp = x->scale_rows(deg_vec1_2);
		x = NULL;
		tmp->materialize_self();
		const detail::mem_matrix_store &mem_in
			= static_cast<const detail::mem_matrix_store &>(tmp->get_data());
		detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
				mat->get_num_rows(), mem_in.get_num_cols(),
				matrix_layout_t::L_COL, mem_in.get_type(), mem_in.get_num_nodes());
		mat->multiply<double>(mem_in, *res);
		dense_matrix::ptr tmp2 = dense_matrix::create(res);
		return tmp2->scale_rows(deg_vec1_2);
	}

	virtual size_t get_num_cols() const {
		return mat->get_num_cols();
	}

	virtual size_t get_num_rows() const {
		return mat->get_num_rows();
	}
};

class SVD_Operator: public spm_function
{
	sparse_matrix::ptr mat;
public:
	SVD_Operator(sparse_matrix::ptr mat) {
		this->mat = mat;
	}

	virtual dense_matrix::ptr run(dense_matrix::ptr x) const {
		assert(x->get_type() == get_scalar_type<double>());
		const detail::mem_matrix_store &mem_in
			= static_cast<const detail::mem_matrix_store &>(x->get_data());
		detail::mem_matrix_store::ptr tmp = detail::mem_matrix_store::create(
				mat->get_num_rows(), mem_in.get_num_cols(),
				matrix_layout_t::L_ROW, mem_in.get_type(), mem_in.get_num_nodes());
		mat->multiply<double>(mem_in, *tmp);
		mat->transpose();

		detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
				mat->get_num_rows(), mem_in.get_num_cols(),
				matrix_layout_t::L_COL, mem_in.get_type(), mem_in.get_num_nodes());
		mat->multiply<double>(*tmp, *res);
		mat->transpose();
		return dense_matrix::create(res);
	}

	virtual size_t get_num_cols() const {
		return mat->get_num_cols();
	}

	virtual size_t get_num_rows() const {
		return mat->get_num_rows();
	}
};

void print_usage()
{
	fprintf(stderr, "eigensolver conf_file matrix_file index_file nev [options]\n");
	fprintf(stderr, "-b block_size\n");
	fprintf(stderr, "-n num_blocks\n");
	fprintf(stderr, "-s solver: Davidson, KrylovSchur, LOBPCG\n");
	fprintf(stderr, "-t tolerance\n");
	fprintf(stderr, "-e: the external memory mode.\n");
	fprintf(stderr, "-o file: output eigenvectors\n");
	fprintf(stderr, "-T type: eigen, SVD, NA_eigen (normalized adjacency)\n");
}

int main (int argc, char *argv[])
{
	int opt;
	int num_opts = 0;
	std::string output_file;
	struct eigen_options opts;
	std::string type = "eigen";
	while ((opt = getopt(argc, argv, "b:n:s:t:eo:T:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'b':
				opts.block_size = atoi(optarg);
				num_opts++;
				break;
			case 'n':
				opts.num_blocks = atoi(optarg);
				num_opts++;
				break;
			case 's':
				opts.solver = optarg;
				num_opts++;
				break;
			case 't':
				opts.tol = atof(optarg);
				num_opts++;
				break;
			case 'e':
				opts.in_mem = false;
				break;
			case 'o':
				output_file = optarg;
				num_opts++;
				break;
			case 'T':
				type = optarg;
				num_opts++;
				break;
			default:
				print_usage();
				abort();
		}
	}

	argv += 1 + num_opts;
	argc -= 1 + num_opts;
	if (argc < 4 || (type == "SVD" && argc < 6)) {
		print_usage();
		exit(1);
	}

	std::string conf_file = argv[0];
	std::string matrix_file = argv[1];
	std::string index_file = argv[2];
	std::string t_matrix_file;
	std::string t_index_file;
	if (type == "SVD") {
		t_matrix_file = argv[3];
		t_index_file = argv[4];
		opts.nev = atoi(argv[5]);
	}
	else
		opts.nev = atoi(argv[3]); // number of eigenvalues for which to solve;

	//
	// Set up the test problem.
	//
	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	// Load index.
	SpM_2d_index::ptr index;
	safs::safs_file idx_f(safs::get_sys_RAID_conf(), index_file);
	if (idx_f.exist())
		index = SpM_2d_index::safs_load(index_file);
	else
		index = SpM_2d_index::load(index_file);
	SpM_2d_index::ptr t_index;
	if (type == "SVD") {
		safs::safs_file t_idx_f(safs::get_sys_RAID_conf(), t_index_file);
		if (t_idx_f.exist())
			t_index = SpM_2d_index::safs_load(t_index_file);
		else
			t_index = SpM_2d_index::load(t_index_file);
	}

	// Load matrix.
	sparse_matrix::ptr mat;
	safs::safs_file mat_f(safs::get_sys_RAID_conf(), matrix_file);
	if (type == "SVD" && mat_f.exist())
		mat = sparse_matrix::create(
				index, safs::create_io_factory(matrix_file, safs::REMOTE_ACCESS),
				t_index, safs::create_io_factory(t_matrix_file,
					safs::REMOTE_ACCESS));
	else if (type == "SVD")
		mat = sparse_matrix::create(
				index, SpM_2d_storage::load(matrix_file, index),
				t_index, SpM_2d_storage::load(t_matrix_file, t_index));
	else if (mat_f.exist())
		mat = sparse_matrix::create(index,
				safs::create_io_factory(matrix_file, safs::REMOTE_ACCESS));
	else
		mat = sparse_matrix::create(index,
				SpM_2d_storage::load(matrix_file, index));

	eigen_res res;
	if (type == "SVD")
		res = compute_eigen(new SVD_Operator(mat), true, opts);
	else if (type == "NA_eigen")
		res = compute_eigen(new NA_eigen_Operator(mat), mat->is_symmetric(), opts);
	else
		res = compute_eigen(new eigen_Operator(mat), mat->is_symmetric(), opts);
	// We only save eigenvectors if they are stored in memory.
	if (!output_file.empty()) {
		printf("Save eigenvectors to %s\n", output_file.c_str());
		if (res.vecs->is_in_mem())
			res.vecs = res.vecs->conv_store(false, -1);
		else
			res.vecs->materialize_self();
		const detail::EM_matrix_store &store
			= dynamic_cast<const detail::EM_matrix_store &>(res.vecs->get_data());
		bool ret = store.set_persistent(output_file);
		assert(ret);
	}

	safs::print_io_summary();

	return 0;
}
