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

#include "eigensolver.h"

using namespace fm;

class FM_Operator: public spm_function
{
	sparse_matrix::ptr mat;
public:
	FM_Operator(sparse_matrix::ptr mat) {
		this->mat = mat;
	}

	virtual void Apply(const block_multi_vector& x,
			block_multi_vector& y) const {
		block_multi_vector::sparse_matrix_multiply<double>(*mat, x, y);

#if 0
		assert((size_t) x.GetGlobalLength() == mat->get_num_cols());
		assert((size_t) y.GetGlobalLength() == mat->get_num_rows());
		mem_vector::ptr in = mem_vector::create(mat->get_num_cols(),
				get_scalar_type<double>());
		mem_vector::ptr out = mem_vector::create(mat->get_num_rows(),
				get_scalar_type<double>());
		for (int i = 0; i < x.GetNumberVecs(); i++) {
			memcpy(in->get_raw_arr(), x.get_ep_mv()[i], in->get_length() * sizeof(double));
			out->reset_data();
			mat->multiply<double>(*in, *out);
			memcpy(y.get_ep_mv()[i], out->get_raw_arr(), out->get_length() * sizeof(double));
		}
		y.sync_ep2fm();
#endif
	}

	virtual size_t get_num_cols() const {
		return mat->get_num_cols();
	}

	virtual size_t get_num_rows() const {
		return mat->get_num_rows();
	}
};

int main (int argc, char *argv[])
{
	if (argc < 5) {
		fprintf(stderr, "eigensolver conf_file matrix_file index_file nev [solver]\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	std::string matrix_file = argv[2];
	std::string index_file = argv[3];

	struct eigen_options opts;
	opts.nev = atoi(argv[4]); // number of eigenvalues for which to solve;
	if (argc >= 6)
		opts.solver = argv[5];
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

	// Load matrix.
	sparse_matrix::ptr mat;
	safs::safs_file mat_f(safs::get_sys_RAID_conf(), matrix_file);
	if (mat_f.exist())
		mat = sparse_matrix::create(index, safs::create_io_factory(
					matrix_file, safs::REMOTE_ACCESS));
	else
		mat = sparse_matrix::create(index,
				SpM_2d_storage::load(matrix_file, index));

	eigen_res res = compute_eigen(new FM_Operator(mat),
			mat->is_symmetric(), opts);

	return 0;
}
