#ifndef __EIGENSOLVER_H__
#define __EIGENSOLVER_H__

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
#include <vector>
#include <string>
#include <memory>

#include "dense_matrix.h"
#include "matrix_store.h"

namespace fm
{

namespace eigen
{

class spm_function
{
public:
	typedef std::shared_ptr<const spm_function> const_ptr;

	/*
	 * We pass a reference of the shared pointer to the run method.
	 * So the user-defined run method can deallocate the input dense matrix.
	 */
	virtual fm::dense_matrix::ptr run(fm::dense_matrix::ptr &x) const = 0;
	virtual size_t get_num_cols() const = 0;
	virtual size_t get_num_rows() const = 0;
};

struct eigen_res
{
	std::vector<double> vals;
	fm::dense_matrix::ptr vecs;
};

struct eigen_options
{
	double tol;
	int num_blocks;
	int max_restarts;
	int max_iters;
	int block_size;
	int nev;
	std::string solver;
	std::string which;
	bool in_mem;

	eigen_options(int nev, std::string solver = "KrylovSchur");
};

/*
 * `func' will be destroyed by this function.
 */
eigen_res compute_eigen(spm_function *func, bool sym,
		struct eigen_options &opts);

}

}

#endif
