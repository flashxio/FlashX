#ifndef __MATRIX_EIGENSOLVER_H__
#define __MATRIX_EIGENSOLVER_H__

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

#include "FG_vector.h"

typedef double ev_float_t;
typedef std::pair<ev_float_t, FG_vector<ev_float_t>::ptr> eigen_pair_t;

class SPMV
{
public:
	virtual void compute(const FG_vector<ev_float_t> &input,
			FG_vector<ev_float_t> &output) = 0;
	virtual size_t get_vector_size() = 0;
};

template<class SparseMatrixType>
class eigen_SPMV: public SPMV
{
	typename SparseMatrixType::ptr A;
public:
	eigen_SPMV(typename SparseMatrixType::ptr A) {
		this->A = A;
		assert(A->get_num_rows() == A->get_num_cols());
	}

	virtual void compute(const FG_vector<ev_float_t> &input,
			FG_vector<ev_float_t> &output) {
		A->multiply(input, output);
	}

	virtual size_t get_vector_size() {
		return A->get_num_rows();
	}
};

template<class SparseMatrixType>
class LS_SPMV: public SPMV
{
	typename SparseMatrixType::ptr A;
	FG_vector<ev_float_t>::ptr tmp;
public:
	LS_SPMV(typename SparseMatrixType::ptr A) {
		this->A = A;
		tmp = FG_vector<ev_float_t>::create(A->get_num_cols());
	}

	virtual void compute(const FG_vector<ev_float_t> &input,
			FG_vector<ev_float_t> &output) {
		A->transpose()->multiply(input, *tmp);
		A->multiply(*tmp, output);
	}

	virtual size_t get_vector_size() {
		return A->get_num_rows();
	}
};

template<class SparseMatrixType>
class RS_SPMV: public SPMV
{
	typename SparseMatrixType::ptr A;
	FG_vector<ev_float_t>::ptr tmp;
public:
	RS_SPMV(typename SparseMatrixType::ptr A) {
		this->A = A;
		tmp = FG_vector<ev_float_t>::create(A->get_num_rows());
	}

	virtual void compute(const FG_vector<ev_float_t> &input,
			FG_vector<ev_float_t> &output) {
		A->multiply(input, *tmp);
		A->transpose()->multiply(*tmp, output);
	}

	virtual size_t get_vector_size() {
		return A->get_num_cols();
	}
};

void eigen_solver(SPMV &spmv, int m, int nv, const std::string &which,
		std::vector<eigen_pair_t> &eigen_pairs);

template<class SparseMatrixType>
void compute_eigen(typename SparseMatrixType::ptr matrix, int m, int nv,
		const std::string &which, std::vector<eigen_pair_t> &eigen_pairs)
{
	eigen_SPMV<SparseMatrixType> spmv(matrix);
	eigen_solver(spmv, m, nv, which, eigen_pairs);
}

template<class SparseMatrixType>
void compute_SVD(typename SparseMatrixType::ptr matrix, int m, int nv,
		const std::string &which, const std::string &type,
		std::vector<eigen_pair_t> &eigen_pairs)
{
	// left-singular vectors of SVD
	if (type == "LS") {
		LS_SPMV<SparseMatrixType> spmv(matrix);
		eigen_solver(spmv, m, nv, which, eigen_pairs);
	}
	// right-singular vectors of SVD
	else if (type == "RS") {
		RS_SPMV<SparseMatrixType> spmv(matrix);
		eigen_solver(spmv, m, nv, which, eigen_pairs);
	}
	else
		assert(0);

	for (int i = 0; i < nv; i++) {
		ev_float_t v = eigen_pairs[i].first;
		eigen_pairs[i].first = sqrt(v);
	}
}

#endif
