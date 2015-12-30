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

#include "matrix_algs.h"

namespace fm
{

namespace alg
{

size_t get_nnz(sparse_matrix::ptr mat)
{
	int num_nodes = matrix_conf.get_num_nodes();
	detail::matrix_store::ptr out_store = detail::matrix_store::create(
			mat->get_num_rows(), 1, matrix_layout_t::L_ROW,
			get_scalar_type<size_t>(), num_nodes, true);
	dense_matrix::ptr one = dense_matrix::create_const<size_t>(1,
			mat->get_num_cols(), 1, matrix_layout_t::L_ROW, num_nodes, true);
	mat->multiply<size_t, size_t>(one->get_raw_store(), out_store);
	dense_matrix::ptr out_deg = dense_matrix::create(out_store);
	scalar_variable::ptr v = out_deg->sum();
	return *(const size_t *) v->get_raw();
}

static dense_matrix::ptr multiply(sparse_matrix::ptr S, dense_matrix::ptr D)
{
	detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
			S->get_num_rows(), D->get_num_cols(), matrix_layout_t::L_ROW,
			D->get_type(), D->get_raw_store()->get_num_nodes());
	S->multiply<double, double>(D->get_raw_store(), res);
	return dense_matrix::create(res);
}

static std::pair<dense_matrix::ptr, dense_matrix::ptr> update_lee(
		sparse_matrix::ptr mat, dense_matrix::ptr W, dense_matrix::ptr H)
{
	double eps = 10e-9;
	// den <- (t(W) %*% W) %*% H
	dense_matrix::ptr D;
	{
		// t(W)
		dense_matrix::ptr tW = W->transpose();
		// t(W) %*% W
		dense_matrix::ptr tmp1 = tW->multiply(*W);
		// t(H)
		dense_matrix::ptr tH = H->transpose();
		// t(H) %*% (t(W) %*% W)
		dense_matrix::ptr tD = tH->multiply(*tmp1);
		D = tD->transpose();
	}

	// H <- fm.pmax2(H * t(tA %*% W), eps) / (den + eps)
	{
		sparse_matrix::ptr tmat = mat->transpose();
		// tA %*% W
		dense_matrix::ptr tmp1 = multiply(tmat, W);
		// t(tA %*% W)
		tmp1 = tmp1->transpose();
		// H * t(tA %*% W)
		tmp1 = tmp1->multiply_ele(*H);
		// fm.pmax2(H * t(tA %*% W), eps)
		tmp1 = tmp1->pmax_scalar(eps);
		// den + eps
		dense_matrix::ptr tmp3 = D->add_scalar(eps);
		H = tmp1->div(*tmp3);
		H->materialize_self();
	}

	// den <- W %*% (H %*% t(H))
	{
		// t(H)
		dense_matrix::ptr tH = H->transpose();
		// H %*% t(H)
		dense_matrix::ptr tmp1 = H->multiply(*tH);
		D = W->multiply(*tmp1);
	}

	// W <- fm.pmax2(W * (A %*% t(H)), eps) / (den + eps)
	{
		// t(H)
		dense_matrix::ptr tH = H->transpose();
		// A %*% t(H)
		dense_matrix::ptr tmp1 = multiply(mat, tH);
		// W * (A %*% t(H))
		tmp1 = W->multiply_ele(*tmp1);
		// fm.pmax2(W * (A %*% t(H)), eps)
		tmp1 = tmp1->pmax_scalar(eps);
		// den + eps
		dense_matrix::ptr tmp2 = D->add_scalar(eps);
		W = tmp1->div(*tmp2);
		W->materialize_self();
	}
	return std::pair<dense_matrix::ptr, dense_matrix::ptr>(W, H);
}

std::pair<dense_matrix::ptr, dense_matrix::ptr> NMF(sparse_matrix::ptr mat,
		size_t k, size_t max_niters, size_t num_in_mem)
{
	if (mat->get_entry_size() > 0) {
		BOOST_LOG_TRIVIAL(error) << "can only handle binary sparse matrix";
		return std::pair<dense_matrix::ptr, dense_matrix::ptr>();
	}
	size_t n = mat->get_num_rows();
	size_t m = mat->get_num_cols();
	int num_nodes = matrix_conf.get_num_nodes();
	dense_matrix::ptr W = dense_matrix::create_randu<double>(0, 1, n, k,
			matrix_layout_t::L_ROW, num_nodes, true);
	dense_matrix::ptr H = dense_matrix::create_randu<double>(0, 1, k, m,
			matrix_layout_t::L_COL, num_nodes, true);
	std::pair<dense_matrix::ptr, dense_matrix::ptr> res(W, H);
	for (size_t i = 0; i < max_niters; i++) {
		struct timeval start, end;
		gettimeofday(&start, NULL);
		res = update_lee(mat, res.first, res.second);
		gettimeofday(&end, NULL);
		printf("iteration %ld takes %.3f seconds\n", i, time_diff(start, end));
	}
	return res;
}

}

}
