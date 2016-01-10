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
#include "combined_matrix_store.h"
#include "cached_matrix_store.h"

namespace fm
{

namespace alg
{

static size_t get_nnz(sparse_matrix::ptr mat)
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
	return scalar_variable::get_val<size_t>(*v);
}

// trace of W %*% H
static double trace_MM(dense_matrix::ptr W, dense_matrix::ptr H)
{
	dense_matrix::ptr tH = H->transpose();
	dense_matrix::ptr tmp = W->multiply_ele(*tH);
	scalar_variable::ptr res = tmp->sum();
	return scalar_variable::get_val<double>(*res);
}

static dense_matrix::ptr multiply(sparse_matrix::ptr S, dense_matrix::ptr D,
		size_t num_in_mem)
{
	size_t k = D->get_num_cols();
	dense_matrix::ptr ret;
	if (num_in_mem >= k * 2) {
		detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
				S->get_num_rows(), D->get_num_cols(), matrix_layout_t::L_ROW,
				D->get_type(), D->get_raw_store()->get_num_nodes());
		S->multiply<double, double>(D->get_raw_store(), res);
		ret = dense_matrix::create(res);
	}
	else if (num_in_mem > k) {
		detail::matrix_store::ptr res = detail::cached_matrix_store::create(
				S->get_num_rows(), D->get_num_cols(),
				D->get_raw_store()->get_num_nodes(), D->get_type(),
				num_in_mem - k, matrix_layout_t::L_COL);
		S->multiply<double, double>(D->get_raw_store(), res);
		ret = dense_matrix::create(res);
	}
	else if (num_in_mem == k) {
		detail::matrix_store::ptr res = detail::matrix_store::create(
				S->get_num_rows(), D->get_num_cols(), matrix_layout_t::L_ROW,
				D->get_type(), -1, false);
		S->multiply<double, double>(D->get_raw_store(), res);
		D->drop_cache();
		ret = dense_matrix::create(res);
	}
	else {
		std::vector<detail::matrix_store::const_ptr> res_mats;
		for (size_t i = 0; i < k; i += num_in_mem) {
			size_t sub_ncol = std::min(num_in_mem, k - i);
			std::vector<off_t> col_idxs(sub_ncol);
			for (size_t j = 0; j < sub_ncol; j++)
				col_idxs[j] = i + j;
			dense_matrix::ptr sub_in = D->get_cols(col_idxs);
			detail::matrix_store::ptr out = detail::matrix_store::create(
					S->get_num_rows(), sub_in->get_num_cols(),
					matrix_layout_t::L_COL, D->get_type(), -1, false);
			S->multiply<double, double>(sub_in->get_raw_store(), out);
			D->drop_cache();
			res_mats.push_back(out);
		}
		ret = dense_matrix::create(detail::combined_matrix_store::create(
					res_mats, matrix_layout_t::L_COL));
	}
	printf("reserve mem for %ld vecs\n",
			detail::get_reserved_bytes() / D->get_num_rows() / D->get_entry_size());
	return ret;
}

// ||A - W %*% H||^2
static double Fnorm(sparse_matrix::ptr A, size_t Annz, dense_matrix::ptr W,
		dense_matrix::ptr H, dense_matrix::ptr tWW, size_t num_in_mem)
{
	// tAW <- t(A) %*% W
	sparse_matrix::ptr tA = A->transpose();
	dense_matrix::ptr tAW = multiply(tA, W, num_in_mem);

	// tHtWW <- t(H) %*% (t(W) %*% W)
	dense_matrix::ptr tH = H->transpose();
	dense_matrix::ptr tHtWW = tH->multiply(*tWW);

	// sumA2 - 2 * trace.MM(tAW, H) + trace.MM(tHtWW, H)
	// TODO H is read twice in this implementation.
	// TODO tAW doesn't need to be materialized.
	return Annz - 2 * trace_MM(tAW, H) + trace_MM(tHtWW, H);
}

struct nmf_state
{
	dense_matrix::ptr W;
	dense_matrix::ptr H;
	dense_matrix::ptr tWW;

	nmf_state(dense_matrix::ptr W, dense_matrix::ptr H,
			dense_matrix::ptr tWW) {
		this->W = W;
		this->H = H;
		this->tWW = tWW;
	}
};

static detail::matrix_store::ptr create_materialize_store(size_t num_rows,
		size_t num_cols, size_t num_in_mem)
{
	int num_nodes = matrix_conf.get_num_nodes();
	if (num_in_mem >= num_cols)
		return detail::matrix_store::create(num_rows, num_cols,
				matrix_layout_t::L_ROW, get_scalar_type<double>(),
				num_nodes, true);
	else
		return detail::cached_matrix_store::create(num_rows, num_cols,
				num_nodes, get_scalar_type<double>(), num_in_mem,
				matrix_layout_t::L_ROW);
}

static nmf_state update_lee(sparse_matrix::ptr mat, dense_matrix::ptr W,
		dense_matrix::ptr H, dense_matrix::ptr tWW, size_t num_in_mem)
{
	double eps = 10e-9;
	// den <- (t(W) %*% W) %*% H
	dense_matrix::ptr D;
	{
		// t(H)
		dense_matrix::ptr tH = H->transpose();
		// t(H) %*% (t(W) %*% W)
		dense_matrix::ptr tD = tH->multiply(*tWW);
		D = tD->transpose();
	}

	// H <- fm.pmax2(H * t(tA %*% W), eps) / (den + eps)
	{
		sparse_matrix::ptr tmat = mat->transpose();
		// tA %*% W
		dense_matrix::ptr tmp1 = multiply(tmat, W, num_in_mem);
		// t(tA %*% W)
		tmp1 = tmp1->transpose();
		// H * t(tA %*% W)
		assert(tmp1->store_layout() == H->store_layout());
		tmp1 = tmp1->multiply_ele(*H);
		// fm.pmax2(H * t(tA %*% W), eps)
		tmp1 = tmp1->pmax_scalar(eps);
		// den + eps
		dense_matrix::ptr tmp3 = D->add_scalar(eps);
		assert(tmp1->store_layout() == tmp3->store_layout());
		H = tmp1->div(*tmp3);
		// mapply_matrix_store materializes a result in a tall matrix buffer.
		// TODO this might be confusing.
		H->set_materialize_level(materialize_level::MATER_FULL,
				create_materialize_store(H->get_num_cols(), H->get_num_rows(),
					num_in_mem));
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
		dense_matrix::ptr tmp1 = multiply(mat, tH, num_in_mem);
		// W * (A %*% t(H))
		assert(W->store_layout() == tmp1->store_layout());
		tmp1 = W->multiply_ele(*tmp1);
		// fm.pmax2(W * (A %*% t(H)), eps)
		tmp1 = tmp1->pmax_scalar(eps);
		// den + eps
		dense_matrix::ptr tmp2 = D->add_scalar(eps);
		assert(tmp1->store_layout() == tmp2->store_layout());
		W = tmp1->div(*tmp2);
		W->set_materialize_level(materialize_level::MATER_FULL,
				create_materialize_store(W->get_num_rows(), W->get_num_cols(),
					num_in_mem));
	}

	dense_matrix::ptr tW = W->transpose();
	tWW = tW->multiply(*W);
	// TODO We have to materialize these matrices here.
	// Otherwise, we'll get segmentation fault. Why?
	W->materialize_self();
	H->materialize_self();
	return nmf_state(W, H, tWW);
}

dense_matrix::ptr create_rand_cached(size_t num_rows, size_t num_cols,
		int num_nodes, size_t num_cached, matrix_layout_t cached_layout)
{
	detail::matrix_store::ptr store = detail::cached_matrix_store::create(
			num_rows, num_cols, num_nodes, get_scalar_type<double>(),
			num_cached, cached_layout);
	store->init_randu(scalar_variable_impl<double>(0),
			scalar_variable_impl<double>(1));
	return dense_matrix::create(store);
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
	size_t nnz = get_nnz(mat);
	int num_nodes = matrix_conf.get_num_nodes();
	dense_matrix::ptr W, H;
	if (num_in_mem >= 2 * k) {
		W = dense_matrix::create_randu<double>(0, 1, n, k,
				matrix_layout_t::L_ROW, num_nodes, true);
		H = dense_matrix::create_randu<double>(0, 1, k, m,
				matrix_layout_t::L_COL, num_nodes, true);
	}
	else if (num_in_mem > k && num_in_mem < 2 * k) {
		W = create_rand_cached(n, k, num_nodes, k, matrix_layout_t::L_ROW);
		H = create_rand_cached(k, m, num_nodes, num_in_mem - k,
				matrix_layout_t::L_ROW);
	}
	else if (num_in_mem == k) {
		W = create_rand_cached(n, k, num_nodes, k, matrix_layout_t::L_ROW);
		H = dense_matrix::create_randu<double>(0, 1, k, m,
				matrix_layout_t::L_COL, num_nodes, false);
	}
	else {
		W = create_rand_cached(n, k, num_nodes, num_in_mem, matrix_layout_t::L_ROW);
		H = dense_matrix::create_randu<double>(0, 1, k, m,
				matrix_layout_t::L_ROW, num_nodes, false);
	}

	dense_matrix::ptr tWW;
	{
		dense_matrix::ptr tW = W->transpose();
		tWW = tW->multiply(*W);
	}
	nmf_state state(W, H, tWW);
	for (size_t i = 0; i < max_niters; i++) {
		struct timeval start, end;
		gettimeofday(&start, NULL);
		state = update_lee(mat, state.W, state.H, state.tWW, num_in_mem);
		gettimeofday(&end, NULL);
		double update_time = time_diff(start, end);
		start = end;
		double dist = Fnorm(mat, nnz, state.W, state.H, state.tWW, num_in_mem);
		gettimeofday(&end, NULL);
		printf("iteration %ld: distance: %f, update time: %.3fs, Fnorm: %.3f\n",
				i, dist, update_time, time_diff(start, end));
	}
	return std::pair<dense_matrix::ptr, dense_matrix::ptr>(state.W, state.H);
}

}

}
