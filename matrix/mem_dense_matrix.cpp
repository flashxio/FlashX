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

#if defined(_OPENMP)
#include <omp.h>
#endif
#include <cblas.h>

#include <atomic>

#include <boost/format.hpp>

#include "log.h"

#include "mem_dense_matrix.h"
#include "generic_type.h"
#include "mem_vector.h"
#include "local_matrix_store.h"
#include "mem_matrix_store.h"
#include "NUMA_dense_matrix.h"

namespace fm
{

mem_dense_matrix::ptr mem_dense_matrix::create(size_t nrow, size_t ncol,
		matrix_layout_t layout, const scalar_type &type, int num_nodes)
{
	detail::matrix_store::ptr store = detail::mem_matrix_store::create(
			nrow, ncol, layout, type, num_nodes);
	return mem_dense_matrix::ptr(new mem_dense_matrix(store));
}

mem_dense_matrix::ptr mem_dense_matrix::create(size_t nrow, size_t ncol,
		matrix_layout_t layout, const scalar_type &type, const set_operate &op,
		int num_nodes)
{
	detail::matrix_store::ptr store = detail::mem_matrix_store::create(
			nrow, ncol, layout, type, num_nodes);
	store->set_data(op);
	return mem_dense_matrix::ptr(new mem_dense_matrix(store));
}

dense_matrix::ptr mem_dense_matrix::get_cols(const std::vector<off_t> &idxs) const
{
	if (store_layout() == matrix_layout_t::L_COL) {
		const detail::mem_matrix_store &store
			= static_cast<const detail::mem_matrix_store &>(get_data());
		return dense_matrix::ptr(new mem_dense_matrix(store.get_cols(idxs)));
	}
	else
		return dense_matrix::ptr();
}

dense_matrix::ptr mem_dense_matrix::get_rows(const std::vector<off_t> &idxs) const
{
	if (store_layout() == matrix_layout_t::L_ROW) {
		const detail::mem_matrix_store &store
			= static_cast<const detail::mem_matrix_store &>(get_data());
		return dense_matrix::ptr(new mem_dense_matrix(store.get_rows(idxs)));
	}
	else
		return dense_matrix::ptr();
}

dense_matrix::ptr mem_dense_matrix::transpose() const
{
	return dense_matrix::ptr(new mem_dense_matrix(get_data().transpose()));
}

dense_matrix::ptr mem_dense_matrix::inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	if (!verify_inner_prod(m, left_op, right_op))
		return dense_matrix::ptr();

	detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
			get_num_rows(), m.get_num_cols(), store_layout(),
			right_op.get_output_type(), get_num_nodes());
	const detail::mem_matrix_store &mem_m
		= dynamic_cast<const detail::mem_matrix_store &>(m.get_data());
	if (is_wide())
		inner_prod_wide(mem_m, left_op, right_op, *res);
	else
		inner_prod_tall(mem_m, left_op, right_op, *res);

	return dense_matrix::ptr(new mem_dense_matrix(res));
}

void mem_dense_matrix::inner_prod_tall(const detail::mem_matrix_store &m,
		const bulk_operate &left_op, const bulk_operate &right_op,
		detail::mem_matrix_store &res) const
{
	// We assume the right matrix is small, so we don't need to partition it.
	detail::local_matrix_store::const_ptr local_right = m.get_portion(0);
	assert(local_right->get_num_rows() == m.get_num_rows()
			&& local_right->get_num_cols() == m.get_num_cols());
	const detail::mem_matrix_store &this_store
		= dynamic_cast<const detail::mem_matrix_store &>(get_data());
	size_t num_chunks = this_store.get_num_portions();
	assert(this_store.get_portion_size().first == res.get_portion_size().first);
#pragma omp parallel for
	for (size_t i = 0; i < num_chunks; i++) {
		detail::local_matrix_store::const_ptr local_store
			= this_store.get_portion(i);
		detail::local_matrix_store::ptr local_res = res.get_portion(i);
		assert(local_store->get_global_start_row()
				== local_res->get_global_start_row());
		assert(local_store->get_global_start_col()
				== local_res->get_global_start_col());
		local_res->reset_data();
		detail::inner_prod(*local_store, *local_right, left_op, right_op,
				*local_res);
	}
}

void mem_dense_matrix::inner_prod_wide(const detail::mem_matrix_store &m,
		const bulk_operate &left_op, const bulk_operate &right_op,
		detail::mem_matrix_store &res) const
{
	size_t nrow = this->get_num_rows();
	int nthreads = get_num_omp_threads();
	std::vector<detail::local_matrix_store::ptr> local_ms(nthreads);

	const detail::mem_matrix_store &this_store
		= dynamic_cast<const detail::mem_matrix_store &>(get_data());
	size_t num_chunks = this_store.get_num_portions();
	assert(this_store.get_portion_size().second == res.get_portion_size().first);
#pragma omp parallel
	{
		detail::local_matrix_store::ptr local_m;
		if (res.store_layout() == matrix_layout_t::L_COL)
			local_m = detail::local_matrix_store::ptr(
					new detail::local_buf_col_matrix_store(0, 0,
						nrow, m.get_num_cols(), right_op.get_output_type()));
		else
			local_m = detail::local_matrix_store::ptr(
					new detail::local_buf_row_matrix_store(0, 0,
						nrow, m.get_num_cols(), right_op.get_output_type()));
		local_m->reset_data();
#pragma omp for
		for (size_t i = 0; i < num_chunks; i++) {
			detail::local_matrix_store::const_ptr local_store
				= this_store.get_portion(i);
			detail::local_matrix_store::const_ptr local_store2 = m.get_portion(i);
			assert(local_store->get_global_start_row()
					== local_store2->get_global_start_col());
			assert(local_store->get_global_start_col()
					== local_store2->get_global_start_row());
			detail::inner_prod(*local_store, *local_store2, left_op, right_op,
					*local_m);
		}
		local_ms[get_omp_thread_num()] = local_m;
	}

	// Aggregate the results from omp threads.
	res.reset_data();
	detail::local_matrix_store::ptr local_res = res.get_portion(0);
	assert(local_res->get_num_rows() == res.get_num_rows()
			&& local_res->get_num_cols() == res.get_num_cols());
	for (int j = 0; j < nthreads; j++)
		detail::mapply2(*local_res, *local_ms[j], right_op, *local_res);
}

scalar_variable::ptr mem_dense_matrix::aggregate(const bulk_operate &op) const
{
	if (!verify_aggregate(op))
		return scalar_variable::ptr();
	scalar_variable::ptr res = op.get_output_type().create_scalar();
	char *raw_arr;

	const detail::mem_matrix_store &this_store
		= dynamic_cast<const detail::mem_matrix_store &>(get_data());
	size_t num_chunks = this_store.get_num_portions();
	if (is_wide()) {
		raw_arr = (char *) malloc(res->get_size() * num_chunks);
#pragma omp parallel for
		for (size_t i = 0; i < num_chunks; i++) {
			detail::local_matrix_store::const_ptr local_store
				= this_store.get_portion(i);
			detail::aggregate(*local_store, op, raw_arr + i * get_entry_size());
		}
	}
	else {
		raw_arr = (char *) malloc(res->get_size() * num_chunks);
#pragma omp parallel for
		for (size_t i = 0; i < num_chunks; i++) {
			detail::local_matrix_store::const_ptr local_store
				= this_store.get_portion(i);
			detail::aggregate(*local_store, op, raw_arr + i * get_entry_size());
		}
	}

	char raw_res[res->get_size()];
	op.runA(num_chunks, raw_arr, raw_res);
	free(raw_arr);
	res->set_raw(raw_res, res->get_size());
	return res;
}

dense_matrix::ptr mem_dense_matrix::mapply2(const dense_matrix &m,
		const bulk_operate &op) const
{
	assert(m.is_in_mem());
	// The same shape and the same data layout.
	if (!verify_mapply2(m, op))
		return dense_matrix::ptr();

	size_t nrow = this->get_num_rows();
	size_t ncol = this->get_num_cols();

	// TODO be careful of NUMA.
	const detail::mem_matrix_store &mem_mat
		= dynamic_cast<const detail::mem_matrix_store &>(m.get_data());
	detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
			nrow, ncol, store_layout(), op.get_output_type(), get_num_nodes());

	const detail::mem_matrix_store &this_store
		= dynamic_cast<const detail::mem_matrix_store &>(get_data());
	size_t num_chunks = this_store.get_num_portions();
	if (is_wide()) {
		assert(this_store.get_portion_size().second
				== mem_mat.get_portion_size().second);
		assert(this_store.get_portion_size().second
				== res->get_portion_size().second);
#pragma omp parallel for
		for (size_t i = 0; i < num_chunks; i++) {
			detail::local_matrix_store::const_ptr local_store
				= this_store.get_portion(i);
			detail::local_matrix_store::const_ptr local_store2
				= mem_mat.get_portion(i);
			detail::local_matrix_store::ptr local_res = res->get_portion(i);
			assert(local_store->get_global_start_col()
					== local_store2->get_global_start_col());
			assert(local_store->get_global_start_col()
					== local_res->get_global_start_col());
			detail::mapply2(*local_store, *local_store2, op, *local_res);
		}
	}
	else {
		assert(this_store.get_portion_size().first
				== mem_mat.get_portion_size().first);
		assert(this_store.get_portion_size().first
				== res->get_portion_size().first);
#pragma omp parallel for
		for (size_t i = 0; i < num_chunks; i++) {
			detail::local_matrix_store::const_ptr local_store
				= this_store.get_portion(i);
			detail::local_matrix_store::const_ptr local_store2
				= mem_mat.get_portion(i);
			detail::local_matrix_store::ptr local_res = res->get_portion(i);
			assert(local_store->get_global_start_row()
					== local_store2->get_global_start_row());
			assert(local_store->get_global_start_row()
					== local_res->get_global_start_row());
			detail::mapply2(*local_store, *local_store2, op, *local_res);
		}
	}
	return mem_dense_matrix::ptr(new mem_dense_matrix(res));
}

dense_matrix::ptr mem_dense_matrix::sapply(const bulk_uoperate &op) const
{
	size_t nrow = this->get_num_rows();
	size_t ncol = this->get_num_cols();
	detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
			nrow, ncol, store_layout(), op.get_output_type(), get_num_nodes());

	const detail::mem_matrix_store &this_store
		= dynamic_cast<const detail::mem_matrix_store &>(get_data());
	size_t num_chunks = this_store.get_num_portions();
	if (is_wide()) {
		assert(this_store.get_portion_size().second
				== res->get_portion_size().second);
#pragma omp parallel for
		for (size_t i = 0; i < num_chunks; i++) {
			detail::local_matrix_store::const_ptr local_store
				= this_store.get_portion(i);
			detail::local_matrix_store::ptr local_res = res->get_portion(i);
			assert(local_res->get_global_start_col()
					== local_store->get_global_start_col());
			detail::sapply(*local_store, op, *local_res);
		}
	}
	else {
		if (this_store.get_portion_size().first != res->get_portion_size().first)
			printf("%ld, %ld\n", this_store.get_portion_size().first,
					res->get_portion_size().first);
		assert(this_store.get_portion_size().first
				== res->get_portion_size().first);
#pragma omp parallel for
		for (size_t i = 0; i < num_chunks; i++) {
			detail::local_matrix_store::const_ptr local_store
				= this_store.get_portion(i);
			detail::local_matrix_store::ptr local_res = res->get_portion(i);
			assert(local_res->get_global_start_row()
					== local_store->get_global_start_row());
			detail::sapply(*local_store, op, *local_res);
		}
	}
	return dense_matrix::ptr(new mem_dense_matrix(res));
}

dense_matrix::ptr mem_dense_matrix::apply(apply_margin margin,
		const arr_apply_operate &op) const
{
	return dense_matrix::ptr();
#if 0
	// Each operation runs on a row
	if (margin == apply_margin::MAR_ROW) {
		size_t out_nrow = nrow;
		size_t out_ncol = op.get_num_out_eles();
		mem_row_dense_matrix::ptr res = mem_row_dense_matrix::create(out_nrow,
				out_ncol, op.get_output_type());
		// We view this row-major matrix as a vector and each row is a subvector.
		mem_vector::ptr res_vec = res->flatten(true);
		// TODO this might be a very large array.
		mem_vector::ptr tmp_vec = mem_vector::create(ncol, get_type());
		for (size_t i = 0; i < nrow; i++) {
			get_row(i, tmp_vec->get_raw_arr());
			bool ret = res_vec->expose_sub_vec(i * res->get_num_cols(),
					res->get_num_cols());
			assert(ret);
			op.run(*tmp_vec, *res_vec);
		}
		return res;
	}
	// Each operation runs on a column
	else {
		size_t out_nrow = op.get_num_out_eles();
		size_t out_ncol = ncol;
		mem_col_dense_matrix::ptr res = mem_col_dense_matrix::create(out_nrow,
				out_ncol, op.get_output_type());
		// We view the input and output matrices as vectors and each column
		// is a subvector.
		mem_vector::ptr in_vec = flatten(false);
		mem_vector::ptr res_vec = res->flatten(false);
		for (size_t i = 0; i < ncol; i++) {
			in_vec->expose_sub_vec(i * get_num_rows(), get_num_rows());
			res_vec->expose_sub_vec(i * res->get_num_rows(), res->get_num_rows());
			op.run(*in_vec, *res_vec);
		}
		return res;
	}
#endif
}

bool mem_dense_matrix::verify_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	if (!m.is_in_mem()) {
		BOOST_LOG_TRIVIAL(error) << "The right matrix isn't in memory";
		return false;
	}
	return dense_matrix::verify_inner_prod(m, left_op, right_op);
}

mem_dense_matrix::ptr mem_dense_matrix::cast(dense_matrix::ptr m)
{
	if (!m->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't cast an EM matrix to mem_dense_matrix";
		return mem_dense_matrix::ptr();
	}
	return std::static_pointer_cast<mem_dense_matrix>(m);
}

#if 0
mem_row_dense_matrix::ptr mem_col_dense_matrix::get_row_store() const
{
	// TODO it works for tall column-wise matrix.
	assert(!is_wide());
	mem_row_dense_matrix::ptr ret = mem_row_dense_matrix::create(get_num_rows(),
			get_num_cols(), get_type());

	size_t entry_size = get_entry_size();
#pragma omp parallel
	{
		std::vector<const char *> col_ptrs(get_num_cols());
#pragma omp for
		for (size_t i = 0; i < get_num_rows(); i++) {
			off_t off = entry_size * i;
			for (size_t j = 0; j < get_num_cols(); j++)
				col_ptrs[j] = this->get_col(j) + off;

			get_type().get_sg().gather(col_ptrs, ret->get_row(i));
		}
	}
	return ret;
}

dense_matrix::ptr mem_col_dense_matrix::gemm(const dense_matrix &Amat,
		const dense_matrix &Bmat, const scalar_variable &alpha,
		const scalar_variable &beta) const
{
	if (Amat.get_num_cols() != Bmat.get_num_rows()
			|| this->get_num_rows() != Amat.get_num_rows()
			|| this->get_num_cols() != Bmat.get_num_cols()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The matrices for gemm don't have compatible dimension";
		return dense_matrix::ptr();
	}
	if (!Amat.is_in_mem() || Amat.store_layout() != matrix_layout_t::L_COL
			|| !Bmat.is_in_mem() || Bmat.store_layout() != matrix_layout_t::L_COL) {
		BOOST_LOG_TRIVIAL(error)
			<< "The A and B matrix needs to be column-major";
		return dense_matrix::ptr();
	}
	if (get_type() != get_scalar_type<double>()
			|| Amat.get_type() != get_scalar_type<double>()
			|| Bmat.get_type() != get_scalar_type<double>()
			|| alpha.get_type() != get_scalar_type<double>()
			|| beta.get_type() != get_scalar_type<double>()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The matrices aren't of double type";
		return dense_matrix::ptr();
	}

	const scalar_variable_impl<double> &d_alpha
		= (const scalar_variable_impl<double> &) alpha;
	const scalar_variable_impl<double> &d_beta
		= (const scalar_variable_impl<double> &) beta;
	mem_col_dense_matrix::ptr ret = mem_col_dense_matrix::create(
			this->get_num_rows(), this->get_num_cols(), get_type());
	ret->copy_from(*this);
	mem_col_dense_matrix::ptr col_Amat
		= dynamic_cast<const mem_col_dense_matrix &>(Amat).get_contig_matrix();
	mem_col_dense_matrix::ptr col_Bmat
		= dynamic_cast<const mem_col_dense_matrix &>(Bmat).get_contig_matrix();
	printf("dgemm\n");
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Amat.get_num_rows(),
			Bmat.get_num_cols(), Amat.get_num_cols(), d_alpha.get(),
			(double *) col_Amat->data.get_raw(), Amat.get_num_rows(),
			(double *) col_Bmat->data.get_raw(), Bmat.get_num_rows(),
			d_beta.get(), (double *) ret->data.get_raw(), ret->get_num_rows());

	return ret;
}

mem_col_dense_matrix::ptr mem_row_dense_matrix::get_col_store() const
{
	// TODO it works for tall row-wise matrix.
	assert(!is_wide());
	mem_col_dense_matrix::ptr ret = mem_col_dense_matrix::create(get_num_rows(),
			get_num_cols(), get_type());

	size_t entry_size = get_entry_size();
#pragma omp parallel
	{
		std::vector<char *> col_ptrs(get_num_cols());
#pragma omp for
		for (size_t i = 0; i < get_num_rows(); i++) {
			off_t off = entry_size * i;
			for (size_t j = 0; j < get_num_cols(); j++)
				col_ptrs[j] = ret->get_col(j) + off;

			get_type().get_sg().scatter(get_row(i), col_ptrs);
		}
	}
	return ret;
}
#endif

}
