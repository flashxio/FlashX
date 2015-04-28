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
#include "mem_worker_thread.h"

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

namespace
{

class inner_prod_tall_task: public thread_task
{
	detail::local_matrix_store::const_ptr local_right;
	detail::local_matrix_store::const_ptr local_store;
	detail::local_matrix_store::ptr local_res;
	const bulk_operate &left_op;
	const bulk_operate &right_op;
public:
	inner_prod_tall_task(detail::local_matrix_store::const_ptr local_store,
			detail::local_matrix_store::const_ptr local_right,
			const bulk_operate &_left_op, const bulk_operate &_right_op,
			detail::local_matrix_store::ptr local_res): left_op(
				_left_op), right_op(_right_op) {
		this->local_right = local_right;
		this->local_store = local_store;
		this->local_res = local_res;
	}
	void run() {
		local_res->reset_data();
		detail::inner_prod(*local_store, *local_right, left_op, right_op,
				*local_res);
	}
};

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
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	for (size_t i = 0; i < num_chunks; i++) {
		detail::local_matrix_store::const_ptr local_store
			= this_store.get_portion(i);
		detail::local_matrix_store::ptr local_res = res.get_portion(i);
		assert(local_store->get_global_start_row()
				== local_res->get_global_start_row());
		assert(local_store->get_global_start_col()
				== local_res->get_global_start_col());
		assert(local_store->get_node_id() == local_res->get_node_id());
		int node_id = local_store->get_node_id();
		// If the local matrix portion is not assigned to any node, 
		// assign the tasks in round robin fashion.
		if (node_id < 0)
			node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id,
				new inner_prod_tall_task(local_store, local_right,
					left_op, right_op, local_res));
	}
	mem_threads->wait4complete();
}

namespace
{

class inner_prod_wide_task: public thread_task
{
	detail::local_matrix_store::const_ptr local_store;
	detail::local_matrix_store::const_ptr local_store2;
	const bulk_operate &left_op;
	const bulk_operate &right_op;
	const detail::mem_matrix_store &res;
	std::vector<detail::local_matrix_store::ptr> &local_ms;
public:
	inner_prod_wide_task(detail::local_matrix_store::const_ptr local_store,
			detail::local_matrix_store::const_ptr local_store2,
			const bulk_operate &_left_op, const bulk_operate &_right_op,
			const detail::mem_matrix_store &_res,
			std::vector<detail::local_matrix_store::ptr> &_local_ms): left_op(
				_left_op), right_op(_right_op), res(_res), local_ms(_local_ms) {
		this->local_store = local_store;
		this->local_store2 = local_store2;
	}

	void run();
};

void inner_prod_wide_task::run()
{
	detail::pool_task_thread *curr
		= dynamic_cast<detail::pool_task_thread *>(thread::get_curr_thread());
	int thread_id = curr->get_pool_thread_id();
	detail::local_matrix_store::ptr local_m = local_ms[thread_id];
	if (local_m == NULL) {
		int node_id = curr->get_node_id();
		if (res.store_layout() == matrix_layout_t::L_COL)
			local_m = detail::local_matrix_store::ptr(
					new detail::local_buf_col_matrix_store(0, 0,
						res.get_num_rows(), res.get_num_cols(),
						right_op.get_output_type(), node_id));
		else
			local_m = detail::local_matrix_store::ptr(
					new detail::local_buf_row_matrix_store(0, 0,
						res.get_num_rows(), res.get_num_cols(),
						right_op.get_output_type(), node_id));
		local_m->reset_data();
		local_ms[thread_id] = local_m;
	}
	detail::inner_prod(*local_store, *local_store2, left_op, right_op,
			*local_m);
}

}

void mem_dense_matrix::inner_prod_wide(const detail::mem_matrix_store &m,
		const bulk_operate &left_op, const bulk_operate &right_op,
		detail::mem_matrix_store &res) const
{
	assert(this->get_num_rows() == res.get_num_rows());
	assert(m.get_num_cols() == res.get_num_cols());

	const detail::mem_matrix_store &this_store
		= dynamic_cast<const detail::mem_matrix_store &>(get_data());
	size_t num_chunks = this_store.get_num_portions();
	assert(this_store.get_portion_size().second == res.get_portion_size().first);
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	int nthreads = mem_threads->get_num_threads();
	std::vector<detail::local_matrix_store::ptr> local_ms(nthreads);
	for (size_t i = 0; i < num_chunks; i++) {
		detail::local_matrix_store::const_ptr local_store
			= this_store.get_portion(i);
		detail::local_matrix_store::const_ptr local_store2 = m.get_portion(i);
		assert(local_store->get_global_start_row()
				== local_store2->get_global_start_col());
		assert(local_store->get_global_start_col()
				== local_store2->get_global_start_row());
		assert(local_store->get_node_id() == local_store2->get_node_id());
		int node_id = local_store->get_node_id();
		// If the local matrix portion is not assigned to any node, 
		// assign the tasks in round robin fashion.
		if (node_id < 0)
			node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id,
				new inner_prod_wide_task(local_store, local_store2,
					left_op, right_op, res, local_ms));
	}
	mem_threads->wait4complete();

	// Aggregate the results from omp threads.
	res.reset_data();
	detail::local_matrix_store::ptr local_res = res.get_portion(0);
	assert(local_res->get_num_rows() == res.get_num_rows()
			&& local_res->get_num_cols() == res.get_num_cols());
	for (int j = 0; j < nthreads; j++)
		detail::mapply2(*local_res, *local_ms[j], right_op, *local_res);
}

namespace
{

class aggregate_task: public thread_task
{
	detail::local_matrix_store::const_ptr local_store;
	const bulk_operate &op;
	char *local_res;
public:
	aggregate_task(detail::local_matrix_store::const_ptr local_store,
			const bulk_operate &_op, char *local_res): op(_op) {
		this->local_store = local_store;
		this->local_res = local_res;
	}

	void run() {
		detail::aggregate(*local_store, op, local_res);
	}
};

}

scalar_variable::ptr mem_dense_matrix::aggregate(const bulk_operate &op) const
{
	if (!verify_aggregate(op))
		return scalar_variable::ptr();
	scalar_variable::ptr res = op.get_output_type().create_scalar();

	const detail::mem_matrix_store &this_store
		= dynamic_cast<const detail::mem_matrix_store &>(get_data());
	size_t num_chunks = this_store.get_num_portions();
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	std::unique_ptr<char[]> raw_arr(new char[res->get_size() * num_chunks]);
	for (size_t i = 0; i < num_chunks; i++) {
		detail::local_matrix_store::const_ptr local_store
			= this_store.get_portion(i);

		int node_id = local_store->get_node_id();
		// If the local matrix portion is not assigned to any node, 
		// assign the tasks in round robin fashion.
		if (node_id < 0)
			node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id,
				new aggregate_task(local_store, op,
					raw_arr.get() + i * op.output_entry_size()));
	}
	mem_threads->wait4complete();

	char raw_res[res->get_size()];
	op.runA(num_chunks, raw_arr.get(), raw_res);
	res->set_raw(raw_res, res->get_size());
	return res;
}

namespace
{

class mapply2_task: public thread_task
{
	detail::local_matrix_store::const_ptr local_store;
	detail::local_matrix_store::const_ptr local_store2;
	detail::local_matrix_store::ptr local_res;
	const bulk_operate &op;
public:
	mapply2_task(detail::local_matrix_store::const_ptr local_store,
			detail::local_matrix_store::const_ptr local_store2,
			const bulk_operate &_op,
			detail::local_matrix_store::ptr local_res): op(_op) {
		this->local_store = local_store;
		this->local_store2 = local_store2;
		this->local_res = local_res;
	}

	void run() {
		assert(local_store->get_global_start_col()
				== local_store2->get_global_start_col());
		assert(local_store->get_global_start_col()
				== local_res->get_global_start_col());
		assert(local_store->get_global_start_row()
				== local_store2->get_global_start_row());
		assert(local_store->get_global_start_row()
				== local_res->get_global_start_row());
		detail::mapply2(*local_store, *local_store2, op, *local_res);
	}
};

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
	}
	else {
		assert(this_store.get_portion_size().first
				== mem_mat.get_portion_size().first);
		assert(this_store.get_portion_size().first
				== res->get_portion_size().first);
	}

	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	for (size_t i = 0; i < num_chunks; i++) {
		detail::local_matrix_store::const_ptr local_store
			= this_store.get_portion(i);
		detail::local_matrix_store::const_ptr local_store2
			= mem_mat.get_portion(i);
		detail::local_matrix_store::ptr local_res = res->get_portion(i);

		int node_id = local_store->get_node_id();
		// If the local matrix portion is not assigned to any node, 
		// assign the tasks in round robin fashion.
		if (node_id < 0)
			node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id,
				new mapply2_task(local_store, local_store2, op, local_res));
	}
	mem_threads->wait4complete();
	return mem_dense_matrix::ptr(new mem_dense_matrix(res));
}

namespace
{

class sapply_task: public thread_task
{
	detail::local_matrix_store::const_ptr local_store;
	detail::local_matrix_store::ptr local_res;
	const bulk_uoperate &op;
public:
	sapply_task(detail::local_matrix_store::const_ptr local_store,
			const bulk_uoperate &_op,
			detail::local_matrix_store::ptr local_res): op(_op) {
		this->local_store = local_store;
		this->local_res = local_res;
	}

	void run() {
		assert(local_store->get_global_start_col()
				== local_res->get_global_start_col());
		assert(local_store->get_global_start_row()
				== local_res->get_global_start_row());
		detail::sapply(*local_store, op, *local_res);
	}
};

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
	if (is_wide())
		assert(this_store.get_portion_size().second
				== res->get_portion_size().second);
	else
		assert(this_store.get_portion_size().first
				== res->get_portion_size().first);

	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	for (size_t i = 0; i < num_chunks; i++) {
		detail::local_matrix_store::const_ptr local_store
			= this_store.get_portion(i);
		detail::local_matrix_store::ptr local_res = res->get_portion(i);

		int node_id = local_store->get_node_id();
		// If the local matrix portion is not assigned to any node, 
		// assign the tasks in round robin fashion.
		if (node_id < 0)
			node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id,
				new sapply_task(local_store, op, local_res));
	}
	mem_threads->wait4complete();
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

namespace
{

class dgemm_task: public thread_task
{
	detail::local_matrix_store::const_ptr local_store;
	detail::local_matrix_store::const_ptr local_Astore;
	detail::local_matrix_store::const_ptr local_Bstore;
	detail::local_matrix_store::ptr local_res;
	double alpha;
	double beta;
public:
	dgemm_task(detail::local_matrix_store::const_ptr local_store,
			detail::local_matrix_store::const_ptr local_Astore,
			detail::local_matrix_store::const_ptr local_Bstore,
			detail::local_matrix_store::ptr local_res,
			double alpha, double beta) {
		this->local_store = local_store;
		this->local_Astore = local_Astore;
		this->local_Bstore = local_Bstore;
		this->local_res = local_res;
		this->alpha = alpha;
		this->beta = beta;
	}

	void run();
};

void dgemm_task::run()
{
	const double *Amat = (const double *) local_Astore->get_raw_arr();
	const double *Bmat = (const double *) local_Bstore->get_raw_arr();
	double *res_mat = (double *) local_res->get_raw_arr();
	if (Amat == NULL) {
		detail::local_matrix_store::ptr tmp(
				new detail::local_buf_col_matrix_store(
					local_Astore->get_global_start_row(),
					local_Astore->get_global_start_col(),
					local_Astore->get_num_rows(), local_Astore->get_num_cols(),
					local_Astore->get_type(), local_Astore->get_node_id()));
		tmp->copy_from(*local_Astore);
		local_Astore = tmp;
		Amat = (const double *) local_Astore->get_raw_arr();
		assert(Amat);
	}
	if (Bmat == NULL) {
		detail::local_matrix_store::ptr tmp(
				new detail::local_buf_col_matrix_store(
					local_Bstore->get_global_start_row(),
					local_Bstore->get_global_start_col(),
					local_Bstore->get_num_rows(), local_Bstore->get_num_cols(),
					local_Bstore->get_type(), local_Bstore->get_node_id()));
		tmp->copy_from(*local_Bstore);
		local_Bstore = tmp;
		Bmat = (const double *) local_Bstore->get_raw_arr();
		assert(Bmat);
	}
	detail::local_matrix_store::ptr orig_res_store;
	if (res_mat == NULL) {
		orig_res_store = local_res;
		detail::local_matrix_store::ptr tmp(
				new detail::local_buf_col_matrix_store(
					local_res->get_global_start_row(),
					local_res->get_global_start_col(),
					local_res->get_num_rows(), local_res->get_num_cols(),
					local_res->get_type(), local_res->get_node_id()));
		if (beta != 0)
			tmp->copy_from(*local_store);
		local_res = tmp;
		res_mat = (double *) local_res->get_raw_arr();
		assert(res_mat);
	}
	else {
		if (beta != 0)
			local_res->copy_from(*local_store);
	}
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			local_Astore->get_num_rows(), local_Bstore->get_num_cols(),
			local_Astore->get_num_cols(), alpha, Amat,
			local_Astore->get_num_rows(), Bmat, local_Bstore->get_num_rows(),
			beta, res_mat, local_res->get_num_rows());
	if (orig_res_store)
		orig_res_store->copy_from(*local_res);
}

}

dense_matrix::ptr mem_dense_matrix::gemm(const dense_matrix &Amat,
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
			|| !Bmat.is_in_mem() || Bmat.store_layout() != matrix_layout_t::L_COL
			|| store_layout() != matrix_layout_t::L_COL) {
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

	const detail::mem_matrix_store &Amat_store
		= static_cast<const detail::mem_matrix_store &>(Amat.get_data());
	const detail::mem_matrix_store &Bmat_store
		= static_cast<const detail::mem_matrix_store &>(Bmat.get_data());
	const detail::mem_matrix_store &this_store
		= dynamic_cast<const detail::mem_matrix_store &>(get_data());
	const scalar_variable_impl<double> &d_alpha
		= static_cast<const scalar_variable_impl<double> &>(alpha);
	const scalar_variable_impl<double> &d_beta
		= static_cast<const scalar_variable_impl<double> &>(beta);
	detail::mem_matrix_store::ptr res_store = detail::mem_matrix_store::create(
			this->get_num_rows(), this->get_num_cols(),
			matrix_layout_t::L_COL, get_type(), get_num_nodes());

	// We assume the right matrix is small, so we don't need to partition it.
	detail::local_matrix_store::const_ptr local_Bstore = Bmat_store.get_portion(0);
	assert(local_Bstore->get_num_rows() == Bmat_store.get_num_rows()
			&& local_Bstore->get_num_cols() == Bmat_store.get_num_cols());
	size_t num_chunks = this_store.get_num_portions();
	assert(this_store.get_portion_size().first == Amat_store.get_portion_size().first);
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	for (size_t i = 0; i < num_chunks; i++) {
		detail::local_matrix_store::const_ptr local_store
			= this_store.get_portion(i);
		detail::local_matrix_store::const_ptr local_Astore
			= Amat_store.get_portion(i);
		detail::local_matrix_store::ptr local_res = res_store->get_portion(i);
		assert(local_store->get_global_start_row()
				== local_res->get_global_start_row());
		assert(local_store->get_global_start_col()
				== local_res->get_global_start_col());
		assert(local_store->get_node_id() == local_res->get_node_id());

		int node_id = local_store->get_node_id();
		// If the local matrix portion is not assigned to any node, 
		// assign the tasks in round robin fashion.
		if (node_id < 0)
			node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id,
				new dgemm_task(local_store, local_Astore, local_Bstore,
					local_res, d_alpha.get(), d_beta.get()));
	}
	mem_threads->wait4complete();

	return dense_matrix::ptr(new mem_dense_matrix(res_store));
}

namespace
{

class scale_col_task: public thread_task
{
	detail::local_matrix_store::const_ptr local_store;
	const mem_vector &vals;
	detail::local_matrix_store::ptr local_res;
public:
	scale_col_task(detail::local_matrix_store::const_ptr local_store,
			const mem_vector &_vals,
			detail::local_matrix_store::ptr local_res): vals(_vals) {
		this->local_store = local_store;
		this->local_res = local_res;
	}

	void run() {
		assert(local_store->get_global_start_col()
				== local_res->get_global_start_col());
		assert(local_store->get_global_start_row()
				== local_res->get_global_start_row());
		detail::scale_cols(*local_store, vals, *local_res);
	}
};

}

dense_matrix::ptr mem_dense_matrix::scale_cols(const mem_vector &vals) const
{
	assert(!is_wide());
	assert(get_num_cols() == vals.get_length());
	size_t nrow = this->get_num_rows();
	size_t ncol = this->get_num_cols();
	detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
			nrow, ncol, store_layout(), get_type(), get_num_nodes());

	const detail::mem_matrix_store &this_store
		= dynamic_cast<const detail::mem_matrix_store &>(get_data());
	size_t num_chunks = this_store.get_num_portions();

	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	for (size_t i = 0; i < num_chunks; i++) {
		detail::local_matrix_store::const_ptr local_store
			= this_store.get_portion(i);
		detail::local_matrix_store::ptr local_res = res->get_portion(i);

		int node_id = local_store->get_node_id();
		// If the local matrix portion is not assigned to any node, 
		// assign the tasks in round robin fashion.
		if (node_id < 0)
			node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id,
				new scale_col_task(local_store, vals, local_res));
	}
	mem_threads->wait4complete();
	return dense_matrix::ptr(new mem_dense_matrix(res));
}

}
