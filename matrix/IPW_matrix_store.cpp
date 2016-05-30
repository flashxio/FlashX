/*
 * Copyright 2016 Open Connectome Project (http://openconnecto.me)
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

#include <cblas.h>

#include "IPW_matrix_store.h"
#include "dense_matrix.h"
#include "local_matrix_store.h"
#include "project_matrix_store.h"

namespace fm
{

namespace detail
{

namespace
{

struct matrix_info
{
	size_t num_rows;
	size_t num_cols;
	matrix_layout_t layout;
};

class combine_op: public detail::portion_mapply_op
{
public:
	combine_op(size_t num_rows, size_t num_cols,
			const scalar_type &type): detail::portion_mapply_op(
				num_rows, num_cols, type) {
	}
	virtual bool has_materialized() const = 0;
	virtual detail::mem_matrix_store::ptr get_combined_result() const = 0;
	virtual void set_require_trans(bool val) = 0;
};

class inner_prod_wide_op: public combine_op
{
	bool require_trans;
	bulk_operate::const_ptr left_op;
	bulk_operate::const_ptr right_op;
	matrix_info out_mat_info;
	std::vector<detail::local_matrix_store::ptr> local_ms;
	std::vector<detail::local_matrix_store::ptr> local_tmps;
public:
	inner_prod_wide_op(bulk_operate::const_ptr left_op,
			bulk_operate::const_ptr right_op, const matrix_info &out_mat_info,
			size_t num_threads): combine_op(0, 0, right_op->get_output_type()) {
		this->left_op = left_op;
		this->right_op = right_op;
		local_ms.resize(num_threads);
		local_tmps.resize(num_threads);
		this->out_mat_info = out_mat_info;
		this->require_trans = false;
	}

	void set_require_trans(bool val) {
		this->require_trans = val;
	}

	virtual bool has_materialized() const {
		bool materialized = false;
		for (size_t i = 0; i < local_ms.size(); i++)
			if (local_ms[i])
				materialized = true;
		return materialized;
	}

	virtual detail::mem_matrix_store::ptr get_combined_result() const;

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const;

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		// We don't need to implement this because we materialize
		// the output matrix immediately.
		assert(0);
		return detail::portion_mapply_op::const_ptr();
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("inner_prod(") + mats[0]->get_name()
			+ "," + mats[1]->get_name() + ")";
	}
};

detail::mem_matrix_store::ptr inner_prod_wide_op::get_combined_result() const
{
	// The first non-empty local matrix.
	detail::local_matrix_store::ptr lmat;
	for (size_t i = 0; i < local_ms.size(); i++)
		if (local_ms[i]) {
			lmat = local_ms[i];
			break;
		}
	assert(lmat);

	// Aggregate the results from omp threads.
	detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
			lmat->get_num_rows(), lmat->get_num_cols(), lmat->store_layout(),
			right_op->get_output_type(), -1);
	detail::local_matrix_store::ptr local_res = res->get_portion(0);
	assert(local_res->get_num_rows() == res->get_num_rows()
			&& local_res->get_num_cols() == res->get_num_cols());
	res->write_portion_async(lmat, 0, 0);

	for (size_t j = 0; j < local_ms.size(); j++) {
		// It's possible that the local matrix store doesn't exist
		// because the input matrix is very small.
		if (local_ms[j] && local_ms[j] != lmat)
			detail::mapply2(*local_res, *local_ms[j], *right_op, *local_res);
	}
	return res;
}

void inner_prod_wide_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins) const
{
	int thread_id = detail::mem_thread_pool::get_curr_thread_id();
	detail::local_matrix_store::ptr local_m = local_ms[thread_id];
	detail::local_matrix_store::ptr local_tmp = local_tmps[thread_id];
	bool is_first = false;
	if (local_m == NULL) {
		is_first = true;
		if (out_mat_info.layout == matrix_layout_t::L_COL) {
			local_m = detail::local_matrix_store::ptr(
					new detail::local_buf_col_matrix_store(0, 0,
						out_mat_info.num_rows, out_mat_info.num_cols,
						right_op->get_output_type(), -1));
			local_tmp = detail::local_matrix_store::ptr(
					new detail::local_buf_col_matrix_store(0, 0,
						out_mat_info.num_rows, out_mat_info.num_cols,
						right_op->get_output_type(), -1));
		}
		else {
			local_m = detail::local_matrix_store::ptr(
					new detail::local_buf_row_matrix_store(0, 0,
						out_mat_info.num_rows, out_mat_info.num_cols,
						right_op->get_output_type(), -1));
			local_tmp = detail::local_matrix_store::ptr(
					new detail::local_buf_row_matrix_store(0, 0,
						out_mat_info.num_rows, out_mat_info.num_cols,
						right_op->get_output_type(), -1));
		}
		local_m->reset_data();
		local_tmp->reset_data();
		assert((size_t) thread_id < local_ms.size());
		const_cast<inner_prod_wide_op *>(this)->local_ms[thread_id] = local_m;
		const_cast<inner_prod_wide_op *>(this)->local_tmps[thread_id] = local_tmp;
	}
	if (!require_trans) {
		assert(ins[0]->get_num_cols() == ins[1]->get_num_rows());
		if (is_first)
			detail::inner_prod(*ins[0], *ins[1], *left_op, *right_op, *local_m);
		else {
			detail::inner_prod(*ins[0], *ins[1], *left_op, *right_op, *local_tmp);
			mapply2(*local_m, *local_tmp, *right_op, *local_m);
		}
	}
	else {
		assert(ins[0]->get_num_rows() == ins[1]->get_num_rows());
		detail::local_matrix_store::const_ptr store
			= std::static_pointer_cast<const detail::local_matrix_store>(
					ins[0]->transpose());
		if (is_first)
			detail::inner_prod(*store, *ins[1], *left_op, *right_op, *local_m);
		else {
			detail::inner_prod(*store, *ins[1], *left_op, *right_op, *local_tmp);
			mapply2(*local_m, *local_tmp, *right_op, *local_m);
		}
	}
}

/*
 * This class accumulates the matrix multiplication results on matrix partitions.
 * For float-point matrices, we require certain precision. When we accumulate
 * the multiply results, we lose precision if we use the original float-point type.
 * As such, we should use a high-precision float-point for accumulation.
 */
class matmul_accumulator
{
public:
	typedef std::shared_ptr<matmul_accumulator> ptr;

	static ptr create(size_t num_rows, size_t num_cols, matrix_layout_t layout,
			const scalar_type &type);

	/*
	 * We accumulate the result from each partition.
	 */
	virtual void add_matrix(const detail::local_matrix_store &mat) = 0;

	/*
	 * We need to combine the results from multiple accumulators to generate
	 * the final result.
	 */
	virtual detail::mem_matrix_store::ptr combine(
			const std::vector<matmul_accumulator::ptr> &accus) const = 0;

	virtual detail::local_matrix_store::ptr get_accu() const = 0;
};

template<class ExposeType, class IntType>
class matmul_accumulator_impl: public matmul_accumulator
{
	detail::local_matrix_store::ptr accu_buf;
public:
	matmul_accumulator_impl(size_t num_rows, size_t num_cols,
			matrix_layout_t layout);
	virtual void add_matrix(const detail::local_matrix_store &mat);
	virtual detail::mem_matrix_store::ptr combine(
			const std::vector<matmul_accumulator::ptr> &accus) const;
	virtual detail::local_matrix_store::ptr get_accu() const {
		return accu_buf;
	}
};

template<class ExposeType, class IntType>
matmul_accumulator_impl<ExposeType, IntType>::matmul_accumulator_impl(
		size_t num_rows, size_t num_cols, matrix_layout_t layout)
{
	if (layout == matrix_layout_t::L_ROW)
		accu_buf = detail::local_matrix_store::ptr(
				new detail::local_buf_row_matrix_store(0, 0,
					num_rows, num_cols, get_scalar_type<IntType>(), -1));
	else
		accu_buf = detail::local_matrix_store::ptr(
				new detail::local_buf_col_matrix_store(0, 0,
					num_rows, num_cols, get_scalar_type<IntType>(), -1));
	accu_buf->reset_data();
}

template<class ExposeType, class IntType>
void matmul_accumulator_impl<ExposeType, IntType>::add_matrix(
		const detail::local_matrix_store &mat)
{
	assert(accu_buf->store_layout() == mat.store_layout());
	const ExposeType *input_arr = reinterpret_cast<const ExposeType *>(
			mat.get_raw_arr());
	IntType *accu_arr = reinterpret_cast<IntType *>(accu_buf->get_raw_arr());
	assert(input_arr && accu_arr);
	size_t num_eles = mat.get_num_rows() * mat.get_num_cols();
	for (size_t i = 0; i < num_eles; i++)
		accu_arr[i] += input_arr[i];
}

template<class ExposeType, class IntType>
detail::mem_matrix_store::ptr matmul_accumulator_impl<ExposeType, IntType>::combine(
			const std::vector<matmul_accumulator::ptr> &accus) const
{
	detail::local_matrix_store::ptr accu_buf = accus[0]->get_accu();
	matrix_layout_t layout = accu_buf->store_layout();
	size_t num_rows = accu_buf->get_num_rows();
	size_t num_cols = accu_buf->get_num_cols();
	detail::mem_matrix_store::ptr accu_mat = detail::mem_matrix_store::create(
				num_rows, num_cols, layout, get_scalar_type<IntType>(), -1);
	detail::mem_matrix_store::ptr expo_mat = detail::mem_matrix_store::create(
				num_rows, num_cols, layout, get_scalar_type<ExposeType>(), -1);
	expo_mat->reset_data();
	accu_mat->reset_data();

	// If there is only one accumulator, we convert the accumulated result
	// to the exposed type and return.
	size_t num_eles = num_rows * num_cols;
	if (accus.size() == 1) {
		ExposeType *expo_arr = reinterpret_cast<ExposeType *>(
				expo_mat->get_raw_arr());
		const IntType *accu_arr = reinterpret_cast<const IntType *>(
				accu_buf->get_raw_arr());
		for (size_t i = 0; i < num_eles; i++)
			expo_arr[i] = accu_arr[i];
		return expo_mat;
	}

	// Otherwise, we need to add the results from multiple accumulators.
	IntType *final_accu = reinterpret_cast<IntType *>(accu_mat->get_raw_arr());
	for (size_t i = 0; i < accus.size(); i++) {
		const IntType *accu_arr = reinterpret_cast<IntType *>(
				accus[i]->get_accu()->get_raw_arr());
		for (size_t j = 0; j < num_eles; j++)
			final_accu[j] += accu_arr[j];
	}

	// And convert the final results to the exposed type.
	ExposeType *expo_arr = reinterpret_cast<ExposeType *>(expo_mat->get_raw_arr());
	for (size_t i = 0; i < num_eles; i++)
		expo_arr[i] = final_accu[i];
	return expo_mat;
}

matmul_accumulator::ptr matmul_accumulator::create(size_t num_rows,
		size_t num_cols, matrix_layout_t layout, const scalar_type &type)
{
	if (type == get_scalar_type<double>())
		return ptr(new matmul_accumulator_impl<double, long double>(
					num_rows, num_cols, layout));
	else
		return ptr(new matmul_accumulator_impl<float, double>(
					num_rows, num_cols, layout));
}

class multiply_wide_op: public combine_op
{
	std::vector<detail::local_matrix_store::ptr> Abufs;
	std::vector<detail::local_matrix_store::ptr> Bbufs;
	std::vector<detail::local_matrix_store::ptr> tmp_bufs;
	std::vector<matmul_accumulator::ptr> res_bufs;
	bool require_trans;
	bool is_sparse;
	size_t out_num_rows;
	size_t out_num_cols;
	matrix_layout_t Alayout;
	matrix_layout_t Blayout;
public:
	multiply_wide_op(size_t num_threads, size_t out_num_rows, size_t out_num_cols,
			matrix_layout_t required_layout, const scalar_type &type,
			bool is_sparse): combine_op(0, 0, type) {
		Abufs.resize(num_threads);
		Bbufs.resize(num_threads);
		res_bufs.resize(num_threads);
		tmp_bufs.resize(num_threads);
		this->out_num_rows = out_num_rows;
		this->out_num_cols = out_num_cols;
		Alayout = required_layout;
		Blayout = required_layout;
		require_trans = false;
		this->is_sparse = is_sparse;
	}

	void set_require_trans(bool val) {
		if (require_trans == val)
			return;

		this->require_trans = val;
		// We need to transpose the A matrix, so we want the data in the A matrix
		// to be organized in the opposite layout to the required.
		if (Alayout == matrix_layout_t::L_COL)
			Alayout = matrix_layout_t::L_ROW;
		else
			Alayout = matrix_layout_t::L_COL;
	}

	virtual bool has_materialized() const {
		bool materialized = false;
		for (size_t i = 0; i < res_bufs.size(); i++)
			if (res_bufs[i])
				materialized = true;
		return materialized;
	}

	virtual detail::mem_matrix_store::ptr get_combined_result() const;

	void run_part_dense(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const;
	void run_sparse(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const;

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const;

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		assert(0);
		return detail::portion_mapply_op::const_ptr();
	}

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 2);
		return std::string("(") + (mats[0]->get_name()
					+ "*") + mats[1]->get_name() + std::string(")");
	}
};

detail::mem_matrix_store::ptr multiply_wide_op::get_combined_result() const
{
	std::vector<matmul_accumulator::ptr> non_empty;
	for (size_t i = 0; i < res_bufs.size(); i++)
		if (res_bufs[i])
			non_empty.push_back(res_bufs[i]);
	assert(non_empty.size() > 0);
	return non_empty[0]->combine(non_empty);
}

template<class T>
void wide_gemm_col(const std::pair<size_t, size_t> &Asize, const T *Amat,
		const std::pair<size_t, size_t> &, const T *Bmat,
		T *res_mat, size_t out_num_rows)
{
	assert(0);
}

template<>
void wide_gemm_col<double>(const std::pair<size_t, size_t> &Asize,
		const double *Amat, const std::pair<size_t, size_t> &Bsize,
		const double *Bmat, double *res_mat, size_t out_num_rows)
{
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			Asize.first, Bsize.second, Asize.second, 1, Amat, Asize.first,
			Bmat, Bsize.first, 0, res_mat, out_num_rows);
}

template<>
void wide_gemm_col<float>(const std::pair<size_t, size_t> &Asize,
		const float *Amat, const std::pair<size_t, size_t> &Bsize,
		const float *Bmat, float *res_mat, size_t out_num_rows)
{
	cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			Asize.first, Bsize.second, Asize.second, 1, Amat, Asize.first,
			Bmat, Bsize.first, 0, res_mat, out_num_rows);
}

template<class T>
void wide_gemm_row(const std::pair<size_t, size_t> &Asize, const T *Amat,
		const std::pair<size_t, size_t> &Bsize, const T *Bmat,
		T *res_mat, size_t out_num_cols)
{
	assert(0);
}

template<>
void wide_gemm_row<double>(const std::pair<size_t, size_t> &Asize,
		const double *Amat, const std::pair<size_t, size_t> &Bsize,
		const double *Bmat, double *res_mat, size_t out_num_cols)
{
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			Asize.first, Bsize.second, Asize.second, 1, Amat, Asize.second,
			Bmat, Bsize.second, 0, res_mat, out_num_cols);
}

template<>
void wide_gemm_row<float>(const std::pair<size_t, size_t> &Asize,
		const float *Amat, const std::pair<size_t, size_t> &Bsize,
		const float *Bmat, float *res_mat, size_t out_num_cols)
{
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			Asize.first, Bsize.second, Asize.second, 1, Amat, Asize.second,
			Bmat, Bsize.second, 0, res_mat, out_num_cols);
}

/*
 * dst_col += B * A[col_idx]
 */
template<class T>
void cblas_axpy(size_t len, T B, const T *Acol, T *dst_col)
{
	assert(0);
}
template<>
void cblas_axpy(size_t len, double B, const double *Acol, double *dst_col)
{
	cblas_daxpy(len, B, Acol, 1, dst_col, 1);
}
template<>
void cblas_axpy(size_t len, float B, const float *Acol, float *dst_col)
{
	cblas_saxpy(len, B, Acol, 1, dst_col, 1);
}

/*
 * This multiplies a dense matrix A with a sparse matrix B.
 */
template<class T>
void multiply_sparse(const detail::local_col_matrix_store &Astore,
		const detail::local_sparse_matrix_store &Bstore,
		detail::local_col_matrix_store &Cstore)
{
	size_t LONG_DIM_LEN = get_long_dim_len(Astore);
	assert(get_scalar_type<T>() == Bstore.get_type());
	std::vector<off_t> col_idxs;
	size_t orig_num_cols = Astore.get_num_cols();
	local_matrix_store::exposed_area orig_left = Astore.get_exposed_area();
	local_col_matrix_store &mutableA = const_cast<local_col_matrix_store &>(
			Astore);
	for (size_t col_idx = 0; col_idx < orig_num_cols; col_idx += LONG_DIM_LEN) {
		size_t llen = std::min(orig_num_cols - col_idx, LONG_DIM_LEN);
		mutableA.resize(orig_left.local_start_row,
				orig_left.local_start_col + col_idx, Astore.get_num_rows(), llen);
		for (size_t i = 0; i < Astore.get_num_cols(); i++) {
			col_idxs.clear();
			off_t Brow_idx = col_idx + i;
			const char *_Brow = Bstore.get_row_nnz(Brow_idx, col_idxs);
			// If the row doesn't have nnz.
			if (_Brow == NULL)
				continue;

			const T *Brow = reinterpret_cast<const T *>(_Brow);
			for (size_t j = 0; j < col_idxs.size(); j++) {
				T *dst_col = reinterpret_cast<T *>(Cstore.get_col(col_idxs[j]));
				const T *Acol = reinterpret_cast<const T *>(Astore.get_col(i));
				cblas_axpy<T>(Astore.get_num_rows(), Brow[j], Acol, dst_col);
			}
		}
	}
	mutableA.restore_size(orig_left);
}

/*
 * This multiplies a dense matrix t(A) with a sparse matrix B.
 */
template<class T>
void multiply_sparse_trans(const detail::local_row_matrix_store &Astore,
		const detail::local_sparse_matrix_store &Bstore,
		detail::local_col_matrix_store &Cstore)
{
	size_t LONG_DIM_LEN = get_long_dim_len(Astore);
	assert(get_scalar_type<T>() == Bstore.get_type());
	std::vector<off_t> col_idxs;
	size_t orig_num_rows = Astore.get_num_rows();
	local_matrix_store::exposed_area orig_left = Astore.get_exposed_area();
	local_row_matrix_store &mutableA = const_cast<local_row_matrix_store &>(
			Astore);
	for (size_t row_idx = 0; row_idx < orig_num_rows; row_idx += LONG_DIM_LEN) {
		size_t llen = std::min(orig_num_rows - row_idx, LONG_DIM_LEN);
		mutableA.resize(orig_left.local_start_row + row_idx,
				orig_left.local_start_col, llen, Astore.get_num_cols());
		for (size_t i = 0; i < Astore.get_num_rows(); i++) {
			col_idxs.clear();
			off_t Brow_idx = row_idx + i;
			const char *_Brow = Bstore.get_row_nnz(Brow_idx, col_idxs);
			// If the row doesn't have nnz.
			if (_Brow == NULL)
				continue;

			const T *Brow = reinterpret_cast<const T *>(_Brow);
			for (size_t j = 0; j < col_idxs.size(); j++) {
				T *dst_col = reinterpret_cast<T *>(Cstore.get_col(col_idxs[j]));
				const T *Arow = reinterpret_cast<const T *>(Astore.get_row(i));
				cblas_axpy<T>(Astore.get_num_cols(), Brow[j], Arow, dst_col);
			}
		}
	}
	mutableA.restore_size(orig_left);
}

void multiply_wide_op::run_sparse(
		const std::vector<detail::local_matrix_store::const_ptr> &ins) const
{
	int thread_id = detail::mem_thread_pool::get_curr_thread_id();
	multiply_wide_op *mutable_this = const_cast<multiply_wide_op *>(this);
	detail::local_matrix_store::const_ptr left = ins[0];
	const detail::local_sparse_matrix_store &Bstore
		= static_cast<const detail::local_sparse_matrix_store &>(*ins[1]);

	if (res_bufs[thread_id] == NULL) {
		mutable_this->tmp_bufs[thread_id] = detail::local_matrix_store::ptr(
				new local_buf_col_matrix_store(0, 0,
					out_num_rows, out_num_cols, get_output_type(), -1));
		mutable_this->res_bufs[thread_id] = matmul_accumulator::create(
				out_num_rows, out_num_cols, Blayout, get_output_type());
	}
	assert(tmp_bufs[thread_id]->store_layout() == Blayout);
	tmp_bufs[thread_id]->reset_data();

	local_col_matrix_store &tmp_col_buf = static_cast<local_col_matrix_store &>(
			*tmp_bufs[thread_id]);
	if (get_output_type() == get_scalar_type<double>() && require_trans) {
		assert(left->store_layout() == matrix_layout_t::L_ROW);
		const detail::local_row_matrix_store &Astore
			= static_cast<const detail::local_row_matrix_store &>(*left);
		multiply_sparse_trans<double>(Astore, Bstore, tmp_col_buf);
	}
	else if (get_output_type() == get_scalar_type<double>()) {
		assert(left->store_layout() == matrix_layout_t::L_COL);
		const detail::local_col_matrix_store &Astore
			= static_cast<const detail::local_col_matrix_store &>(*left);
		multiply_sparse<double>(Astore, Bstore, tmp_col_buf);
	}
	else if (require_trans) {
		assert(get_output_type() == get_scalar_type<float>());
		assert(left->store_layout() == matrix_layout_t::L_ROW);
		const detail::local_row_matrix_store &Astore
			= static_cast<const detail::local_row_matrix_store &>(*left);
		multiply_sparse_trans<float>(Astore, Bstore, tmp_col_buf);
	}
	else {
		assert(get_output_type() == get_scalar_type<float>());
		assert(left->store_layout() == matrix_layout_t::L_COL);
		const detail::local_col_matrix_store &Astore
			= static_cast<const detail::local_col_matrix_store &>(*left);
		multiply_sparse<float>(Astore, Bstore, tmp_col_buf);
	}
	res_bufs[thread_id]->add_matrix(*tmp_bufs[thread_id]);
}

void multiply_wide_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins) const
{
	assert(ins.size() == 2);
	// If one of the matrix is a sparse matrix, some of its portions might
	// be empty.
	if (ins[0] == NULL || ins[1] == NULL)
		return;

	size_t LONG_DIM_LEN = get_long_dim_len(*ins[0], *ins[1]);

	size_t long_dim = ins[1]->get_num_rows();
	if (long_dim <= LONG_DIM_LEN) {
		run_part_dense(ins);
		return;
	}

	// We need to reset the temporary matrix that stores the product every time,
	// which can be expense relative to the computation of the matrix
	// multiplication.
	if (is_sparse) {
		run_sparse(ins);
		return;
	}

	local_matrix_store::exposed_area orig_A = ins[0]->get_exposed_area();
	local_matrix_store::exposed_area orig_B = ins[1]->get_exposed_area();
	local_matrix_store &mutableA = const_cast<local_matrix_store &>(*ins[0]);
	local_matrix_store &mutableB = const_cast<local_matrix_store &>(*ins[1]);
	for (size_t row_idx = 0; row_idx < long_dim; row_idx += LONG_DIM_LEN) {
		size_t llen = std::min(long_dim - row_idx, LONG_DIM_LEN);
		if (require_trans)
			mutableA.resize(orig_A.local_start_row + row_idx,
					orig_A.local_start_col, llen, mutableA.get_num_cols());
		else
			mutableA.resize(orig_A.local_start_row,
					orig_A.local_start_col + row_idx, mutableA.get_num_rows(),
					llen);
		mutableB.resize(orig_B.local_start_row + row_idx,
				orig_B.local_start_col, llen, mutableB.get_num_cols());
		run_part_dense(ins);
	}
	mutableA.restore_size(orig_A);
	mutableB.restore_size(orig_B);
}

void multiply_wide_op::run_part_dense(
		const std::vector<detail::local_matrix_store::const_ptr> &ins) const
{
	int thread_id = detail::mem_thread_pool::get_curr_thread_id();

	multiply_wide_op *mutable_this = const_cast<multiply_wide_op *>(this);
	detail::local_matrix_store::const_ptr Astore = ins[0];
	const void *Amat = Astore->get_raw_arr();
	if (Amat == NULL || Astore->store_layout() != Alayout) {
		if (Abufs[thread_id] == NULL
				|| Astore->get_num_rows() != Abufs[thread_id]->get_num_rows()
				|| Astore->get_num_cols() != Abufs[thread_id]->get_num_cols()) {
			if (Alayout == matrix_layout_t::L_ROW)
				mutable_this->Abufs[thread_id] = detail::local_matrix_store::ptr(
						new fm::detail::local_buf_row_matrix_store(0, 0,
							Astore->get_num_rows(), Astore->get_num_cols(),
							Astore->get_type(), -1));
			else
				mutable_this->Abufs[thread_id] = detail::local_matrix_store::ptr(
						new fm::detail::local_buf_col_matrix_store(0, 0,
							Astore->get_num_rows(), Astore->get_num_cols(),
							Astore->get_type(), -1));
		}
		Abufs[thread_id]->copy_from(*Astore);
		Amat = Abufs[thread_id]->get_raw_arr();
	}
	assert(Amat);

	detail::local_matrix_store::const_ptr Bstore = ins[1];
	const void *Bmat = Bstore->get_raw_arr();
	if (Bmat == NULL || Bstore->store_layout() != Blayout) {
		if (Bbufs[thread_id] == NULL
				|| Bstore->get_num_rows() != Bbufs[thread_id]->get_num_rows()
				|| Bstore->get_num_cols() != Bbufs[thread_id]->get_num_cols()) {
			if (Blayout == matrix_layout_t::L_COL)
				mutable_this->Bbufs[thread_id] = detail::local_matrix_store::ptr(
						new fm::detail::local_buf_col_matrix_store(0, 0,
							Bstore->get_num_rows(), Bstore->get_num_cols(),
							Bstore->get_type(), -1));
			else
				mutable_this->Bbufs[thread_id] = detail::local_matrix_store::ptr(
						new fm::detail::local_buf_row_matrix_store(0, 0,
							Bstore->get_num_rows(), Bstore->get_num_cols(),
							Bstore->get_type(), -1));
		}
		Bbufs[thread_id]->copy_from(*Bstore);
		Bmat = Bbufs[thread_id]->get_raw_arr();
	}
	assert(Bmat);

	if (res_bufs[thread_id] == NULL) {
		if (Blayout == matrix_layout_t::L_COL)
			mutable_this->tmp_bufs[thread_id] = detail::local_matrix_store::ptr(
					new fm::detail::local_buf_col_matrix_store(0, 0,
						out_num_rows, out_num_cols, get_output_type(), -1));
		else
			mutable_this->tmp_bufs[thread_id] = detail::local_matrix_store::ptr(
					new fm::detail::local_buf_row_matrix_store(0, 0,
						out_num_rows, out_num_cols, get_output_type(), -1));
		tmp_bufs[thread_id]->reset_data();
		mutable_this->res_bufs[thread_id] = matmul_accumulator::create(
				out_num_rows, out_num_cols, Blayout, get_output_type());
	}
	assert(tmp_bufs[thread_id]->store_layout() == Blayout);
	void *tmp_mat = tmp_bufs[thread_id]->get_raw_arr();
	std::pair<size_t, size_t> Asize, Bsize;
	if (require_trans) {
		assert(Alayout != Blayout);
		Asize.first = Astore->get_num_cols();
		Asize.second = Astore->get_num_rows();
	}
	else {
		assert(Alayout == Blayout);
		Asize.first = Astore->get_num_rows();
		Asize.second = Astore->get_num_cols();
	}
	Bsize.first = Bstore->get_num_rows();
	Bsize.second = Bstore->get_num_cols();
	assert(out_num_rows == Asize.first && out_num_cols == Bsize.second);
	if (get_output_type() == get_scalar_type<double>()) {
		const double *t_Amat = reinterpret_cast<const double *>(Amat);
		const double *t_Bmat = reinterpret_cast<const double *>(Bmat);
		double *t_tmp_mat = reinterpret_cast<double *>(tmp_mat);
		if (Blayout == matrix_layout_t::L_COL)
			wide_gemm_col<double>(Asize, t_Amat, Bsize, t_Bmat, t_tmp_mat,
					out_num_rows);
		else
			wide_gemm_row<double>(Asize, t_Amat, Bsize, t_Bmat, t_tmp_mat,
					out_num_cols);
	}
	else {
		const float *t_Amat = reinterpret_cast<const float *>(Amat);
		const float *t_Bmat = reinterpret_cast<const float *>(Bmat);
		float *t_tmp_mat = reinterpret_cast<float *>(tmp_mat);
		if (Blayout == matrix_layout_t::L_COL)
			wide_gemm_col<float>(Asize, t_Amat, Bsize, t_Bmat, t_tmp_mat,
					out_num_rows);
		else
			wide_gemm_row<float>(Asize, t_Amat, Bsize, t_Bmat, t_tmp_mat,
					out_num_cols);
	}
	res_bufs[thread_id]->add_matrix(*tmp_bufs[thread_id]);
}

}

IPW_matrix_store::IPW_matrix_store(matrix_store::const_ptr left,
		matrix_store::const_ptr right, bulk_operate::const_ptr left_op,
		bulk_operate::const_ptr right_op,
		matrix_layout_t layout): virtual_matrix_store(
			left->get_num_rows(), left->get_num_cols(),
			left->is_in_mem() && right->is_in_mem(), left->get_type())
{
	this->left_mat = left;
	this->right_mat = right;

	size_t nthreads = detail::mem_thread_pool::get_global_num_threads();
	bool use_blas = left_op == NULL;
	if (use_blas) {
		assert(!left->is_sparse());
		this->left_op = NULL;
		assert(left->get_type() == get_scalar_type<double>()
				|| left->get_type() == get_scalar_type<float>());
		assert(right->get_type() == get_scalar_type<double>()
				|| right->get_type() == get_scalar_type<float>());
		this->right_op = bulk_operate::conv2ptr(get_type().get_basic_ops().get_add());

		matrix_layout_t required_layout;
		if (right->is_sparse())
			required_layout = matrix_layout_t::L_COL;
		// If both input matrices have the same data layout, it's easy.
		else if (left->store_layout() == right->store_layout())
			required_layout = left->store_layout();
		// If they are different, we convert the smaller matrix.
		else if (left->get_num_rows() * left->get_num_cols()
				> right->get_num_rows() * right->get_num_cols())
			required_layout = left->store_layout();
		else
			required_layout = right->store_layout();

		if (layout == matrix_layout_t::L_NONE)
			this->layout = required_layout;
		else
			this->layout = layout;

		assert(left->get_type() == right->get_type());
		portion_op = std::shared_ptr<portion_mapply_op>(
				new multiply_wide_op(nthreads, left->get_num_rows(),
					right->get_num_cols(), required_layout, left->get_type(),
					right->is_sparse()));
	}
	else {
		assert(!left->is_sparse());
		assert(!right->is_sparse());
		this->left_op = left_op;
		this->right_op = right_op;

		if (layout != matrix_layout_t::L_NONE)
			this->layout = layout;
		// If the left matrix is in col-major, we prefer the output matrix
		// is also in col-major. It helps computation in local matrices.
		else if (left->store_layout() == matrix_layout_t::L_COL)
			this->layout = matrix_layout_t::L_COL;
		else
			this->layout = matrix_layout_t::L_ROW;

		matrix_info info;
		info.num_rows = left->get_num_rows();
		info.num_cols = right->get_num_cols();
		info.layout = this->layout;
		portion_op = std::shared_ptr<portion_mapply_op>(new inner_prod_wide_op(
					left_op, right_op, info, nthreads));
	}
}

matrix_store::ptr IPW_matrix_store::get_combine_res() const
{
	// Aggregate the results from omp threads.
	detail::matrix_store::ptr res = std::static_pointer_cast<const combine_op>(
				portion_op)->get_combined_result();
	if (this->layout == res->store_layout())
		return res;
	else {
		// Otherwise, we need to convert the matrix layout.
		detail::mem_matrix_store::ptr tmp = detail::mem_matrix_store::create(
				res->get_num_rows(), res->get_num_cols(), layout,
				res->get_type(), -1);
		detail::local_matrix_store::const_ptr lres = res->get_portion(0);
		tmp->write_portion_async(lres, 0, 0);
		return tmp;
	}
}

bool IPW_matrix_store::has_materialized() const
{
	return std::static_pointer_cast<combine_op>(portion_op)->has_materialized();
}

void IPW_matrix_store::materialize_self() const
{
	if (!has_materialized()) {
		std::static_pointer_cast<combine_op>(portion_op)->set_require_trans(true);
		// This computes the partial aggregation result.
		std::vector<detail::matrix_store::const_ptr> ins(2);
		ins[0] = left_mat->transpose();
		ins[1] = right_mat;
		__mapply_portion(ins, portion_op, layout);
		std::static_pointer_cast<combine_op>(portion_op)->set_require_trans(false);
	}
}

matrix_store::const_ptr IPW_matrix_store::materialize(bool in_mem,
		int num_nodes) const
{
	materialize_self();
	return get_combine_res();
}

vec_store::const_ptr IPW_matrix_store::get_col_vec(off_t idx) const
{
	matrix_store::const_ptr ret = materialize(true, -1);
	return ret->get_col_vec(idx);
}

vec_store::const_ptr IPW_matrix_store::get_row_vec(off_t idx) const
{
	matrix_store::const_ptr ret = materialize(true, -1);
	return ret->get_row_vec(idx);
}

matrix_store::const_ptr IPW_matrix_store::get_cols(
		const std::vector<off_t> &idxs) const
{
	matrix_store::const_ptr ret = materialize(true, -1);
	return ret->get_cols(idxs);
}

matrix_store::const_ptr IPW_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	matrix_store::const_ptr ret = materialize(true, -1);
	return ret->get_rows(idxs);
}

static int get_node_id(const local_matrix_store &left,
		const local_matrix_store &right)
{
	// If both matrices are stored in NUMA memory, the portion must be
	// stored on the same NUMA node. Otherwise, we need to return
	// the node Id from the matrix stored in NUMA.
	if (left.get_node_id() < 0)
		return right.get_node_id();
	else
		return left.get_node_id();
}

namespace
{

/*
 * The role of these two matrices is to materialize the underlying local matrix
 * piece by piece so that we can keep data in the CPU cache when computing
 * aggregation.
 */

class lmaterialize_col_matrix_store: public lvirtual_col_matrix_store
{
	std::vector<local_matrix_store::const_ptr> parts;
	local_matrix_store &mutable_left_part;
	local_matrix_store &mutable_right_part;
	portion_mapply_op::const_ptr portion_op;
public:
	lmaterialize_col_matrix_store(local_matrix_store::const_ptr left_part,
			local_matrix_store::const_ptr right_part, const scalar_type &type,
			portion_mapply_op::const_ptr portion_op): lvirtual_col_matrix_store(
				left_part->get_global_start_row(),
				left_part->get_global_start_col(),
				left_part->get_num_rows(), left_part->get_num_cols(), type,
				fm::detail::get_node_id(*left_part, *right_part)),
			mutable_left_part(const_cast<local_matrix_store &>(*left_part)),
			mutable_right_part(const_cast<local_matrix_store &>(*right_part)) {
		this->portion_op = portion_op;
		parts.resize(2);
		parts[0] = left_part;
		parts[1] = right_part;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		assert(local_start_row == 0);
		assert(local_num_rows == mutable_left_part.get_num_rows());
		mutable_left_part.resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
		// We need to resize the portion of the right matrix accordingly.
		mutable_right_part.resize(local_start_col, 0, local_num_cols,
				mutable_right_part.get_num_cols());
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		mutable_left_part.reset_size();
		mutable_right_part.reset_size();
		local_matrix_store::reset_size();
	}

	using lvirtual_col_matrix_store::get_raw_arr;
	virtual const char *get_raw_arr() const {
		assert(0);
		return NULL;
	}

	using lvirtual_col_matrix_store::transpose;
	virtual matrix_store::const_ptr transpose() const {
		assert(0);
		return matrix_store::const_ptr();
	}

	using lvirtual_col_matrix_store::get_col;
	virtual const char *get_col(size_t col) const {
		assert(0);
		return NULL;
	}

	virtual local_matrix_store::const_ptr get_portion(
			size_t local_start_row, size_t local_start_col, size_t num_rows,
			size_t num_cols) const {
		assert(0);
		return local_matrix_store::const_ptr();
	}

	virtual void materialize_self() const {
		portion_op->run(parts);
	}
};

class lmaterialize_row_matrix_store: public lvirtual_row_matrix_store
{
	std::vector<local_matrix_store::const_ptr> parts;
	local_matrix_store &mutable_left_part;
	local_matrix_store &mutable_right_part;
	portion_mapply_op::const_ptr portion_op;
public:
	lmaterialize_row_matrix_store(local_matrix_store::const_ptr left_part,
			local_matrix_store::const_ptr right_part, const scalar_type &type,
			portion_mapply_op::const_ptr portion_op): lvirtual_row_matrix_store(
				left_part->get_global_start_row(),
				left_part->get_global_start_col(),
				left_part->get_num_rows(), left_part->get_num_cols(), type,
				fm::detail::get_node_id(*left_part, *right_part)),
			mutable_left_part(const_cast<local_matrix_store &>(*left_part)),
			mutable_right_part(const_cast<local_matrix_store &>(*right_part)) {
		this->portion_op = portion_op;
		parts.resize(2);
		parts[0] = left_part;
		parts[1] = right_part;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		assert(local_start_row == 0);
		assert(local_num_rows == mutable_left_part.get_num_rows());
		mutable_left_part.resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
		// We need to resize the portion of the right matrix accordingly.
		mutable_right_part.resize(local_start_col, 0, local_num_cols,
				mutable_right_part.get_num_cols());
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		mutable_left_part.reset_size();
		mutable_right_part.reset_size();
		local_matrix_store::reset_size();
	}

	using lvirtual_row_matrix_store::get_raw_arr;
	virtual const char *get_raw_arr() const {
		assert(0);
		return NULL;
	}

	using lvirtual_row_matrix_store::transpose;
	virtual matrix_store::const_ptr transpose() const {
		assert(0);
		return matrix_store::const_ptr();
	}

	using lvirtual_row_matrix_store::get_row;
	virtual const char *get_row(size_t row) const {
		assert(0);
		return NULL;
	}

	virtual local_matrix_store::const_ptr get_portion(
			size_t local_start_row, size_t local_start_col, size_t num_rows,
			size_t num_cols) const {
		assert(0);
		return local_matrix_store::const_ptr();
	}

	virtual void materialize_self() const {
		portion_op->run(parts);
	}
};

}

static local_matrix_store::const_ptr create_lmaterialize_matrix(
		local_matrix_store::const_ptr left_part,
		local_matrix_store::const_ptr right_part, const scalar_type &type,
		portion_mapply_op::const_ptr portion_op)
{
	if (left_part->store_layout() == matrix_layout_t::L_ROW)
		return local_matrix_store::const_ptr(new lmaterialize_row_matrix_store(
					left_part, right_part, type, portion_op));
	else
		return local_matrix_store::const_ptr(new lmaterialize_col_matrix_store(
					left_part, right_part, type, portion_op));
}

local_matrix_store::const_ptr IPW_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	assert(start_row == 0);
	assert(num_rows == left_mat->get_num_rows());
	local_matrix_store::const_ptr left_part = left_mat->get_portion(start_row,
			start_col, num_rows, num_cols);
	local_matrix_store::const_ptr right_part = right_mat->get_portion(
			start_col, 0, num_cols, right_mat->get_num_cols());
	assert(left_part->get_num_cols() == right_part->get_num_rows());
	return create_lmaterialize_matrix(left_part, right_part, get_type(),
			portion_op);
}

local_matrix_store::const_ptr IPW_matrix_store::get_portion(size_t id) const
{
	local_matrix_store::const_ptr left_part = left_mat->get_portion(id);
	local_matrix_store::const_ptr right_part = right_mat->get_portion(id);
	assert(left_part->get_num_cols() == right_part->get_num_rows());
	return create_lmaterialize_matrix(left_part, right_part, get_type(),
			portion_op);
}

namespace
{

class collect_portion_compute: public portion_compute
{
	size_t num_EM_parts;
	size_t num_reads;
	portion_compute::ptr orig_compute;
public:
	typedef std::shared_ptr<collect_portion_compute> ptr;

	collect_portion_compute(portion_compute::ptr orig_compute) {
		this->num_EM_parts = 0;
		this->num_reads = 0;
		this->orig_compute = orig_compute;
	}

	void set_EM_count(size_t num_EM_parts) {
		this->num_EM_parts = num_EM_parts;
	}

	virtual void run(char *buf, size_t size) {
		num_reads++;
		if (num_reads == num_EM_parts) {
			orig_compute->run(NULL, 0);
			// This only runs once.
			// Let's remove all user's portion compute to indicate that it has
			// been invoked.
			orig_compute = NULL;
		}
	}
};

}

async_cres_t IPW_matrix_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, std::shared_ptr<portion_compute> compute) const
{
	assert(start_row == 0);
	assert(num_rows == left_mat->get_num_rows());
	collect_portion_compute::ptr new_compute(new collect_portion_compute(
				compute));
	async_cres_t left_ret = left_mat->get_portion_async(start_row, start_col,
			num_rows, num_cols, new_compute);
	async_cres_t right_ret = right_mat->get_portion_async(start_col, 0,
			num_cols, right_mat->get_num_cols(), new_compute);
	assert(left_ret.second->get_num_cols() == right_ret.second->get_num_rows());
	if (left_ret.first && right_ret.first)
		return async_cres_t(true, create_lmaterialize_matrix(left_ret.second,
					right_ret.second, get_type(), portion_op));
	else {
		new_compute->set_EM_count(!left_ret.first + !right_ret.first);
		return async_cres_t(false, create_lmaterialize_matrix(left_ret.second,
					right_ret.second, get_type(), portion_op));
	}
}

matrix_store::const_ptr IPW_matrix_store::transpose() const
{
	// TODO do we need this?
	assert(0);
	return matrix_store::const_ptr();
}

std::vector<safs::io_interface::ptr> IPW_matrix_store::create_ios() const
{
	std::vector<safs::io_interface::ptr> ret;
	if (!left_mat->is_in_mem()) {
		const EM_object *obj = dynamic_cast<const EM_object *>(left_mat.get());
		std::vector<safs::io_interface::ptr> tmp = obj->create_ios();
		ret.insert(ret.end(), tmp.begin(), tmp.end());
	}
	if (!right_mat->is_in_mem()) {
		const EM_object *obj = dynamic_cast<const EM_object *>(right_mat.get());
		std::vector<safs::io_interface::ptr> tmp = obj->create_ios();
		ret.insert(ret.end(), tmp.begin(), tmp.end());
	}
	return ret;
}

std::unordered_map<size_t, size_t> IPW_matrix_store::get_underlying_mats() const
{
	std::unordered_map<size_t, size_t> final_res = left_mat->get_underlying_mats();
	std::unordered_map<size_t, size_t> right = right_mat->get_underlying_mats();
	for (auto it = right.begin(); it != right.end(); it++) {
		auto to_it = final_res.find(it->first);
		if (to_it == final_res.end())
			final_res.insert(std::pair<size_t, size_t>(it->first, it->second));
	}
	return final_res;
}

std::string IPW_matrix_store::get_name() const
{
	std::vector<matrix_store::const_ptr> mats(2);
	mats[0] = left_mat;
	mats[1] = right_mat;
	return portion_op->to_string(mats);
}

}

}
