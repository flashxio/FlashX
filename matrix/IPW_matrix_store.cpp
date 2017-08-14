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
#include "materialize.h"
#include "local_matrix_store.h"
#include "project_matrix_store.h"

namespace fm
{

namespace detail
{

namespace
{

static inline matrix_layout_t opposite_layout(matrix_layout_t layout)
{
	return layout == matrix_layout_t::L_ROW
		? matrix_layout_t::L_COL : matrix_layout_t::L_ROW;
}

struct matrix_info
{
	size_t num_rows;
	size_t num_cols;
	matrix_layout_t layout;
};

class combine_op: public detail::portion_mapply_op
{
	const size_t data_id;
public:
	combine_op(size_t num_rows, size_t num_cols,
			const scalar_type &type): detail::portion_mapply_op(
				num_rows, num_cols, type), data_id(matrix_store::mat_counter++) {
	}
	size_t get_data_id() const {
		return data_id;
	}
	virtual bool has_materialized() const = 0;
	virtual detail::mem_matrix_store::const_ptr get_combined_result() const = 0;
};

class inner_prod_wide_op: public combine_op
{
	bulk_operate::const_ptr left_op;
	bulk_operate::const_ptr right_op;
	matrix_info out_mat_info;
	matrix_layout_t left_layout;
	matrix_layout_t right_layout;
	std::vector<detail::local_matrix_store::ptr> local_ms;
	std::vector<detail::local_matrix_store::ptr> local_tmps;
public:
	inner_prod_wide_op(bulk_operate::const_ptr left_op,
			bulk_operate::const_ptr right_op, const matrix_info &out_mat_info,
			matrix_layout_t left_layout, matrix_layout_t right_layout,
			size_t num_threads): combine_op(0, 0, right_op->get_output_type()) {
		this->left_op = left_op;
		this->right_op = right_op;
		local_ms.resize(num_threads);
		local_tmps.resize(num_threads);
		this->out_mat_info = out_mat_info;
		this->left_layout = left_layout;
		this->right_layout = right_layout;
		if (left_layout == matrix_layout_t::L_COL)
			this->out_mat_info.layout = matrix_layout_t::L_COL;
		else
			this->out_mat_info.layout = matrix_layout_t::L_ROW;
	}

	virtual bool has_materialized() const {
		bool materialized = false;
		for (size_t i = 0; i < local_ms.size(); i++)
			if (local_ms[i])
				materialized = true;
		return materialized;
	}

	virtual detail::mem_matrix_store::const_ptr get_combined_result() const;

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const;

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		matrix_info info;
		info.num_rows = out_mat_info.num_cols;
		info.num_cols = out_mat_info.num_rows;
		return detail::portion_mapply_op::const_ptr(
				new inner_prod_wide_op(left_op, right_op, info,
					opposite_layout(right_layout), opposite_layout(left_layout),
					local_ms.size()));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 2);
		return std::string("inner_prod(") + mats[0]->get_name()
			+ "," + mats[1]->get_name() + ")";
	}
};

detail::mem_matrix_store::const_ptr inner_prod_wide_op::get_combined_result() const
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
			detail::mapply2(*local_res, *local_ms[j], *right_op,
					part_dim_t::PART_NONE, *local_res);
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

	assert(ins[0]->get_num_rows() == ins[1]->get_num_rows());
	// We always transpose the left matrix to make it a tall matrix.
	detail::local_matrix_store::const_ptr store
		= std::static_pointer_cast<const detail::local_matrix_store>(
				ins[0]->transpose());
	if (is_first)
		detail::inner_prod_wide(*store, *ins[1], *left_op, *right_op, *local_m);
	else {
		detail::inner_prod_wide(*store, *ins[1], *left_op, *right_op, *local_tmp);
		// We don't need to further partition the result matrix when
		// summing them up.
		mapply2(*local_m, *local_tmp, *right_op, part_dim_t::PART_NONE,
				*local_m);
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

class part_mul_res
{
	std::vector<matmul_accumulator::ptr> res_bufs;
	size_t num_rows;
	size_t num_cols;
	matrix_layout_t layout;
	const scalar_type &type;
public:
	typedef std::shared_ptr<part_mul_res> ptr;

	part_mul_res(int num_threads, size_t num_rows, size_t num_cols,
			matrix_layout_t layout, const scalar_type &_type): type(_type) {
		res_bufs.resize(num_threads);
		this->num_rows = num_rows;
		this->num_cols = num_cols;
		this->layout = layout;
	}

	void acc_part_res(int thread, const local_matrix_store &part) {
		if (res_bufs[thread] == NULL)
			res_bufs[thread] = matmul_accumulator::create(num_rows, num_cols,
					layout, type);
		if (layout == part.store_layout()) {
			assert(num_rows == part.get_num_rows());
			assert(num_cols == part.get_num_cols());
			res_bufs[thread]->add_matrix(part);
		}
		else {
			assert(num_rows == part.get_num_cols());
			assert(num_cols == part.get_num_rows());
			local_matrix_store::const_ptr tpart
				= std::static_pointer_cast<const local_matrix_store>(
						part.transpose());
			res_bufs[thread]->add_matrix(*tpart);
		}
	}

	bool has_materialized() const {
		bool materialized = false;
		for (size_t i = 0; i < res_bufs.size(); i++)
			if (res_bufs[i])
				materialized = true;
		return materialized;
	}

	detail::mem_matrix_store::ptr get_combined_result() const {
		std::vector<matmul_accumulator::ptr> non_empty;
		for (size_t i = 0; i < res_bufs.size(); i++) {
			if (res_bufs[i])
				non_empty.push_back(res_bufs[i]);
		}
		assert(non_empty.size() > 0);
		return non_empty[0]->combine(non_empty);
	}
};

class multiply_wide_op: public combine_op
{
	// The number of times we have accumulated the computation results
	// on a portion in the tmp_buf. This is only useful for sparse matrix
	// multiplication.
	std::vector<size_t> num_tmp_accs;
	std::vector<detail::local_matrix_store::ptr> tmp_bufs;
	part_mul_res::ptr part_res;
	bool is_sparse;
	// The output matrix will be symmetric.
	bool is_sym;
	size_t out_num_rows;
	size_t out_num_cols;
	matrix_layout_t Alayout;
	matrix_layout_t Blayout;
public:
	multiply_wide_op(size_t num_threads, size_t out_num_rows, size_t out_num_cols,
			matrix_layout_t required_layout, const scalar_type &type,
			bool is_sparse, bool is_sym): combine_op(0, 0, type) {
		tmp_bufs.resize(num_threads);
		this->out_num_rows = out_num_rows;
		this->out_num_cols = out_num_cols;
		Alayout = required_layout;
		Blayout = required_layout;
		this->is_sparse = is_sparse;
		this->num_tmp_accs.resize(num_threads);
		this->is_sym = is_sym;
		this->part_res = part_mul_res::ptr(new part_mul_res(num_threads,
					out_num_rows, out_num_cols, Blayout, type));
	}

	virtual bool has_materialized() const {
		bool ret = part_res->has_materialized();
		// If this is a sparse matrix, we should check num_tmp_accs
		// because some results may be buffered in tmp_bufs and haven't
		// been accumulated to num_tmp_accs yet.
		if (!ret && is_sparse) {
			for (size_t i = 0; i < num_tmp_accs.size(); i++)
				if (num_tmp_accs[i] > 0)
					return true;
		}
		return ret;
	}

	virtual detail::mem_matrix_store::const_ptr get_combined_result() const;

	void run_part_dense(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const;
	void run_part_sparse(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const;

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const;

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		multiply_wide_op *ret = new multiply_wide_op(*this);
		for (size_t i = 0; i < ret->num_tmp_accs.size(); i++) {
			ret->num_tmp_accs[i] = 0;
			ret->tmp_bufs[i] = NULL;
		}
		ret->out_num_rows = out_num_cols;
		ret->out_num_cols = out_num_rows;
		ret->Alayout = opposite_layout(Blayout);
		ret->Blayout = opposite_layout(Alayout);
		// TODO currently, we only optimize the multiplication of a dense
		// matrix with a sparse matrix. If we transpose the operation,
		// it becomes the multiplication of a sparse matrix and a dense
		// matrix and our optimization wont' work any more. So we disable
		// the optimization here after transpose.
		ret->is_sparse = false;
		return std::shared_ptr<portion_mapply_op>(ret);
	}

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 2);
		return std::string("(") + (mats[0]->get_name()
					+ "*") + mats[1]->get_name() + std::string(")");
	}
};

detail::mem_matrix_store::const_ptr multiply_wide_op::get_combined_result() const
{
	multiply_wide_op *mutable_this = const_cast<multiply_wide_op *>(this);
	for (size_t i = 0; i < num_tmp_accs.size(); i++) {
		// We should accumulate the results in the tmp buffer.
		if (num_tmp_accs[i] > 0) {
			part_res->acc_part_res(i, *tmp_bufs[i]);
			tmp_bufs[i]->reset_data();
			mutable_this->num_tmp_accs[i] = 0;
		}
	}
	detail::mem_matrix_store::ptr ret = part_res->get_combined_result();
	if (is_sym)
		// For self-crossprod, we store data in the upper triangle,
		// we need to copy the data to the lower triangle.
		ret->symmetrize(true);
	// The combined result should have the same layout as B.
	// If not, the partial result was created for the transpose of the matrix.
	if (ret->store_layout() != Blayout)
		return std::static_pointer_cast<const detail::mem_matrix_store>(
				ret->transpose());
	else
		return ret;
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

template<class T>
void wide_syrk_col(const std::pair<size_t, size_t> &Asize, const T *Amat,
		T *res_mat)
{
	assert(0);
}

template<>
void wide_syrk_col<double>(const std::pair<size_t, size_t> &Asize,
		const double *Amat, double *res_mat)
{
	size_t n = Asize.first;
	size_t k = Asize.second;
	size_t lda = Asize.first;
	size_t ldc = Asize.first;
	cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans, n, k, 1, Amat, lda,
			0, res_mat, ldc);
}

template<>
void wide_syrk_col<float>(const std::pair<size_t, size_t> &Asize,
		const float *Amat, float *res_mat)
{
	size_t n = Asize.first;
	size_t k = Asize.second;
	size_t lda = Asize.first;
	size_t ldc = Asize.first;
	cblas_ssyrk(CblasColMajor, CblasUpper, CblasNoTrans, n, k, 1, Amat, lda,
			0, res_mat, ldc);
}

template<class T>
void wide_syrk_row(const std::pair<size_t, size_t> &Asize, const T *Amat,
		T *res_mat)
{
	assert(0);
}

template<>
void wide_syrk_row<double>(const std::pair<size_t, size_t> &Asize,
		const double *Amat, double *res_mat)
{
	size_t n = Asize.first;
	size_t k = Asize.second;
	size_t lda = Asize.second;
	size_t ldc = Asize.first;
	cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans, n, k, 1, Amat, lda,
			0, res_mat, ldc);
}

template<>
void wide_syrk_row<float>(const std::pair<size_t, size_t> &Asize,
		const float *Amat, float *res_mat)
{
	size_t n = Asize.first;
	size_t k = Asize.second;
	size_t lda = Asize.second;
	size_t ldc = Asize.first;
	cblas_ssyrk(CblasRowMajor, CblasUpper, CblasNoTrans, n, k, 1, Amat, lda,
			0, res_mat, ldc);
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
 * This multiplies a dense matrix t(A) with a sparse matrix B.
 */
template<class T>
void multiply_sparse_trans(const detail::local_row_matrix_store &Astore,
		const detail::lsparse_row_matrix_store &Bstore,
		detail::local_col_matrix_store &Cstore)
{
	assert(get_scalar_type<T>() == Bstore.get_type());
	std::vector<sparse_project_matrix_store::nz_idx> Bidxs;
	const char * _Brows = Bstore.get_rows_nnz(0, Bstore.get_num_rows(), Bidxs);
	// If the sparse submatrix doesn't have non-zero values.
	if (_Brows == NULL)
		return;
	const T *Brows = reinterpret_cast<const T *>(_Brows);
	for (auto it = Bidxs.begin(); it != Bidxs.end(); it++) {
		T B = Brows[it - Bidxs.begin()];
		const T *Arow = reinterpret_cast<const T *>(Astore.get_row(it->row_idx));
		T *dst_col = reinterpret_cast<T *>(Cstore.get_col(it->col_idx));
		cblas_axpy<T>(Astore.get_num_cols(), B, Arow, dst_col);
	}
}

void multiply_wide_op::run_part_sparse(
		const std::vector<detail::local_matrix_store::const_ptr> &ins) const
{
	// We always make the left matrix a dense matrix, so there must be
	// two matrices as input.
	assert(ins.size() == 2);
	int thread_id = detail::mem_thread_pool::get_curr_thread_id();
	multiply_wide_op *mutable_this = const_cast<multiply_wide_op *>(this);
	detail::local_matrix_store::const_ptr left = ins[0];
	assert(ins[1]->store_layout() == matrix_layout_t::L_ROW);
	const detail::lsparse_row_matrix_store &Bstore
		= static_cast<const detail::lsparse_row_matrix_store &>(*ins[1]);

	if (tmp_bufs[thread_id] == NULL) {
		mutable_this->tmp_bufs[thread_id] = detail::local_matrix_store::ptr(
				new local_buf_col_matrix_store(0, 0,
					out_num_rows, out_num_cols, get_output_type(), -1));
		tmp_bufs[thread_id]->reset_data();
	}
	assert(tmp_bufs[thread_id]->store_layout() == Blayout);

	local_col_matrix_store &tmp_col_buf = static_cast<local_col_matrix_store &>(
			*tmp_bufs[thread_id]);
	if (get_output_type() == get_scalar_type<double>()) {
		assert(left->store_layout() == matrix_layout_t::L_ROW);
		const detail::local_row_matrix_store &Astore
			= static_cast<const detail::local_row_matrix_store &>(*left);
		multiply_sparse_trans<double>(Astore, Bstore, tmp_col_buf);
	}
	else {
		assert(get_output_type() == get_scalar_type<float>());
		assert(left->store_layout() == matrix_layout_t::L_ROW);
		const detail::local_row_matrix_store &Astore
			= static_cast<const detail::local_row_matrix_store &>(*left);
		multiply_sparse_trans<float>(Astore, Bstore, tmp_col_buf);
	}

	// We should accumulate results in res_bufs, which offers high precision
	// for accumulation. This is a tradeoff between precision and computation
	// overhead.
	size_t thres = mem_matrix_store::CHUNK_SIZE / std::max(
			ins[0]->get_num_rows(), ins[0]->get_num_cols());
	mutable_this->num_tmp_accs[thread_id]++;
	if (num_tmp_accs[thread_id] > thres) {
		part_res->acc_part_res(thread_id, *tmp_bufs[thread_id]);
		tmp_bufs[thread_id]->reset_data();
		mutable_this->num_tmp_accs[thread_id] = 0;
	}
}

void multiply_wide_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins) const
{
	assert(ins.size() == 2 || ins.size() == 1);
	// If one of the matrix is a sparse matrix, some of its portions might
	// be empty.
	for (size_t i = 0; i < ins.size(); i++)
		if (ins[i] == NULL)
			return;

	/*
	 * We can't determine the partition size for dense matrix multiplication
	 * in the same way as other computations, because dense matrix multiplication
	 * is implemented with BLAS. To get better performance from BLAS, we should
	 * use larger partition sizes. For crossprod, we even use a larger block
	 * size to partition the input matrices. As such, the short dimension
	 * is usually larger. We want to the partition dimension length have
	 * a larger size (at least 4 times larger than the short dimension).
	 */
	size_t LONG_DIM_LEN;
	size_t short_dim_len;
	if (ins.size() == 2) {
		LONG_DIM_LEN = get_long_dim_len(*ins[0], *ins[1]);
		size_t short_dim1 = std::min(ins[0]->get_num_rows(),
				ins[0]->get_num_cols());
		size_t short_dim2 = std::min(ins[1]->get_num_rows(),
				ins[1]->get_num_cols());
		short_dim_len = std::max(short_dim1, short_dim2);
	}
	else {
		LONG_DIM_LEN = get_long_dim_len(*ins[0]);
		short_dim_len = std::min(ins[0]->get_num_rows(),
				ins[0]->get_num_cols());
	}
	if (!is_sparse)
		LONG_DIM_LEN = std::max(short_dim_len * 4, LONG_DIM_LEN);

	// ins[0] stores the transpose of the left matrix, so we should
	// check the number of rows.
	size_t long_dim = ins[0]->get_num_rows();
	if (long_dim <= LONG_DIM_LEN) {
		if (is_sparse)
			run_part_sparse(ins);
		else
			run_part_dense(ins);
		return;
	}

	local_matrix_store::exposed_area orig_A = ins[0]->get_exposed_area();
	local_matrix_store::exposed_area orig_B;
	if (ins.size() == 2)
		orig_B = ins[1]->get_exposed_area();
	local_matrix_store &mutableA = const_cast<local_matrix_store &>(*ins[0]);
	bool resize_success = true;
	for (size_t row_idx = 0; row_idx < long_dim; row_idx += LONG_DIM_LEN) {
		size_t llen = std::min(long_dim - row_idx, LONG_DIM_LEN);
		resize_success = mutableA.resize(orig_A.local_start_row + row_idx,
				orig_A.local_start_col, llen, mutableA.get_num_cols());
		if (resize_success && ins.size() == 2)
			resize_success = const_cast<local_matrix_store &>(*ins[1]).resize(
					orig_B.local_start_row + row_idx, orig_B.local_start_col,
					llen, ins[1]->get_num_cols());
		// If we resize both matrices, we perform computation on the resized
		// matrix.
		if (resize_success) {
			if (is_sparse)
				run_part_sparse(ins);
			else
				run_part_dense(ins);
		}
		else
			break;
	}
	mutableA.restore_size(orig_A);
	if (ins.size() == 2)
		const_cast<local_matrix_store &>(*ins[1]).restore_size(orig_B);
	// If the resize on one of the matrices failed, we perform the computation
	// here.
	if (!resize_success) {
		if (is_sparse)
			run_part_sparse(ins);
		else
			run_part_dense(ins);
	}
}

void multiply_wide_op::run_part_dense(
		const std::vector<detail::local_matrix_store::const_ptr> &ins) const
{
	int thread_id = detail::mem_thread_pool::get_curr_thread_id();

	detail::local_matrix_store::ptr Atmp, Btmp;
	multiply_wide_op *mutable_this = const_cast<multiply_wide_op *>(this);
	detail::local_matrix_store::const_ptr Astore = ins[0];
	const void *Amat = Astore->get_raw_arr();
	// We will transpose A later. So, if A has the expected layout, we actually
	// need to convert the layout of A.
	if (Amat == NULL || Astore->store_layout() == Alayout) {
		// If we expect a col-major layout, we should have the data stored in
		// row-major order.
		if (Alayout == matrix_layout_t::L_COL)
			Atmp = detail::local_matrix_store::ptr(
					new fm::detail::local_buf_row_matrix_store(0, 0,
						Astore->get_num_rows(), Astore->get_num_cols(),
						Astore->get_type(), -1));
		else
			Atmp = detail::local_matrix_store::ptr(
					new fm::detail::local_buf_col_matrix_store(0, 0,
						Astore->get_num_rows(), Astore->get_num_cols(),
						Astore->get_type(), -1));
		Atmp->copy_from(*Astore);
		Amat = Atmp->get_raw_arr();
	}
	assert(Amat);

	detail::local_matrix_store::const_ptr Bstore;
	const void *Bmat = NULL;
	if (ins.size() == 2) {
		Bstore = ins[1];
		Bmat = Bstore->get_raw_arr();
		if (Bmat == NULL || Bstore->store_layout() != Blayout) {
			if (Blayout == matrix_layout_t::L_COL)
				Btmp = detail::local_matrix_store::ptr(
						new fm::detail::local_buf_col_matrix_store(0, 0,
							Bstore->get_num_rows(), Bstore->get_num_cols(),
							Bstore->get_type(), -1));
			else
				Btmp = detail::local_matrix_store::ptr(
						new fm::detail::local_buf_row_matrix_store(0, 0,
							Bstore->get_num_rows(), Bstore->get_num_cols(),
							Bstore->get_type(), -1));
			Btmp->copy_from(*Bstore);
			Bmat = Btmp->get_raw_arr();
		}
	}

	if (tmp_bufs[thread_id] == NULL) {
		if (Blayout == matrix_layout_t::L_COL)
			mutable_this->tmp_bufs[thread_id] = detail::local_matrix_store::ptr(
					new fm::detail::local_buf_col_matrix_store(0, 0,
						out_num_rows, out_num_cols, get_output_type(), -1));
		else
			mutable_this->tmp_bufs[thread_id] = detail::local_matrix_store::ptr(
					new fm::detail::local_buf_row_matrix_store(0, 0,
						out_num_rows, out_num_cols, get_output_type(), -1));
		tmp_bufs[thread_id]->reset_data();
	}
	assert(tmp_bufs[thread_id]->store_layout() == Alayout);
	void *tmp_mat = tmp_bufs[thread_id]->get_raw_arr();
	// We have two matrices as input.
	if (Bstore) {
		std::pair<size_t, size_t> Asize, Bsize;
		Asize.first = Astore->get_num_cols();
		Asize.second = Astore->get_num_rows();
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
	}
	// There is only one matrix as input. The right matrix is the transpose
	// of the first matrix.
	else {
		std::pair<size_t, size_t> Asize;
		Asize.first = Astore->get_num_cols();
		Asize.second = Astore->get_num_rows();
		assert(out_num_rows == Asize.first && out_num_cols == Asize.first);
		if (get_output_type() == get_scalar_type<double>()) {
			const double *t_Amat = reinterpret_cast<const double *>(Amat);
			double *t_tmp_mat = reinterpret_cast<double *>(tmp_mat);
			if (Astore->store_layout() == matrix_layout_t::L_ROW)
				wide_syrk_col<double>(Asize, t_Amat, t_tmp_mat);
			else
				wide_syrk_row<double>(Asize, t_Amat, t_tmp_mat);
		}
		else {
			const float *t_Amat = reinterpret_cast<const float *>(Amat);
			float *t_tmp_mat = reinterpret_cast<float *>(tmp_mat);
			if (Astore->store_layout() == matrix_layout_t::L_ROW)
				wide_syrk_col<float>(Asize, t_Amat, t_tmp_mat);
			else
				wide_syrk_row<float>(Asize, t_Amat, t_tmp_mat);
		}
	}
	part_res->acc_part_res(thread_id, *tmp_bufs[thread_id]);
	for (size_t i = 0; i < ins.size(); i++)
		ins[i]->complete();
}

}

IPW_matrix_store::IPW_matrix_store(matrix_store::const_ptr left,
		matrix_store::const_ptr right, bulk_operate::const_ptr left_op,
		bulk_operate::const_ptr right_op, matrix_layout_t layout): sink_store(
			left->get_num_rows(), right->get_num_cols(),
			left->is_in_mem() && right->is_in_mem(), left->get_type())
{
	this->left_mat = left;
	this->right_mat = right;

	size_t nthreads = detail::mem_thread_pool::get_global_num_threads();
	bool use_blas = left_op == NULL;
	if (use_blas && (left->get_type() == get_scalar_type<double>()
				|| left->get_type() == get_scalar_type<float>())
			&& (right->get_type() == get_scalar_type<double>()
				|| right->get_type() == get_scalar_type<float>())) {
		// If the two matrices share data, this operation is self-crossprod.
		// We will optimize for this case because self-crossprod is a very
		// common operation in machine learning.
		if (this->left_mat->share_data(*this->right_mat))
			this->right_mat = NULL;

		this->left_op = NULL;
		this->right_op = bulk_operate::conv2ptr(get_type().get_basic_ops().get_add());

		matrix_layout_t required_layout;
		if (right->is_sparse())
			required_layout = matrix_layout_t::L_COL;
		// If both input matrices have the same data layout, it's easy.
		else if (left->store_layout() == right->store_layout())
			required_layout = left->store_layout();
		// If they are different, we convert the smaller matrix.
		else if (left->get_num_rows() * left->get_num_cols()
				>= right->get_num_rows() * right->get_num_cols())
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
					// We only get benefit if the right matrix is sparse
					// and it's stored in row major.
					right->is_sparse()
					&& right->store_layout() == matrix_layout_t::L_ROW,
					// If we don't use right_mat, we'll output a symmetric matrix.
					this->right_mat == NULL));
	}
	else {
		this->right_mat = right;

		if (left_op) {
			this->left_op = left_op;
			this->right_op = right_op;
		}
		else {
			left_op = bulk_operate::conv2ptr(
					left->get_type().get_basic_ops().get_multiply());
			right_op = bulk_operate::conv2ptr(
					left->get_type().get_basic_ops().get_add());
			assert(left->get_type() == right->get_type());
			this->left_op = left_op;
			this->right_op = right_op;
		}
		assert(left_op);
		assert(right_op);

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
		portion_op = std::shared_ptr<portion_mapply_op>(new inner_prod_wide_op(
					left_op, right_op, info, left->store_layout(),
					right->store_layout(), nthreads));
	}
	this->underlying = get_underlying_mats();
}

IPW_matrix_store::IPW_matrix_store(matrix_store::const_ptr left,
		matrix_store::const_ptr right, bulk_operate::const_ptr left_op,
		bulk_operate::const_ptr right_op, matrix_layout_t layout,
		std::shared_ptr<const portion_mapply_op> portion_op): sink_store(
			left->get_num_rows(), right->get_num_cols(),
			left->is_in_mem() && right->is_in_mem(), left->get_type())
{
	this->left_mat = left;
	this->right_mat = right;
	this->left_op = left_op;
	this->right_op = right_op;
	this->portion_op = portion_op;
	this->layout = layout;
	this->underlying = get_underlying_mats();
}

matrix_store::const_ptr IPW_matrix_store::get_combine_res() const
{
	// Aggregate the results from omp threads.
	detail::matrix_store::const_ptr res = std::static_pointer_cast<const combine_op>(
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

size_t IPW_matrix_store::get_data_id() const
{
	return std::static_pointer_cast<const combine_op>(portion_op)->get_data_id();
}

bool IPW_matrix_store::has_materialized() const
{
	return std::static_pointer_cast<const combine_op>(portion_op)->has_materialized();
}

void IPW_matrix_store::materialize_self() const
{
	if (!has_materialized()) {
		// This computes the partial aggregation result.
		std::vector<detail::matrix_store::const_ptr> ins;
		ins.push_back(left_mat->transpose());
		if (right_mat)
			ins.push_back(right_mat);
		__mapply_portion(ins, portion_op, layout);
	}
}

matrix_store::const_ptr IPW_matrix_store::materialize(bool in_mem,
		int num_nodes) const
{
	materialize_self();
	return get_combine_res();
}

std::unordered_map<size_t, size_t> IPW_matrix_store::get_underlying_mats() const
{
	if (has_materialized())
		return std::unordered_map<size_t, size_t>();
	if (!this->underlying.empty())
		return this->underlying;

	std::unordered_map<size_t, size_t> final_res = left_mat->get_underlying_mats();
	if (right_mat) {
		std::unordered_map<size_t, size_t> right
			= right_mat->get_underlying_mats();
		for (auto it = right.begin(); it != right.end(); it++) {
			auto to_it = final_res.find(it->first);
			if (to_it == final_res.end())
				final_res.insert(std::pair<size_t, size_t>(it->first,
							it->second));
		}
	}
	return final_res;
}

static int get_node_id(local_matrix_store::const_ptr left,
		local_matrix_store::const_ptr right)
{
	// If both matrices are stored in NUMA memory, the portion must be
	// stored on the same NUMA node. Otherwise, we need to return
	// the node Id from the matrix stored in NUMA.
	if (left->get_node_id() < 0 && right)
		return right->get_node_id();
	else
		return left->get_node_id();
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
	portion_mapply_op::const_ptr portion_op;

	local_matrix_store &get_mutable_part(size_t i) {
		return const_cast<local_matrix_store &>(*parts[i]);
	}
public:
	lmaterialize_col_matrix_store(local_matrix_store::const_ptr left_part,
			local_matrix_store::const_ptr right_part, const scalar_type &type,
			portion_mapply_op::const_ptr portion_op): lvirtual_col_matrix_store(
				left_part->get_global_start_row(),
				left_part->get_global_start_col(),
				left_part->get_num_rows(), left_part->get_num_cols(), type,
				fm::detail::get_node_id(left_part, right_part)) {
		this->portion_op = portion_op;
		parts.push_back(left_part);
		if (right_part)
			parts.push_back(right_part);
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		assert(local_start_col == 0);
		assert(local_num_cols == parts[0]->get_num_cols());

		local_matrix_store::exposed_area orig_left = parts[0]->get_exposed_area();
		bool success = get_mutable_part(0).resize(local_start_row,
				local_start_col, local_num_rows, local_num_cols);
		if (!success)
			return false;

		// We need to resize the portion of the right matrix accordingly.
		if (parts.size() > 1) {
			success = get_mutable_part(1).resize(local_start_row, 0,
					local_num_rows, parts[1]->get_num_cols());
			if (!success) {
				get_mutable_part(0).restore_size(orig_left);
				return false;
			}
		}
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		for (size_t i = 0; i < parts.size(); i++)
			get_mutable_part(i).reset_size();
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
	portion_mapply_op::const_ptr portion_op;

	local_matrix_store &get_mutable_part(size_t i) {
		return const_cast<local_matrix_store &>(*parts[i]);
	}
public:
	lmaterialize_row_matrix_store(local_matrix_store::const_ptr left_part,
			local_matrix_store::const_ptr right_part, const scalar_type &type,
			portion_mapply_op::const_ptr portion_op): lvirtual_row_matrix_store(
				left_part->get_global_start_row(),
				left_part->get_global_start_col(),
				left_part->get_num_rows(), left_part->get_num_cols(), type,
				fm::detail::get_node_id(left_part, right_part)) {
		this->portion_op = portion_op;
		parts.push_back(left_part);
		if (right_part)
			parts.push_back(right_part);
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		assert(local_start_col == 0);
		assert(local_num_cols == parts[0]->get_num_cols());
		local_matrix_store::exposed_area orig_left = parts[0]->get_exposed_area();
		bool success = get_mutable_part(0).resize(local_start_row,
				local_start_col, local_num_rows, local_num_cols);
		if (!success)
			return false;

		// We need to resize the portion of the right matrix accordingly.
		if (parts.size() > 1) {
			success = get_mutable_part(1).resize(local_start_row, 0,
					local_num_rows, parts[1]->get_num_cols());
			if (!success) {
				get_mutable_part(0).restore_size(orig_left);
				return false;
			}
		}
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		for (size_t i = 0; i < parts.size(); i++)
			get_mutable_part(i).reset_size();
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

matrix_store::const_ptr IPW_matrix_store::transpose() const
{
	if (has_materialized())
		return get_combine_res()->transpose();

	if (right_mat) {
		matrix_store::const_ptr tleft = right_mat->transpose();
		matrix_store::const_ptr tright = left_mat->transpose();
		assert(layout == matrix_layout_t::L_ROW || layout == matrix_layout_t::L_COL);
		matrix_layout_t tlayout;
		if (layout == matrix_layout_t::L_ROW)
			tlayout = matrix_layout_t::L_COL;
		else
			tlayout = matrix_layout_t::L_ROW;
		return matrix_store::const_ptr(new IPW_matrix_store(tleft, tright,
					left_op, right_op, tlayout, portion_op->transpose()));
	}
	else
		return matrix_store::const_ptr(new IPW_matrix_store(*this));
}

std::string IPW_matrix_store::get_name() const
{
	std::vector<matrix_store::const_ptr> mats(2);
	mats[0] = left_mat;
	if (right_mat)
		mats[1] = right_mat;
	else
		mats[1] = left_mat;
	return portion_op->to_string(mats);
}

class IPW_compute_store: public sink_compute_store, public EM_object
{
	matrix_store::const_ptr left_mat;
	matrix_store::const_ptr right_mat;
	bulk_operate::const_ptr left_op;
	bulk_operate::const_ptr right_op;
	std::shared_ptr<const combine_op> portion_op;
	matrix_layout_t layout;
public:
	IPW_compute_store(matrix_store::const_ptr left_mat,
			matrix_store::const_ptr right_mat, bulk_operate::const_ptr left_op,
			bulk_operate::const_ptr right_op,
			portion_mapply_op::const_ptr portion_op, matrix_layout_t layout);
	virtual size_t get_data_id() const {
		return portion_op->get_data_id();
	}
	using virtual_matrix_store::get_portion;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t id) const;
	using virtual_matrix_store::get_portion_async;
	virtual async_cres_t get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, std::shared_ptr<portion_compute> compute) const;

	virtual int get_portion_node_id(size_t id) const {
		// If both matrices are stored in NUMA memory, the portion must be
		// stored on the same NUMA node. Otherwise, we need to return
		// the node Id from the matrix stored in NUMA.
		if (left_mat->get_num_nodes() > 0)
			return left_mat->get_portion_node_id(id);
		else if (right_mat)
			return right_mat->get_portion_node_id(id);
		else
			return -1;
	}

	virtual std::pair<size_t, size_t> get_portion_size() const {
		if (right_mat)
			assert(left_mat->get_portion_size().first
					== right_mat->get_portion_size().first);
		return left_mat->get_portion_size();
	}

	virtual int get_num_nodes() const {
		if (left_mat->get_num_nodes() > 0)
			return left_mat->get_num_nodes();
		else if (right_mat)
			return right_mat->get_num_nodes();
		else
			return -1;
	}

	virtual matrix_layout_t store_layout() const {
		// TODO what is the right layout?
		return layout;
	}

	std::string get_name() const {
		std::vector<matrix_store::const_ptr> mats(2);
		mats[0] = left_mat;
		if (right_mat)
			mats[1] = right_mat;
		else
			mats[1] = left_mat;
		return portion_op->to_string(mats);
	}

	virtual std::vector<safs::io_interface::ptr> create_ios() const;
	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const;
};

static inline bool is_in_mem(matrix_store::const_ptr left_mat,
		matrix_store::const_ptr right_mat)
{
	if (right_mat)
		return left_mat->is_in_mem() && right_mat->is_in_mem();
	else
		return left_mat->is_in_mem();
}

IPW_compute_store::IPW_compute_store(matrix_store::const_ptr left_mat,
		matrix_store::const_ptr right_mat, bulk_operate::const_ptr left_op,
		bulk_operate::const_ptr right_op,
		portion_mapply_op::const_ptr portion_op,
		matrix_layout_t layout): sink_compute_store(left_mat->get_num_rows(),
			left_mat->get_num_cols(), fm::detail::is_in_mem(left_mat, right_mat),
			left_mat->get_type())
{
	this->left_mat = left_mat;
	this->right_mat = right_mat;
	// We always transpose the left matrix, so IPW compute matrix is always tall.
	if (right_mat)
		assert(left_mat->get_num_rows() == right_mat->get_num_rows());
	this->left_op = left_op;
	this->right_op = right_op;
	this->portion_op = std::dynamic_pointer_cast<const combine_op>(portion_op);
	assert(this->portion_op);
	this->layout = layout;
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

local_matrix_store::const_ptr IPW_compute_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	assert(start_col == 0);
	assert(num_cols == left_mat->get_num_cols());
	local_matrix_store::const_ptr left_part = left_mat->get_portion(start_row,
			start_col, num_rows, num_cols);
	local_matrix_store::const_ptr right_part;
	if (right_mat)
		right_part = right_mat->get_portion(start_row, 0, num_rows,
				right_mat->get_num_cols());
	assert(left_part->get_num_rows() == right_part->get_num_rows());
	return create_lmaterialize_matrix(left_part, right_part, get_type(),
			portion_op);
}

local_matrix_store::const_ptr IPW_compute_store::get_portion(size_t id) const
{
	local_matrix_store::const_ptr left_part = left_mat->get_portion(id);
	local_matrix_store::const_ptr right_part;
	if (right_mat) {
		right_part = right_mat->get_portion(id);
		assert(left_part->get_num_rows() == right_part->get_num_rows());
	}
	return create_lmaterialize_matrix(left_part, right_part, get_type(),
			portion_op);
}

async_cres_t IPW_compute_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, std::shared_ptr<portion_compute> compute) const
{
	assert(start_col == 0);
	assert(num_cols == left_mat->get_num_cols());
	collect_portion_compute::ptr new_compute(new collect_portion_compute(
				compute));
	async_cres_t left_ret = left_mat->get_portion_async(start_row, start_col,
			num_rows, num_cols, new_compute);
	async_cres_t right_ret;
	if (right_mat) {
		right_ret = right_mat->get_portion_async(start_row, 0,
				num_rows, right_mat->get_num_cols(), new_compute);
		assert(left_ret.second->get_num_rows()
				== right_ret.second->get_num_rows());
	}
	else
		right_ret.first = true;
	if (left_ret.first && right_ret.first)
		return async_cres_t(true, create_lmaterialize_matrix(left_ret.second,
					right_ret.second, get_type(), portion_op));
	else {
		new_compute->set_EM_count(!left_ret.first + !right_ret.first);
		return async_cres_t(false, create_lmaterialize_matrix(left_ret.second,
					right_ret.second, get_type(), portion_op));
	}
}

std::vector<safs::io_interface::ptr> IPW_compute_store::create_ios() const
{
	std::vector<safs::io_interface::ptr> ret;
	if (!left_mat->is_in_mem()) {
		const EM_object *obj = dynamic_cast<const EM_object *>(left_mat.get());
		std::vector<safs::io_interface::ptr> tmp = obj->create_ios();
		ret.insert(ret.end(), tmp.begin(), tmp.end());
	}
	if (right_mat && !right_mat->is_in_mem()) {
		const EM_object *obj = dynamic_cast<const EM_object *>(right_mat.get());
		std::vector<safs::io_interface::ptr> tmp = obj->create_ios();
		ret.insert(ret.end(), tmp.begin(), tmp.end());
	}
	return ret;
}

std::unordered_map<size_t, size_t> IPW_compute_store::get_underlying_mats() const
{
	std::unordered_map<size_t, size_t> final_res = left_mat->get_underlying_mats();
	if (right_mat) {
		std::unordered_map<size_t, size_t> right = right_mat->get_underlying_mats();
		for (auto it = right.begin(); it != right.end(); it++) {
			auto to_it = final_res.find(it->first);
			if (to_it == final_res.end())
				final_res.insert(std::pair<size_t, size_t>(it->first, it->second));
		}
	}
	return final_res;
}

std::vector<virtual_matrix_store::const_ptr> IPW_matrix_store::get_compute_matrices() const
{
	// If the IPW matrix has been materialized, we don't need to do
	// anything.
	if (has_materialized())
		return std::vector<virtual_matrix_store::const_ptr>();
	else
		return std::vector<virtual_matrix_store::const_ptr>(1,
				virtual_matrix_store::const_ptr(new IPW_compute_store(
						// We always transpose the left matrix.
						left_mat->transpose(), right_mat, left_op, right_op,
						portion_op, layout)));
}

}

}
