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
	virtual const std::vector<detail::local_matrix_store::ptr> &get_partial_results() const = 0;
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

	const std::vector<detail::local_matrix_store::ptr> &get_partial_results() const {
		return local_ms;
	}

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

class multiply_wide_op: public combine_op
{
	std::vector<detail::local_matrix_store::ptr> Abufs;
	std::vector<detail::local_matrix_store::ptr> Bbufs;
	std::vector<detail::local_matrix_store::ptr> res_bufs;
	bool require_trans;
	size_t out_num_rows;
	size_t out_num_cols;
	matrix_layout_t Alayout;
	matrix_layout_t Blayout;
public:
	multiply_wide_op(size_t num_threads, size_t out_num_rows, size_t out_num_cols,
			matrix_layout_t required_layout, const scalar_type &type): combine_op(
				0, 0, type) {
		Abufs.resize(num_threads);
		Bbufs.resize(num_threads);
		res_bufs.resize(num_threads);
		this->out_num_rows = out_num_rows;
		this->out_num_cols = out_num_cols;
		Alayout = required_layout;
		Blayout = required_layout;
		require_trans = false;
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

	virtual const std::vector<detail::local_matrix_store::ptr> &get_partial_results(
			) const {
		return res_bufs;
	}

	void run_part(
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
			Bmat, Bsize.first, 1, res_mat, out_num_rows);
}

template<>
void wide_gemm_col<float>(const std::pair<size_t, size_t> &Asize,
		const float *Amat, const std::pair<size_t, size_t> &Bsize,
		const float *Bmat, float *res_mat, size_t out_num_rows)
{
	cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			Asize.first, Bsize.second, Asize.second, 1, Amat, Asize.first,
			Bmat, Bsize.first, 1, res_mat, out_num_rows);
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
			Bmat, Bsize.second, 1, res_mat, out_num_cols);
}

template<>
void wide_gemm_row<float>(const std::pair<size_t, size_t> &Asize,
		const float *Amat, const std::pair<size_t, size_t> &Bsize,
		const float *Bmat, float *res_mat, size_t out_num_cols)
{
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			Asize.first, Bsize.second, Asize.second, 1, Amat, Asize.second,
			Bmat, Bsize.second, 1, res_mat, out_num_cols);
}

void multiply_wide_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins) const
{
	assert(ins.size() == 2);
	size_t LONG_DIM_LEN = 1024;
	size_t long_dim = ins[1]->get_num_rows();
	if (long_dim <= LONG_DIM_LEN) {
		run_part(ins);
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
		run_part(ins);
	}
	mutableA.restore_size(orig_A);
	mutableB.restore_size(orig_B);
}

void multiply_wide_op::run_part(
		const std::vector<detail::local_matrix_store::const_ptr> &ins) const
{
	int thread_id = detail::mem_thread_pool::get_curr_thread_id();

	detail::local_matrix_store::const_ptr Astore = ins[0];
	const void *Amat = Astore->get_raw_arr();
	if (Amat == NULL || Astore->store_layout() != Alayout) {
		if (Abufs[thread_id] == NULL
				|| Astore->get_num_rows() != Abufs[thread_id]->get_num_rows()
				|| Astore->get_num_cols() != Abufs[thread_id]->get_num_cols()) {
			if (Alayout == matrix_layout_t::L_ROW)
				const_cast<multiply_wide_op *>(this)->Abufs[thread_id]
					= detail::local_matrix_store::ptr(
							new fm::detail::local_buf_row_matrix_store(0, 0,
								Astore->get_num_rows(), Astore->get_num_cols(),
								Astore->get_type(), -1));
			else
				const_cast<multiply_wide_op *>(this)->Abufs[thread_id]
					= detail::local_matrix_store::ptr(
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
				const_cast<multiply_wide_op *>(this)->Bbufs[thread_id]
					= detail::local_matrix_store::ptr(
							new fm::detail::local_buf_col_matrix_store(0, 0,
								Bstore->get_num_rows(), Bstore->get_num_cols(),
								Bstore->get_type(), -1));
			else
				const_cast<multiply_wide_op *>(this)->Bbufs[thread_id]
					= detail::local_matrix_store::ptr(
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
			const_cast<multiply_wide_op *>(this)->res_bufs[thread_id]
				= detail::local_matrix_store::ptr(
						new fm::detail::local_buf_col_matrix_store(0, 0,
							out_num_rows, out_num_cols, get_output_type(), -1));
		else
			const_cast<multiply_wide_op *>(this)->res_bufs[thread_id]
				= detail::local_matrix_store::ptr(
						new fm::detail::local_buf_row_matrix_store(0, 0,
							out_num_rows, out_num_cols, get_output_type(), -1));
		res_bufs[thread_id]->reset_data();
	}
	assert(res_bufs[thread_id]->store_layout() == Blayout);
	void *res_mat = res_bufs[thread_id]->get_raw_arr();
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
		double *t_res_mat = reinterpret_cast<double *>(res_mat);
		if (Blayout == matrix_layout_t::L_COL)
			wide_gemm_col<double>(Asize, t_Amat, Bsize, t_Bmat, t_res_mat,
					out_num_rows);
		else
			wide_gemm_row<double>(Asize, t_Amat, Bsize, t_Bmat, t_res_mat,
					out_num_cols);
	}
	else {
		const float *t_Amat = reinterpret_cast<const float *>(Amat);
		const float *t_Bmat = reinterpret_cast<const float *>(Bmat);
		float *t_res_mat = reinterpret_cast<float *>(res_mat);
		if (Blayout == matrix_layout_t::L_COL)
			wide_gemm_col<float>(Asize, t_Amat, Bsize, t_Bmat, t_res_mat,
					out_num_rows);
		else
			wide_gemm_row<float>(Asize, t_Amat, Bsize, t_Bmat, t_res_mat,
					out_num_cols);
	}
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
		this->left_op = NULL;
		assert(left->get_type() == get_scalar_type<double>()
				|| left->get_type() == get_scalar_type<float>());
		assert(right->get_type() == get_scalar_type<double>()
				|| right->get_type() == get_scalar_type<float>());
		this->right_op = bulk_operate::conv2ptr(get_type().get_basic_ops().get_add());

		matrix_layout_t required_layout;
		// If both input matrices have the same data layout, it's easy.
		if (left->store_layout() == right->store_layout())
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
					right->get_num_cols(), required_layout, left->get_type()));
	}
	else {
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
	std::vector<detail::local_matrix_store::ptr> local_ms
		= std::static_pointer_cast<const inner_prod_wide_op>(
				portion_op)->get_partial_results();
	matrix_layout_t local_layout = matrix_layout_t::L_NONE;
	for (size_t i = 0; i < local_ms.size(); i++) {
		if (local_ms[i]) {
			local_layout = local_ms[i]->store_layout();
			break;
		}
	}
	assert(local_layout != matrix_layout_t::L_NONE);

	// Aggregate the results from omp threads.
	detail::matrix_store::ptr res = detail::matrix_store::create(
			left_mat->get_num_rows(), right_mat->get_num_cols(), local_layout,
			right_op->get_output_type(), -1, true);
	res->reset_data();
	detail::local_matrix_store::ptr local_res = res->get_portion(0);
	assert(local_res->get_num_rows() == res->get_num_rows()
			&& local_res->get_num_cols() == res->get_num_cols());
	for (size_t j = 0; j < local_ms.size(); j++) {
		// It's possible that the local matrix store doesn't exist
		// because the input matrix is very small.
		if (local_ms[j])
			detail::mapply2(*local_res, *local_ms[j], *right_op, *local_res);
	}
	if (this->layout == local_layout)
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
		new_compute->set_EM_count(left_ret.first + right_ret.first);
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
