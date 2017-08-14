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

#include <cblas.h>

#include "log.h"

#include "dense_matrix.h"
#include "bulk_operate.h"
#include "NUMA_dense_matrix.h"
#include "EM_dense_matrix.h"
#include "generic_type.h"
#include "one_val_matrix_store.h"
#include "local_matrix_store.h"
#include "virtual_matrix_store.h"
#include "mapply_matrix_store.h"
#include "vector.h"
#include "matrix_stats.h"
#include "local_mem_buffer.h"
#include "factor.h"
#include "cached_matrix_store.h"
#include "agg_matrix_store.h"
#include "IPW_matrix_store.h"
#include "groupby_matrix_store.h"
#include "block_matrix.h"
#include "col_vec.h"
#include "sink_matrix.h"
#include "data_frame.h"
#include "project_matrix_store.h"
#include "set_data_matrix_store.h"
#include "set_rc_matrix_store.h"
#include "generic_hashtable.h"
#include "local_vec_store.h"
#include "cum_matrix.h"

namespace fm
{

namespace detail
{

void portion_mapply_op::run(
		const std::vector<std::shared_ptr<const local_matrix_store> > &ins) const
{
	BOOST_LOG_TRIVIAL(error)
		<< "It doesn't support running on only input matrices";
	assert(0);
}

void portion_mapply_op::run(
		const std::vector<std::shared_ptr<const local_matrix_store> > &ins,
		local_matrix_store &out) const
{
	BOOST_LOG_TRIVIAL(error)
		<< "It doesn't support running on input matrices and output one matrix";
	assert(0);
}

void portion_mapply_op::run(
		const std::vector<std::shared_ptr<const local_matrix_store> > &ins,
		const std::vector<std::shared_ptr<local_matrix_store> > &outs) const
{
	BOOST_LOG_TRIVIAL(error)
		<< "It doesn't support running on input matrices and output multiple matrices";
	assert(0);
}

}

vector::ptr dense_matrix::conv2vec() const
{
	materialize_self();
	auto ret = store->conv2vec();
	if (ret)
		return vector::create(ret);
	else
		return vector::ptr();
}

bool dense_matrix::verify_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	if (this->get_entry_size() != left_op.left_entry_size()
			|| m.get_entry_size() != left_op.right_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The left operator isn't compatible with input matrices";
		return false;
	}

	if (left_op.output_entry_size() != right_op.left_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The type of the left operator doesn't match the right operator";
		return false;
	}

	if (right_op.left_entry_size() != right_op.right_entry_size()
			|| right_op.left_entry_size() != right_op.output_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The input and output of the right operator has different types";
		return false;
	}

	if (get_num_cols() != m.get_num_rows()) {
		BOOST_LOG_TRIVIAL(error) << "The matrix size doesn't match";
		return false;
	}
	return true;
}

bool dense_matrix::verify_mapply2(const dense_matrix &m,
			const bulk_operate &op) const
{
	if (this->get_num_rows() != m.get_num_rows()
			|| this->get_num_cols() != m.get_num_cols()) {
		BOOST_LOG_TRIVIAL(error)
			<< "two matrices in mapply2 don't have the same shape";
		return false;
	}

	if (get_entry_size() != op.left_entry_size()
			|| m.get_entry_size() != op.right_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "the element type in the matrices isn't compatible with the operator";
		return false;
	}

	return true;
}

bool dense_matrix::verify_apply(matrix_margin margin, const arr_apply_operate &op) const
{
	if (get_entry_size() != op.input_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "the element type in the matrices isn't compatible with the operator";
		return false;
	}

	return true;
}

namespace
{

class double_square: public bulk_uoperate
{
public:
	virtual void runA(size_t num_eles, const void *in_arr,
			void *out_arr) const {
		long double *t_out_arr = (long double *) out_arr;
		const double *t_in_arr = (const double *) in_arr;
		for (size_t i = 0; i < num_eles; i++)
			t_out_arr[i]
				= ((long double) t_in_arr[i]) * ((long double) t_in_arr[i]);
	}
	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<double>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<long double>();
	}
	virtual std::string get_name() const {
		return "dsqrt";
	}
};

class sum_agg: public bulk_operate
{
public:
	virtual void runAgg(size_t num_eles, const void *left_arr1,
			void *output) const {
		const long double *t_input = (const long double *) left_arr1;
		long double *t_output = (long double *) output;
		t_output[0] = t_input[0];
		for (size_t i = 1; i < num_eles; i++)
			t_output[0] += t_input[i];
	}

	virtual void runAA(size_t num_eles, const void *left_arr,
			const void *right_arr, void *output_arr) const {
		assert(0);
	}

	virtual void runAE(size_t num_eles, const void *left_arr,
			const void *right, void *output_arr) const {
		assert(0);
	}

	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr, void *output_arr) const {
		assert(0);
	}
	virtual void runCum(size_t num_eles, const void *left_arr,
			const void *prev, void *output) const {
		assert(0);
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<long double>();
	}

	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<long double>();
	}

	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<long double>();
	}
	virtual std::string get_name() const {
		return "ldsum";
	}
};

}

double dense_matrix::norm2() const
{
	detail::matrix_stats.inc_multiplies(get_num_rows() * get_num_cols());
	double ret = 0;
	if (get_type() == get_scalar_type<double>()) {
		dense_matrix::ptr sq_mat
			= this->sapply(bulk_uoperate::const_ptr(new double_square()));
		assert(sq_mat->get_type() == get_scalar_type<long double>());
		scalar_variable::ptr res = sq_mat->aggregate(
				bulk_operate::const_ptr(new sum_agg()));
		assert(res->get_type() == get_scalar_type<long double>());
		ret = sqrtl(*(long double *) res->get_raw());
	}
	else {
		const bulk_uoperate *op = get_type().get_basic_uops().get_op(
				basic_uops::op_idx::SQ);
		dense_matrix::ptr sq_mat = this->sapply(bulk_uoperate::conv2ptr(*op));
		scalar_variable::ptr res = sq_mat->aggregate(bulk_operate::conv2ptr(
					sq_mat->get_type().get_basic_ops().get_add()));
		res->get_type().get_basic_uops().get_op(
				basic_uops::op_idx::SQRT)->runA(1, res->get_raw(), &ret);
	}
	return ret;
}

namespace
{

class multiply_tall_op: public detail::portion_mapply_op
{
	detail::mem_matrix_store::const_ptr Bmem;
	detail::local_matrix_store::const_ptr Bstore;
	std::vector<detail::local_matrix_store::ptr> Abufs;
	std::vector<detail::local_matrix_store::ptr> res_bufs;
public:
	multiply_tall_op(detail::mem_matrix_store::const_ptr Bmem,
			size_t num_threads, size_t out_num_rows,
			size_t out_num_cols): detail::portion_mapply_op(
				out_num_rows, out_num_cols, Bmem->get_type()) {
		// We need B matrix stores data in contiguous memory.
		if (Bmem->get_raw_arr()) {
			this->Bmem = Bmem;
			this->Bstore = Bmem->get_portion(0);
		}
		else {
			detail::mem_matrix_store::ptr copy = detail::mem_matrix_store::create(
					Bmem->get_num_rows(), Bmem->get_num_cols(),
					Bmem->store_layout(), Bmem->get_type(), -1);
			detail::local_matrix_store::ptr copy_portion = copy->get_portion(0);
			detail::local_matrix_store::const_ptr Bstore = Bmem->get_portion(0);
			copy_portion->copy_from(*Bstore);
			this->Bmem = copy;
			this->Bstore = copy_portion;
		}
		assert(Bstore->get_raw_arr());
		assert(this->Bstore->get_num_rows() == Bmem->get_num_rows());
		assert(this->Bstore->get_num_cols() == Bmem->get_num_cols());

		Abufs.resize(num_threads);
		res_bufs.resize(num_threads);
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		// We can only resize the number of the rows in this operation.
		return local_num_cols == get_out_num_cols();
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	virtual detail::portion_mapply_op::const_ptr transpose() const;
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("(") + (mats[0]->get_name()
					+ "%*%") + Bmem->get_name() + std::string(")");
	}
};

class t_multiply_tall_op: public detail::portion_mapply_op
{
	multiply_tall_op op;
public:
	t_multiply_tall_op(const multiply_tall_op &_op): detail::portion_mapply_op(
			_op.get_out_num_cols(), _op.get_out_num_rows(),
			_op.get_output_type()), op(_op) {
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		assert(ins.size() == 1);
		std::vector<fm::detail::local_matrix_store::const_ptr> t_ins(ins.size());
		t_ins[0] = std::static_pointer_cast<const fm::detail::local_matrix_store>(
				ins[0]->transpose());
		fm::detail::local_matrix_store::ptr t_out
			= std::static_pointer_cast<fm::detail::local_matrix_store>(
					out.transpose());
		op.run(t_ins, *t_out);
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		return op.is_resizable(local_start_col, local_start_row,
				local_num_cols, local_num_rows);
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		return fm::detail::portion_mapply_op::const_ptr(
				new multiply_tall_op(op));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return op.to_string(mats);
	}
};

detail::portion_mapply_op::const_ptr multiply_tall_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new t_multiply_tall_op(*this));
}

void multiply_tall_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	detail::local_matrix_store::const_ptr Astore = ins[0];
	detail::matrix_stats.inc_multiplies(
			Astore->get_num_rows() * Astore->get_num_cols() * Bstore->get_num_cols());

	int thread_id = detail::mem_thread_pool::get_curr_thread_id();
	// TODO Maybe I should store the pairs in the vector directly.
	std::pair<detail::local_matrix_store::ptr, detail::local_matrix_store::ptr> bufs(
			Abufs[thread_id], res_bufs[thread_id]);
	detail::matrix_tall_multiply(*Astore, *Bstore, out, bufs);

	multiply_tall_op *mutable_this = const_cast<multiply_tall_op *>(this);
	if (bufs.first != Abufs[thread_id])
		mutable_this->Abufs[thread_id] = bufs.first;
	if (bufs.second != res_bufs[thread_id])
		mutable_this->res_bufs[thread_id] = bufs.second;
}

dense_matrix::ptr blas_multiply_tall(const dense_matrix &m1,
		const dense_matrix &m2, matrix_layout_t out_layout)
{
	if (m2.get_num_cols() > detail::mem_matrix_store::CHUNK_SIZE) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't multiply a tall matrix with a wide matrix";
		return dense_matrix::ptr();
	}

	if (out_layout == matrix_layout_t::L_NONE)
		out_layout = m1.store_layout();

	assert(m1.get_type() == m2.get_type());
	detail::matrix_store::const_ptr right = m2.get_raw_store();
	if (out_layout != m2.store_layout()) {
		dense_matrix::ptr tmp = m2.conv2(out_layout);
		right = tmp->get_raw_store();
	}
	// TODO the right matrix might be a sparse matrix. conv_store
	// will turn it into a dense matrix.
	// We should optimize for the right sparse matrix later.
	if (right->is_virtual() || !right->is_in_mem() || right->get_num_nodes() > 0) {
		dense_matrix::ptr tmp = dense_matrix::create(right);
		tmp = tmp->conv_store(true, -1);
		right = tmp->get_raw_store();
	}
	assert(right->store_layout() == out_layout);
	assert(!right->is_virtual());
	assert(right->get_num_nodes() == -1);
	assert(right->is_in_mem());

	std::vector<detail::matrix_store::const_ptr> ins(1);
	ins[0] = m1.get_raw_store();
	multiply_tall_op::const_ptr mapply_op(new multiply_tall_op(
				detail::mem_matrix_store::cast(right),
				detail::mem_thread_pool::get_global_num_threads(),
				m1.get_num_rows(), m2.get_num_cols()));
	return dense_matrix::create(__mapply_portion_virtual(ins, mapply_op,
				out_layout));
}

dense_matrix::ptr blas_multiply_wide(const dense_matrix &m1,
		const dense_matrix &m2, matrix_layout_t out_layout)
{
	detail::matrix_stats.inc_multiplies(
			m1.get_num_rows() * m1.get_num_cols() * m2.get_num_cols());

	assert(m1.get_type() == m2.get_type());

	detail::matrix_store::ptr res = detail::IPW_matrix_store::create(
				m1.get_raw_store(), m2.get_raw_store(), NULL, NULL, out_layout);
	return dense_matrix::create(res);
}

}

dense_matrix::ptr dense_matrix::multiply_sparse_combined(
		const dense_matrix &mat, matrix_layout_t out_layout) const
{
	detail::combined_matrix_store::const_ptr store
		= std::dynamic_pointer_cast<const detail::combined_matrix_store>(
				mat.get_raw_store());
	assert(store);

	std::vector<dense_matrix::ptr> res_mats(store->get_num_mats());
	for (size_t i = 0; i < store->get_num_mats(); i++) {
		dense_matrix::ptr right = dense_matrix::create(store->get_mat(i));
		res_mats[i] = multiply(*right, out_layout);
	}
	materialize(res_mats, false);
	return dense_matrix::cbind(res_mats);
}

static inline bool is_floating_point(const scalar_type &type)
{
	return type == get_scalar_type<float>() || type == get_scalar_type<double>();
}

dense_matrix::ptr dense_matrix::multiply(const dense_matrix &mat,
		matrix_layout_t out_layout) const
{
	if (get_num_cols() != mat.get_num_rows()) {
		BOOST_LOG_TRIVIAL(error)
			<< "Matrix multiply: #cols of mat1 should be the same as #rows of mat2";
		return dense_matrix::ptr();
	}

	// If input matrices are sink matrices, we materialize them first.
	if (store->is_sink())
		materialize_self();
	if (mat.store->is_sink())
		mat.materialize_self();

	// We should treat a sparse matrix differently to improve performance of
	// matrix multiplication.
	// TODO right now this optimization is only useful for a wide matrix
	// times a tall matrix.
	if (is_floating_point(get_type()) && is_floating_point(mat.get_type())
			&& mat.get_data().is_sparse() && this->is_wide() && !mat.is_wide()) {
		detail::sparse_project_matrix_store::const_ptr store
			= std::dynamic_pointer_cast<const detail::sparse_project_matrix_store>(
					mat.get_raw_store());

		if (store) {
			if (out_layout == matrix_layout_t::L_NONE)
				out_layout = matrix_layout_t::L_COL;
			dense_matrix::ptr tmp = conv2(matrix_layout_t::L_COL);
			// TODO we assume the left matrix is a wide matrix.
			return dense_matrix::create(detail::IPW_matrix_store::create(
						tmp->get_raw_store(), store, NULL, NULL, out_layout));
		}
		else
			return multiply_sparse_combined(mat, out_layout);
	}
	else if (is_floating_point(get_type())) {
		assert(get_type() == mat.get_type());
		size_t long_dim1 = std::max(get_num_rows(), get_num_cols());
		size_t long_dim2 = std::max(mat.get_num_rows(), mat.get_num_cols());
		// We prefer to perform computation on the larger matrix.
		// If the matrix in the right operand is larger, we transpose
		// the entire computation.
		if (long_dim2 > long_dim1) {
			dense_matrix::ptr t_mat1 = this->transpose();
			dense_matrix::ptr t_mat2 = mat.transpose();
			matrix_layout_t t_layout = out_layout;
			if (t_layout == matrix_layout_t::L_ROW)
				t_layout = matrix_layout_t::L_COL;
			else if (t_layout == matrix_layout_t::L_COL)
				t_layout = matrix_layout_t::L_ROW;
			dense_matrix::ptr t_res = t_mat2->multiply(*t_mat1, t_layout);
			if (t_res)
				return t_res->transpose();
			else
				return dense_matrix::ptr();
		}

		if (is_wide())
			return blas_multiply_wide(*this, mat, out_layout);
		else
			return blas_multiply_tall(*this, mat, out_layout);
	}
	else {
		bulk_operate::const_ptr multiply = bulk_operate::conv2ptr(
				get_type().get_basic_ops().get_multiply());
		bulk_operate::const_ptr add = bulk_operate::conv2ptr(
				get_type().get_basic_ops().get_add());
		return inner_prod(mat, multiply, add, out_layout);
	}
}

namespace
{

class apply_scalar_op: public bulk_uoperate
{
	scalar_variable::const_ptr var;
	bulk_operate::const_ptr op;
public:
	apply_scalar_op(scalar_variable::const_ptr var, bulk_operate::const_ptr op) {
		this->var = var;
		this->op = op;
	}
	virtual void runA(size_t num_eles, const void *in_arr, void *out_arr) const {
		op->runAE(num_eles, in_arr, var->get_raw(), out_arr);
	}
	virtual const scalar_type &get_input_type() const {
		return op->get_left_type();
	}
	virtual const scalar_type &get_output_type() const {
		return op->get_output_type();
	}
	virtual std::string get_name() const {
		return "apply_scalar";
	}
};

}

dense_matrix::ptr dense_matrix::apply_scalar(
		scalar_variable::const_ptr var, bulk_operate::const_ptr op) const
{
	if (get_type() != var->get_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't multiply a scalar of incompatible type";
		return dense_matrix::ptr();
	}
	return sapply(bulk_uoperate::const_ptr(new apply_scalar_op(var, op)));
}

bool dense_matrix::materialize_self() const
{
	if (!store->is_virtual())
		return true;

	detail::matrix_store::const_ptr tmp;
	try {
		auto vstore = detail::virtual_matrix_store::cast(store);
		detail::sink_store::materialize_matrices(vstore);
		tmp = vstore->materialize(store->is_in_mem(), store->get_num_nodes());
	} catch (std::exception &e) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("fail to materialize: %1%") % e.what();
	}
	if (tmp == NULL)
		return false;

	const_cast<dense_matrix *>(this)->store = tmp;
	return true;
}

std::vector<detail::virtual_matrix_store::const_ptr> dense_matrix::get_compute_matrices() const
{
	// If the dense matrix isn't virtual.
	if (!is_virtual())
		return std::vector<detail::virtual_matrix_store::const_ptr>();

	detail::sink_store::const_ptr sink
		= std::dynamic_pointer_cast<const detail::sink_store>(get_raw_store());
	// If the dense matrix isn't a sink matrix.
	if (sink == NULL) {
		detail::virtual_matrix_store::const_ptr vmat
			= std::dynamic_pointer_cast<const detail::virtual_matrix_store>(
					get_raw_store());
		assert(vmat);
		return std::vector<detail::virtual_matrix_store::const_ptr>(1, vmat);
	}
	else
		return sink->get_compute_matrices();
}

void dense_matrix::set_materialize_level(materialize_level level,
		detail::matrix_store::ptr materialize_buf)
{
	const detail::virtual_matrix_store *tmp
		= dynamic_cast<const detail::virtual_matrix_store *>(store.get());
	// If the matrix isn't a virtual matrix, we don't need to materialize it.
	if (tmp == NULL)
		return;

	detail::virtual_matrix_store *tmp1
		= const_cast<detail::virtual_matrix_store *>(tmp);
	tmp1->set_materialize_level(level, materialize_buf);
}

///////////////////////// Scale rows and columns //////////////////////////////

namespace
{

class mapply_col_op: public detail::portion_mapply_op
{
	detail::mem_matrix_store::const_ptr vals;
	bulk_operate::const_ptr op;
public:
	mapply_col_op(detail::mem_matrix_store::const_ptr vals, bulk_operate::const_ptr op,
			size_t out_num_rows, size_t out_num_cols): detail::portion_mapply_op(
				out_num_rows, out_num_cols, op->get_output_type()) {
		this->vals = vals;
		this->op = op;
		if (vals)
			assert(vals->get_num_cols() == 1);
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		// We can only resize the number of the cols in this operation.
		return local_num_rows == get_out_num_rows();
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;
	virtual portion_mapply_op::const_ptr transpose() const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() >= 1);
		std::string ret = std::string("mapply_col(") + mats[0]->get_name();
		if (mats.size() == 1)
			ret += ", vec";
		else
			ret += mats[1]->get_name();
		return ret + std::string(")");
	}
};

class mapply_row_op: public detail::portion_mapply_op
{
	detail::mem_matrix_store::const_ptr vals;
	bulk_operate::const_ptr op;
public:
	mapply_row_op(detail::mem_matrix_store::const_ptr vals, bulk_operate::const_ptr op,
			size_t out_num_rows, size_t out_num_cols): detail::portion_mapply_op(
				out_num_rows, out_num_cols, op->get_output_type()) {
		this->vals = vals;
		this->op = op;
		if (vals)
			assert(vals->get_num_cols() == 1);
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		// We can only resize the number of the rows in this operation.
		return local_num_cols == get_out_num_cols();
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;
	virtual portion_mapply_op::const_ptr transpose() const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() >= 1);
		std::string ret = std::string("mapply_row(") + mats[0]->get_name();
		if (mats.size() == 1)
			ret += ", vec";
		else
			ret += mats[1]->get_name();
		return ret + std::string(")");
	}
};

detail::portion_mapply_op::const_ptr mapply_col_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new mapply_row_op(vals, op,
				get_out_num_cols(), get_out_num_rows()));
}

void mapply_row_op::run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() >= 1);
	detail::matrix_stats.inc_multiplies(
			ins[0]->get_num_rows() * ins[0]->get_num_cols());

	assert(ins[0]->get_global_start_col() == out.get_global_start_col());
	assert(ins[0]->get_global_start_row() == out.get_global_start_row());
	// This is a tall matrix. We divide the matrix horizontally.
	if (ins.size() == 1) {
		assert(vals);
		assert(ins[0]->get_global_start_col() == 0
				&& vals->get_num_rows() == ins[0]->get_num_cols());
		// vals is a column vector.
		const char *col = vals->get_col(0);
		assert(col);
		local_cref_vec_store lvals(col, 0, vals->get_num_rows(),
				vals->get_type(), -1);
		detail::mapply_rows(*ins[0], lvals, *op, out);
	}
	// We divide the matrix vertically.
	else {
		assert(ins.size() == 2);
		assert(ins[1]->get_num_rows() == 1);
		assert(ins[1]->store_layout() == matrix_layout_t::L_ROW);
		detail::local_row_matrix_store::const_ptr row_in1
			= std::static_pointer_cast<const detail::local_row_matrix_store>(
					ins[1]);
		const char *arr = row_in1->get_row(0);
		assert(arr);
		local_cref_vec_store lvals(arr, ins[1]->get_global_start_col(),
				ins[1]->get_num_cols(), ins[1]->get_type(), -1);
		detail::mapply_rows(*ins[0], lvals, *op, out);
	}
}

void mapply_col_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() >= 1);
	detail::matrix_stats.inc_multiplies(
			ins[0]->get_num_rows() * ins[0]->get_num_cols());

	assert(ins[0]->get_global_start_col() == out.get_global_start_col());
	assert(ins[0]->get_global_start_row() == out.get_global_start_row());
	// This is a wide matrix. We divide the matrix vertically.
	if (ins.size() == 1) {
		assert(vals);
		assert(ins[0]->get_global_start_row() == 0
				&& vals->get_num_rows() == ins[0]->get_num_rows());
		// If we use get_raw_arr, it may not work with NUMA vector.
		const char *col = vals->get_col(0);
		assert(col);
		local_cref_vec_store lvals(col, 0, vals->get_num_rows(),
				vals->get_type(), -1);
		detail::mapply_cols(*ins[0], lvals, *op, out);
	}
	// We divide the tall matrix horizontally.
	else {
		assert(ins.size() == 2);
		assert(ins[1]->get_num_cols() == 1);
		assert(ins[1]->store_layout() == matrix_layout_t::L_COL);
		detail::local_col_matrix_store::const_ptr col_in1
			= std::static_pointer_cast<const detail::local_col_matrix_store>(
					ins[1]);
		const char *arr = col_in1->get_col(0);
		assert(arr);
		local_cref_vec_store lvals(arr, ins[1]->get_global_start_row(),
				ins[1]->get_num_rows(), ins[1]->get_type(), -1);
		detail::mapply_cols(*ins[0], lvals, *op, out);
	}
}

detail::portion_mapply_op::const_ptr mapply_row_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new mapply_col_op(vals, op,
				get_out_num_cols(), get_out_num_rows()));
}

}

dense_matrix::ptr dense_matrix::mapply_cols(col_vec::const_ptr vals,
		bulk_operate::const_ptr op) const
{
	if (get_num_rows() != vals->get_length()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The vector's length needs to equal to #rows";
		return dense_matrix::ptr();
	}
	if (get_type() != op->get_left_type()
			|| vals->get_type() != op->get_right_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "mapply_rows: the input type is different from what the bulk operator expects";
		return dense_matrix::ptr();
	}

	// If input matrices are sink matrices, we materialize them first.
	if (store->is_sink())
		materialize_self();
	if (vals->store->is_sink())
		vals->materialize_self();

	std::vector<detail::matrix_store::const_ptr> ins(1);
	ins[0] = this->get_raw_store();
	mapply_col_op::const_ptr mapply_op;
	// If this is a wide matrix or a square matrix.
	if (is_wide() || get_num_rows() == get_num_cols()) {
		dense_matrix::ptr mem_mat = vals->conv_store(true, -1);
		mapply_op = mapply_col_op::const_ptr(new mapply_col_op(
					detail::mem_matrix_store::cast(mem_mat->get_raw_store()),
					op, get_num_rows(), get_num_cols()));
	}
	// If this is a tall matrix, the input vector may also be stored on
	// disks. We should give it as the input of mapply_portion.
	else {
		ins.push_back(vals->get_raw_store());
		mapply_op = mapply_col_op::const_ptr(new mapply_col_op(
					NULL, op, get_num_rows(), get_num_cols()));
	}
	detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
			mapply_op, this->store_layout());
	return dense_matrix::create(ret);
}

dense_matrix::ptr dense_matrix::mapply_rows(col_vec::const_ptr vals,
		bulk_operate::const_ptr op) const
{
	if (get_num_cols() != vals->get_length()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The vector's length needs to equal to #columns";
		return dense_matrix::ptr();
	}
	if (get_type() != op->get_left_type()
			|| vals->get_type() != op->get_right_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "mapply_rows: the input type is different from what the bulk operator expects";
		return dense_matrix::ptr();
	}

	// If input matrices are sink matrices, we materialize them first.
	if (store->is_sink())
		materialize_self();
	if (vals->store->is_sink())
		vals->materialize_self();

	// We should handle mapply_row and mapply_col separately, so we don't
	// mess up the case of square matrices.
	std::vector<detail::matrix_store::const_ptr> ins(1);
	ins[0] = this->get_raw_store();
	mapply_row_op::const_ptr mapply_op;
	// If this is a wide matrix, the input vector may also be stored on
	// disks. We should give it as the input of mapply_portion.
	if (is_wide()) {
		ins.push_back(vals->get_raw_store()->transpose());
		mapply_op = mapply_row_op::const_ptr(new mapply_row_op(
					NULL, op, get_num_rows(), get_num_cols()));
	}
	// If this is a tall matrix or a square matrix.
	else {
		dense_matrix::ptr mem_mat = vals->conv_store(true, -1);
		mapply_op = mapply_row_op::const_ptr(new mapply_row_op(
					detail::mem_matrix_store::cast(mem_mat->get_raw_store()),
					op, get_num_rows(), get_num_cols()));
	}
	detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
			mapply_op, this->store_layout());
	return dense_matrix::create(ret);
}

//////////////////////////// Cast the element types ///////////////////////////

dense_matrix::ptr dense_matrix::cast_ele_type(const scalar_type &type,
		bool forced) const
{
	// If they have the same type, we just return the matrix itself
	// regardless of `forced'.
	if (get_type() == type)
		return clone();
	if (!require_cast(get_type(), type) && !forced)
		return clone();
	else
		return sapply(bulk_uoperate::conv2ptr(get_type().get_type_cast(type)));
}

///////////////////////////////////// mapply2 /////////////////////////////////

namespace
{

class mapply2_op: public detail::portion_mapply_op
{
	bulk_operate::const_ptr op;
public:
	mapply2_op(bulk_operate::const_ptr op, size_t out_num_rows,
			size_t out_num_cols): detail::portion_mapply_op(out_num_rows,
				out_num_cols, op->get_output_type()) {
		this->op = op;
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;
	virtual portion_mapply_op::const_ptr transpose() const {
		return portion_mapply_op::const_ptr(new mapply2_op(op,
					get_out_num_cols(), get_out_num_rows()));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 2);
		return op->get_name() + std::string("(") + mats[0]->get_name()
			+ ", " + mats[1]->get_name() + ")";
	}
};

void mapply2_op::run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 2);
	assert(ins[0]->get_global_start_col() == ins[1]->get_global_start_col());
	assert(ins[0]->get_global_start_col() == out.get_global_start_col());
	assert(ins[0]->get_global_start_row() == ins[1]->get_global_start_row());
	assert(ins[0]->get_global_start_row() == out.get_global_start_row());
	detail::part_dim_t dim = is_wide()
		? detail::part_dim_t::PART_DIM2 : detail::part_dim_t::PART_DIM1;
	detail::mapply2(*ins[0], *ins[1], *op, dim, out);
}

}

dense_matrix::ptr dense_matrix::mapply2(const dense_matrix &m,
		bulk_operate::const_ptr op) const
{
	// The same shape and the same data layout.
	if (!verify_mapply2(m, *op))
		return dense_matrix::ptr();

	// If input matrices are sink matrices, we materialize them first.
	if (store->is_sink() || m.store->is_sink()) {
		std::vector<detail::matrix_store::const_ptr> stores(2);
		stores[0] = get_raw_store();
		stores[1] = m.get_raw_store();
		detail::portion_mapply_op::const_ptr portion_op(new mapply2_op(op,
					get_num_rows(), get_num_cols()));
		return dense_matrix::create(detail::mapply_sink_store::create(
					stores, portion_op));
	}

	std::vector<detail::matrix_store::const_ptr> ins(2);
	ins[0] = this->get_raw_store();
	if (this->store_layout() == m.store_layout())
		ins[1] = m.get_raw_store();
	else {
		dense_matrix::ptr m1 = m.conv2(this->store_layout());
		ins[1] = m1->get_raw_store();
	}
	mapply2_op::const_ptr mapply_op(new mapply2_op(op, get_num_rows(),
				get_num_cols()));
	return dense_matrix::create(__mapply_portion_virtual(ins, mapply_op,
				this->store_layout()));
}

static std::pair<dense_matrix::ptr, dense_matrix::ptr> match_type(
		const dense_matrix &m1, const dense_matrix &m2)
{
	dense_matrix::ptr new_m1 = m1.clone();
	dense_matrix::ptr new_m2 = m2.clone();
	// If these two matrices don't have the same element type, we need to
	// cast one of them to match the other.
	if (m1.get_type() != m2.get_type()) {
		// The element type is ordered based on the type size.
		if ((int) m1.get_type().get_type() > (int) m2.get_type().get_type())
			new_m2 = new_m2->cast_ele_type(new_m1->get_type());
		else
			new_m1 = new_m1->cast_ele_type(new_m2->get_type());
	}
	return std::pair<dense_matrix::ptr, dense_matrix::ptr>(new_m1, new_m2);
}

dense_matrix::ptr dense_matrix::mapply2(const dense_matrix &m,
		basic_ops::op_idx idx) const
{
	if (this->get_type() == m.get_type()) {
		const bulk_operate &op = *get_type().get_basic_ops().get_op(idx);
		return this->mapply2(m, bulk_operate::conv2ptr(op));
	}
	else {
		auto match_res = match_type(*this, m);
		dense_matrix::ptr new_this = match_res.first;
		dense_matrix::ptr new_m = match_res.second;
		const bulk_operate &op = *new_this->get_type().get_basic_ops().get_op(idx);
		return new_this->mapply2(*new_m, bulk_operate::conv2ptr(op));
	}
}

namespace
{

class sapply_op: public detail::portion_mapply_op
{
	bulk_uoperate::const_ptr op;
public:
	sapply_op(bulk_uoperate::const_ptr op, size_t out_num_rows,
			size_t out_num_cols): detail::portion_mapply_op(out_num_rows,
				out_num_cols, op->get_output_type()) {
		this->op = op;
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;
	virtual portion_mapply_op::const_ptr transpose() const {
		return portion_mapply_op::const_ptr(new sapply_op(op, get_out_num_cols(),
					get_out_num_rows()));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return op->get_name() + std::string("(") + mats[0]->get_name() + ")";
	}
};

void sapply_op::run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 1);
	assert(ins[0]->get_global_start_col() == out.get_global_start_col());
	assert(ins[0]->get_global_start_row() == out.get_global_start_row());
	detail::part_dim_t dim = is_wide()
		? detail::part_dim_t::PART_DIM2 : detail::part_dim_t::PART_DIM1;
	detail::sapply(*ins[0], *op, dim, out);
}

}

dense_matrix::ptr dense_matrix::sapply(bulk_uoperate::const_ptr op) const
{
	if (get_type() != op->get_input_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The matrix type doesn't match with the sapply input type";
		return dense_matrix::ptr();
	}

	// If the input matrix is a sink matrix, we materialize it first.
	if (store->is_sink()) {
		std::vector<detail::matrix_store::const_ptr> stores(1);
		stores[0] = get_raw_store();
		detail::portion_mapply_op::const_ptr portion_op(new sapply_op(op,
					get_num_rows(), get_num_cols()));
		return dense_matrix::create(detail::mapply_sink_store::create(
					stores, portion_op));
	}

	std::vector<detail::matrix_store::const_ptr> ins(1);
	ins[0] = this->get_raw_store();
	sapply_op::const_ptr mapply_op(new sapply_op(op, get_num_rows(),
				get_num_cols()));
	detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
			mapply_op, this->store_layout());
	return dense_matrix::create(ret);
}

dense_matrix::dense_matrix(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, int num_nodes, bool in_mem,
			safs::safs_file_group::ptr group)
{
	store = detail::matrix_store::ptr(new detail::one_val_matrix_store(
				type.create_scalar(), nrow, ncol, layout, num_nodes));
}

dense_matrix::ptr dense_matrix::create(detail::matrix_store::const_ptr store)
{
	if (store->get_num_rows() == 0 || store->get_num_cols() == 0) {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't create a matrix with 0 rows/cols";
		return dense_matrix::ptr();
	}

	return dense_matrix::ptr(new dense_matrix(store));
}

dense_matrix::ptr dense_matrix::create_const(scalar_variable::ptr val,
		size_t nrow, size_t ncol, matrix_layout_t layout, int num_nodes,
		bool in_mem, safs::safs_file_group::ptr group)
{
	if (nrow == 0 || ncol == 0) {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't create a matrix with 0 rows/cols";
		return dense_matrix::ptr();
	}

	size_t long_dim = std::max(nrow, ncol);
	size_t short_dim = std::min(nrow, ncol);
	// We don't want to create a small block matrix.
	if (long_dim < detail::EM_matrix_store::CHUNK_SIZE
			|| short_dim <= matrix_conf.get_block_size()) {
		detail::matrix_store::ptr store(new detail::one_val_matrix_store(
					val, nrow, ncol, layout, num_nodes));
		return dense_matrix::ptr(new dense_matrix(store));
	}
	else
		return block_matrix::create_layout(val, nrow, ncol, layout,
				matrix_conf.get_block_size(), num_nodes, in_mem, group);
}

namespace
{

/*
 * In this case, we copy the entire vector to the destination array each time.
 */
class repeat_operate: public set_operate
{
	detail::mem_col_matrix_store::const_ptr vec;
public:
	repeat_operate(detail::mem_col_matrix_store::const_ptr vec) {
		this->vec = vec;
	}

	virtual void set(void *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		assert(num_eles == vec->get_num_rows());
		memcpy(arr, vec->get_raw_arr(), num_eles * vec->get_entry_size());
	}

	virtual const scalar_type &get_type() const {
		return vec->get_type();
	}
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr(new repeat_operate(vec));
	}
};

/*
 * In this case, we repeat one of the elements in the vector in the destination
 * array each time.
 */
class repeat_ele_operate: public set_operate
{
	detail::mem_col_matrix_store::const_ptr vec;
	// Whether or not to repeat an element in a row.
	bool in_row;

	void repeat_ele(void *dst, const void *src, size_t num_eles) const {
		for (size_t i = 0; i < num_eles; i++)
			memcpy(((char *) dst) + i * vec->get_entry_size(), src,
					vec->get_entry_size());
	}
public:
	repeat_ele_operate(detail::mem_col_matrix_store::const_ptr vec, bool in_row) {
		this->vec = vec;
		this->in_row = in_row;
	}

	virtual void set(void *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		if (in_row)
			repeat_ele(arr, vec->get(row_idx, 0), num_eles);
		else
			repeat_ele(arr, vec->get(col_idx, 0), num_eles);
	}

	virtual const scalar_type &get_type() const {
		return vec->get_type();
	}
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr(new repeat_ele_operate(vec, !in_row));
	}
};

}

dense_matrix::ptr dense_matrix::create_repeat(col_vec::ptr vec, size_t nrow,
		size_t ncol, matrix_layout_t layout, bool byrow, int num_nodes)
{
	if (nrow == 0 || ncol == 0) {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't create a matrix with 0 rows/cols";
		return dense_matrix::ptr();
	}

	size_t long_dim = std::max(nrow, ncol);
	size_t short_dim = std::min(nrow, ncol);
	// We don't want to create a small block matrix.
	if (long_dim > detail::EM_matrix_store::CHUNK_SIZE
			&& short_dim > matrix_conf.get_block_size())
		return block_matrix::create_repeat_layout(vec, nrow, ncol, layout,
				matrix_conf.get_block_size(), byrow, num_nodes);

	if (byrow && vec->get_length() != ncol) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't repeat a vector whose length doesn't match matrix width";
		return dense_matrix::ptr();
	}
	else if (!byrow && vec->get_length() != nrow) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't repeat a vector whose length doesn't match matrix height";
		return dense_matrix::ptr();
	}

	// This is a tall matrix and we repeat the vector by rows.
	if (byrow && nrow > ncol) {
		dense_matrix::ptr tmp = vec->conv_store(true, -1);
		detail::mem_col_matrix_store::const_ptr vec
			= std::dynamic_pointer_cast<const detail::mem_col_matrix_store>(
					tmp->get_raw_store());
		assert(vec);
		set_operate::const_ptr row_op(new repeat_operate(vec));
		set_operate::const_ptr col_op(new repeat_ele_operate(vec, false));
		return dense_matrix::create(detail::set_data_matrix_store::create(
					row_op, col_op, nrow, ncol, layout, num_nodes));
	}
	// This is a wide matrix and we repeat the vector by rows.
	else if (byrow) {
		detail::matrix_store::const_ptr store = vec->get_raw_store();
		assert(store->get_num_cols() == 1);
		std::vector<detail::matrix_store::const_ptr> rows(nrow,
				store->transpose());
		return dense_matrix::create(detail::combined_matrix_store::create(rows,
					layout));
	}
	// This is a tall matrix and we repeat the vector by columns.
	else if (nrow > ncol) {
		detail::matrix_store::const_ptr store = vec->get_raw_store();
		assert(store->get_num_cols() == 1);
		std::vector<detail::matrix_store::const_ptr> cols(ncol, store);
		return dense_matrix::create(detail::combined_matrix_store::create(cols,
					layout));
	}
	// This is a wide matrix and we repeat the vector by columns.
	else {
		dense_matrix::ptr tmp = vec->conv_store(true, -1);
		detail::mem_col_matrix_store::const_ptr vec
			= std::dynamic_pointer_cast<const detail::mem_col_matrix_store>(
					tmp->get_raw_store());
		assert(vec);
		set_operate::const_ptr row_op(new repeat_ele_operate(vec, true));
		set_operate::const_ptr col_op(new repeat_operate(vec));
		return dense_matrix::create(detail::set_data_matrix_store::create(
					row_op, col_op, nrow, ncol, layout, num_nodes));
	}
}

dense_matrix::ptr dense_matrix::create_seq(scalar_variable::ptr start,
		scalar_variable::ptr stride, size_t nrow, size_t ncol,
		matrix_layout_t layout, bool byrow, int num_nodes, bool in_mem,
		safs::safs_file_group::ptr group)
{
	if (nrow == 0 || ncol == 0) {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't create a matrix with 0 rows/cols";
		return dense_matrix::ptr();
	}

	size_t long_dim = std::max(nrow, ncol);
	size_t short_dim = std::min(nrow, ncol);
	const scalar_type &type = start->get_type();
	// We don't want to create a small block matrix.
	if (long_dim < detail::EM_matrix_store::CHUNK_SIZE
			|| short_dim <= matrix_conf.get_block_size()) {
		auto row_op = type.get_set_seq(*start, *stride, nrow, ncol, byrow,
				matrix_layout_t::L_ROW);
		auto col_op = type.get_set_seq(*start, *stride, nrow, ncol, byrow,
				matrix_layout_t::L_COL);
		detail::matrix_store::ptr store = detail::set_data_matrix_store::create(
				row_op, col_op, nrow, ncol, layout, num_nodes);
		return dense_matrix::ptr(new dense_matrix(store));
	}
	else
		return block_matrix::create_seq_layout(start, stride, nrow, ncol, layout,
				matrix_conf.get_block_size(), byrow, num_nodes, in_mem, group);
}

dense_matrix::ptr dense_matrix::create(size_t nrow, size_t ncol,
		matrix_layout_t layout, const scalar_type &type, const set_operate &op,
		int num_nodes, bool in_mem, safs::safs_file_group::ptr group)
{
	if (nrow == 0 || ncol == 0) {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't create a matrix with 0 rows/cols";
		return dense_matrix::ptr();
	}

	size_t long_dim = std::max(nrow, ncol);
	size_t short_dim = std::min(nrow, ncol);
	// We don't want to create a small block matrix.
	if (long_dim < detail::EM_matrix_store::CHUNK_SIZE
			|| short_dim <= matrix_conf.get_block_size()) {
		detail::matrix_store::ptr store = detail::matrix_store::create(
				nrow, ncol, layout, type, num_nodes, in_mem, group);
		if (store == NULL)
			return dense_matrix::ptr();

		store->set_data(op);
		return dense_matrix::ptr(new dense_matrix(store));
	}
	else
		return block_matrix::create_layout(nrow, ncol, layout,
				matrix_conf.get_block_size(), type, op, num_nodes,
				in_mem, group);
}

dense_matrix::ptr dense_matrix::create(data_frame::const_ptr df)
{
	std::vector<detail::matrix_store::const_ptr> mats(df->get_num_vecs());
	for (size_t i = 0; i < mats.size(); i++) {
		detail::vec_store::const_ptr v = df->get_vec(i);
		mats[i] = v->conv2mat(v->get_length(), 1, false);
	}
	return dense_matrix::create(detail::combined_matrix_store::create(mats,
				matrix_layout_t::L_COL));
}

dense_matrix::ptr dense_matrix::transpose() const
{
	return dense_matrix::ptr(new dense_matrix(get_data().transpose()));
}

////////////////////////////// Get rows/cols //////////////////////////////////

namespace
{

/*
 * This gets rows from a wide matrix.
 */
class get_rows_op: public detail::portion_mapply_op
{
	std::vector<off_t> row_idxs;
public:
	get_rows_op(const std::vector<off_t> &row_idxs, size_t num_cols,
			const scalar_type &type): detail::portion_mapply_op(row_idxs.size(),
				num_cols, type) {
		this->row_idxs = row_idxs;
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		// We can only resize the number of the cols in this operation.
		return local_num_rows == get_out_num_rows();
	}

	virtual portion_mapply_op::const_ptr transpose() const;

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("get_rows(") + mats[0]->get_name() + ")";
	}
};

/*
 * This gets cols from a tall matrix.
 */
class get_cols_op: public detail::portion_mapply_op
{
	std::vector<off_t> col_idxs;
public:
	get_cols_op(const std::vector<off_t> &col_idxs, size_t num_rows,
			const scalar_type &type): detail::portion_mapply_op(num_rows,
				col_idxs.size(), type) {
		this->col_idxs = col_idxs;
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		// We can only resize the number of the rows in this operation.
		return local_num_cols == get_out_num_cols();
	}

	virtual portion_mapply_op::const_ptr transpose() const;

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("get_cols(") + mats[0]->get_name() + ")";
	}
};

void get_rows_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 1);
	assert(out.get_num_rows() == row_idxs.size());
	assert(ins[0]->store_layout() == out.store_layout());
	// If the data is layout in row major, we can get the rows easily.
	if (out.store_layout() == matrix_layout_t::L_ROW) {
		const detail::local_row_matrix_store &row_in
			= static_cast<const detail::local_row_matrix_store &>(*ins[0]);
		detail::local_row_matrix_store &row_out
			= static_cast<detail::local_row_matrix_store &>(out);
		for (size_t i = 0; i < row_idxs.size(); i++)
			memcpy(row_out.get_row(i), row_in.get_row(row_idxs[i]),
					out.get_num_cols() * out.get_entry_size());
	}
	// Otherwise, we need to get the wanted elements individually.
	else {
		detail::local_col_matrix_store &col_out
			= static_cast<detail::local_col_matrix_store &>(out);
		const scatter_gather &sg = out.get_type().get_sg();
		std::vector<const char *> col_src(row_idxs.size());
		for (size_t i = 0; i < out.get_num_cols(); i++) {
			for (size_t j = 0; j < col_src.size(); j++)
				col_src[j] = ins[0]->get(row_idxs[j], i);
			sg.gather(col_src, col_out.get_col(i));
		}
	}
}

void get_cols_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 1);
	assert(out.get_num_cols() == col_idxs.size());
	assert(ins[0]->store_layout() == out.store_layout());
	// If the data is layout in col major, we can get cols easily.
	if (out.store_layout() == matrix_layout_t::L_COL) {
		const detail::local_col_matrix_store &col_in
			= static_cast<const detail::local_col_matrix_store &>(*ins[0]);
		detail::local_col_matrix_store &col_out
			= static_cast<detail::local_col_matrix_store &>(out);
		for (size_t i = 0; i < col_idxs.size(); i++)
			memcpy(col_out.get_col(i), col_in.get_col(col_idxs[i]),
					out.get_num_rows() * out.get_entry_size());
	}
	// Otherwise, we need to get the wanted elements indivudally.
	else {
		detail::local_row_matrix_store &row_out
			= static_cast<detail::local_row_matrix_store &>(out);
		const scatter_gather &sg = out.get_type().get_sg();
		std::vector<const char *> row_src(col_idxs.size());
		for (size_t i = 0; i < out.get_num_rows(); i++) {
			for (size_t j = 0; j < row_src.size(); j++)
				row_src[j] = ins[0]->get(i, col_idxs[j]);
			sg.gather(row_src, row_out.get_row(i));
		}
	}
}

detail::portion_mapply_op::const_ptr get_rows_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new get_cols_op(row_idxs,
				get_out_num_cols(), get_output_type()));
}

detail::portion_mapply_op::const_ptr get_cols_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new get_rows_op(col_idxs,
				get_out_num_rows(), get_output_type()));
}

}

scalar_variable::ptr dense_matrix::get(off_t row_idx, off_t col_idx) const
{
	// This is a little dangerous to do. The matrix might be an EM matrix.
	// Pulling data to memory might fill all memory.
	// Normally, I assume people only access individual elements in a large
	// matrix.
	if (is_virtual() || !is_in_mem())
		move_store(true, -1);

	auto mem_store = std::dynamic_pointer_cast<const detail::mem_matrix_store>(
			get_raw_store());
	if (mem_store == NULL)
		return scalar_variable::ptr();
	scalar_variable::ptr tmp = get_type().create_scalar();
	memcpy(tmp->get_raw(), mem_store->get(row_idx, col_idx),
			get_type().get_size());
	return tmp;
}

dense_matrix::ptr dense_matrix::get_col(off_t idx) const
{
	if (idx < 0 || (size_t) idx >= get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "the col index is out of bound";
		return dense_matrix::ptr();
	}

	dense_matrix::ptr t = transpose();
	assert(t);
	return t->get_row(idx);
}

dense_matrix::ptr dense_matrix::get_row(off_t idx) const
{
	if (idx < 0 || (size_t) idx >= get_num_rows()) {
		BOOST_LOG_TRIVIAL(error) << "the row index is out of bound";
		return dense_matrix::ptr();
	}

	std::vector<off_t> idx_vec(1, idx);
	// Different matrix stores have their own best way of getting rows.
	// If a matrix store can't get rows in an optimal way, we can fall back
	// to the default solution.
	detail::matrix_store::const_ptr ret = get_data().get_rows(idx_vec);
	if (ret)
		return col_vec::create(ret);
	// This is the default solution. We construct a virtual matrix that
	// represents the required rows. However, the default solution only
	// works for wide matrices.
	else if (is_wide()) {
		std::vector<detail::matrix_store::const_ptr> ins(1, store);
		get_rows_op::const_ptr op(new get_rows_op(idx_vec, get_num_cols(),
					get_type()));
		ret = __mapply_portion_virtual(ins, op, store_layout());
		assert(ret);
		return col_vec::create(ret);
	}
	else
		return dense_matrix::ptr();
}

dense_matrix::ptr dense_matrix::get_cols(const std::vector<off_t> &idxs) const
{
	if (idxs.empty()) {
		BOOST_LOG_TRIVIAL(error) << "cannot get 0 cols";
		return dense_matrix::ptr();
	}

	for (size_t i = 0; i < idxs.size(); i++)
		if (idxs[i] < 0 || (size_t) idxs[i] >= get_num_cols()) {
			BOOST_LOG_TRIVIAL(error) << "the col index is out of bound";
			return dense_matrix::ptr();
		}

	dense_matrix::ptr t = transpose();
	assert(t);
	dense_matrix::ptr ret = t->get_rows(idxs);
	if (ret == NULL)
		return dense_matrix::ptr();
	return ret->transpose();
}

dense_matrix::ptr dense_matrix::get_rows(const std::vector<off_t> &idxs) const
{
	if (idxs.empty()) {
		BOOST_LOG_TRIVIAL(error) << "cannot get 0 rows";
		return dense_matrix::ptr();
	}

	for (size_t i = 0; i < idxs.size(); i++)
		if (idxs[i] < 0 || (size_t) idxs[i] >= get_num_rows()) {
			BOOST_LOG_TRIVIAL(error) << "the row index is out of bound";
			return dense_matrix::ptr();
		}

	// Different matrix stores have their own best way of getting rows.
	// If a matrix store can't get rows in an optimal way, we can fall back
	// to the default solution.
	detail::matrix_store::const_ptr ret = get_data().get_rows(idxs);
	if (ret)
		return dense_matrix::ptr(new dense_matrix(ret));
	// This is the default solution. We construct a virtual matrix that
	// represents the required rows. However, the default solution only
	// works for wide matrices.
	else if (is_wide()) {
		std::vector<detail::matrix_store::const_ptr> ins(1, store);
		get_rows_op::const_ptr op(new get_rows_op(idxs, get_num_cols(),
					get_type()));
		detail::matrix_store::ptr res = __mapply_portion_virtual(ins, op,
				store_layout());
		assert(res);
		return dense_matrix::create(res);
	}
	else
		return dense_matrix::ptr();
}

dense_matrix::ptr dense_matrix::get_cols(size_t start, size_t end,
		size_t step) const
{
	if (end > get_num_cols() || start >= get_num_cols() || start >= end) {
		BOOST_LOG_TRIVIAL(error) << "column index is out of the range";
		return dense_matrix::ptr();
	}

	std::vector<off_t> col_idxs((end - start) / step);
	for (size_t i = 0; i < col_idxs.size(); i++)
		col_idxs[i] = start + i * step;
	// TODO we need optimize this.
	return get_cols(col_idxs);
}

dense_matrix::ptr dense_matrix::get_rows(size_t start, size_t end,
		size_t step) const
{
	if (end > get_num_rows() || start >= get_num_rows() || start >= end) {
		BOOST_LOG_TRIVIAL(error) << "row index is out of the range";
		return dense_matrix::ptr();
	}

	std::vector<off_t> row_idxs((end - start) / step);
	for (size_t i = 0; i < row_idxs.size(); i++)
		row_idxs[i] = start + i * step;
	// TODO we need optimize this.
	return get_rows(row_idxs);
}

namespace
{

/*
 * This gets rows from a small matrix and outputs a tall matrix.
 */
class repeat_rows_op: public detail::portion_mapply_op
{
	detail::mem_row_matrix_store::const_ptr store;
	volatile bool success;
public:
	repeat_rows_op(detail::mem_row_matrix_store::const_ptr store,
			size_t tot_nrow): detail::portion_mapply_op(tot_nrow,
				store->get_num_cols(), store->get_type()) {
		this->store = store;
		this->success = true;
		assert(tot_nrow > store->get_num_cols());
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		// We can only resize the number of the rows in this operation.
		return local_num_cols == get_out_num_cols();
	}

	virtual portion_mapply_op::const_ptr transpose() const;

	virtual bool is_success() const {
		return success;
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("repeat_rows(") + mats[0]->get_name() + ")";
	}
};

/*
 * This gets cols from a small matrix and outputs a wide matrix.
 */
class repeat_cols_op: public detail::portion_mapply_op
{
	detail::mem_col_matrix_store::const_ptr store;
	volatile bool success;
public:
	repeat_cols_op(detail::mem_col_matrix_store::const_ptr store,
			size_t tot_ncol): detail::portion_mapply_op(store->get_num_rows(),
				tot_ncol, store->get_type()) {
		this->store = store;
		this->success = true;
		assert(store->get_num_rows() < tot_ncol);
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		// We can only resize the number of the cols in this operation.
		return local_num_rows == get_out_num_rows();
	}

	virtual portion_mapply_op::const_ptr transpose() const;

	virtual bool is_success() const {
		return success;
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("repeat_cols(") + mats[0]->get_name() + ")";
	}
};

void repeat_rows_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	// The operation might fail in one portion.
	// If it does, we don't need to perform computation on other portions.
	if (!success)
		return;

	assert(ins[0]->get_num_cols() == 1);
	assert(out.get_num_cols() == store->get_num_cols());
	assert(out.store_layout() == matrix_layout_t::L_ROW);
	detail::local_row_matrix_store &row_out
		= static_cast<detail::local_row_matrix_store &>(out);
	const size_t *idx
		= reinterpret_cast<const size_t *>(ins[0]->get_raw_arr());
	for (size_t i = 0; i < ins[0]->get_num_rows(); i++) {
		if (idx[i] >= store->get_num_rows()) {
			BOOST_LOG_TRIVIAL(error) << boost::format(
					"index (%1%) exceeds the number of rows (%2%) in the matrix")
				% idx[i] % store->get_num_rows();
			const_cast<repeat_rows_op *>(this)->success = false;
			break;
		}
		memcpy(row_out.get_row(i), store->get_row(idx[i]),
				out.get_num_cols() * out.get_entry_size());
	}
}

void repeat_cols_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	// The operation might fail in one portion.
	// If it does, we don't need to perform computation on other portions.
	if (!success)
		return;

	assert(ins[0]->get_num_rows() == 1);
	assert(out.get_num_rows() == store->get_num_rows());
	assert(out.store_layout() == matrix_layout_t::L_COL);
	detail::local_col_matrix_store &col_out
		= static_cast<detail::local_col_matrix_store &>(out);
	const size_t *idx
		= reinterpret_cast<const size_t *>(ins[0]->get_raw_arr());
	for (size_t i = 0; i < ins[0]->get_num_cols(); i++) {
		if (idx[i] >= store->get_num_cols()) {
			BOOST_LOG_TRIVIAL(error)
				<< "the index exceeds the number of cols in the matrix";
			const_cast<repeat_cols_op *>(this)->success = false;
			break;
		}
		memcpy(col_out.get_col(i), store->get_col(idx[i]),
				out.get_num_rows() * out.get_entry_size());
	}
}

detail::portion_mapply_op::const_ptr repeat_rows_op::transpose() const
{
	detail::matrix_store::const_ptr t = store->transpose();
	detail::mem_col_matrix_store::const_ptr col_store
		= std::dynamic_pointer_cast<const detail::mem_col_matrix_store>(t);
	assert(col_store);
	return portion_mapply_op::const_ptr(new repeat_cols_op(col_store,
				get_out_num_rows()));
}

detail::portion_mapply_op::const_ptr repeat_cols_op::transpose() const
{
	detail::matrix_store::const_ptr t = store->transpose();
	detail::mem_row_matrix_store::const_ptr row_store
		= std::dynamic_pointer_cast<const detail::mem_row_matrix_store>(t);
	assert(row_store);
	return portion_mapply_op::const_ptr(new repeat_rows_op(row_store,
				get_out_num_cols()));
}

/*
 * This portion operation is designed to get rows from a tall matrix.
 * We can assume the matrix is stored in row major order.
 */
class get_row_portion_op: public detail::portion_mapply_op
{
	detail::matrix_store::ptr res;
	detail::matrix_append::ptr append;
	// this determines the size of a portion read from the input matrix
	// each time. We need this to determine the relative location of
	// the selected rows in the output matrix.
	size_t portion_size;
public:
	typedef std::shared_ptr<const get_row_portion_op> const_ptr;

	get_row_portion_op(detail::matrix_store::ptr res,
			size_t portion_size): portion_mapply_op(res->get_num_rows(),
				res->get_num_cols(), res->get_type()) {
		this->res = res;
		this->append = detail::matrix_append::create(res);
		this->portion_size = portion_size;
	}

	detail::matrix_store::const_ptr get_result() const {
		append->flush();
		if (res->is_wide())
			res->resize(res->get_num_rows(),
					append->get_written_eles() / res->get_num_rows());
		else
			res->resize(append->get_written_eles() / res->get_num_cols(),
					res->get_num_cols());
		return res;
	}

	// We don't need this method.
	virtual detail::portion_mapply_op::const_ptr transpose() const {
		return detail::portion_mapply_op::const_ptr();
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return std::string("get_row_bool(") + mats[0]->get_name() + ")";
	}

	/*
	 * Give a hint if this operation is aggregation, so we can optimize
	 * the backend accordingly. When this is an aggregation operation,
	 * the second `run' method has to be implemented.
	 */
	virtual bool is_agg() const {
		return false;
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		return true;
	}
};

void get_row_portion_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins) const
{
	// The first input matrix is the matrix where we get rows from.
	// The second input matrix is a boolean vector that indicates which
	// rows should be selected.
	assert(ins.size() == 2);
	assert(ins[1]->get_num_cols() == 1);
	assert(ins[1]->get_type() == get_scalar_type<bool>());
	assert(ins[0]->store_layout() == matrix_layout_t::L_ROW);
	assert(ins[1]->store_layout() == matrix_layout_t::L_COL);
	detail::local_row_matrix_store::const_ptr data_store
		= std::static_pointer_cast<const detail::local_row_matrix_store>(ins[0]);
	detail::local_col_matrix_store::const_ptr idx_store
		= std::static_pointer_cast<const detail::local_col_matrix_store>(ins[1]);
	const bool *bool_idxs
		= reinterpret_cast<const bool *>(idx_store->get_col(0));
	std::vector<const char *> selected_rows;
	for (size_t i = 0; i < idx_store->get_num_rows(); i++)
		if (bool_idxs[i])
			selected_rows.push_back(data_store->get_row(i));
	if (selected_rows.empty())
		// If there isn't data we need to write, we still need to inform
		// matrix append of data from this portion. Otherwise, data from portions
		// behind it won't be appended.
		append->write_async(NULL,
				// Here we generate contiguous sequence numbers.
				data_store->get_global_start_row() / portion_size);
	else {
		// We don't want the data in the local matrix to be cached because
		// the size of this matrix is irregular.
		detail::local_row_matrix_store::ptr out(new detail::local_buf_row_matrix_store(
					0, 0, selected_rows.size(), data_store->get_num_cols(),
					data_store->get_type(), -1, false));
		for (size_t i = 0; i < selected_rows.size(); i++)
			memcpy(out->get_row(i), selected_rows[i],
					out->get_num_cols() * out->get_type().get_size());
		assert(data_store->get_global_start_row() % portion_size == 0);
		append->write_async(out,
				// Here we generate contiguous sequence numbers.
				data_store->get_global_start_row() / portion_size);
	}
}

}

dense_matrix::ptr dense_matrix::get_cols(col_vec::ptr idxs) const
{
	if (idxs->get_length() == 0) {
		BOOST_LOG_TRIVIAL(error) << "cannot get 0 cols";
		return dense_matrix::ptr();
	}

	dense_matrix::ptr t = transpose();
	assert(t);
	dense_matrix::ptr ret = t->get_rows(idxs);
	if (ret == NULL)
		return dense_matrix::ptr();
	return ret->transpose();
}

static std::vector<detail::mem_row_matrix_store::const_ptr> split_mat_vertial(
		detail::matrix_store::const_ptr store, size_t block_size)
{
	size_t num_blocks = div_ceil<size_t>(store->get_num_cols(), block_size);
	std::vector<detail::mem_row_matrix_store::const_ptr> blocks(num_blocks);
	for (size_t i = 0; i < num_blocks; i++) {
		size_t sub_num_cols = std::min(block_size,
				store->get_num_cols() - i * block_size);
		detail::mem_row_matrix_store::ptr row_store
			= detail::mem_row_matrix_store::create(store->get_num_rows(),
					sub_num_cols, store->get_type());
		detail::local_matrix_store::const_ptr src_part
			= store->get_portion(0, i * block_size, store->get_num_rows(),
					sub_num_cols);
		detail::local_matrix_store::ptr dst_part = row_store->get_portion(0, 0,
				row_store->get_num_rows(), row_store->get_num_cols());
		dst_part->copy_from(*src_part);
		blocks[i] = row_store;
	}
	return blocks;
}

dense_matrix::ptr dense_matrix::get_rows_bool(col_vec::ptr idxs) const
{
	if (get_num_rows() != idxs->get_length()) {
		BOOST_LOG_TRIVIAL(error)
			<< "get rows: #rows doesn't match length of boolean vector";
		return dense_matrix::ptr();
	}

	// If the matrix is wide, we can explicitly pick the rows.
	// But we need to figure out first which rows we need to pick.
	if (is_wide()) {
		std::vector<bool> bool_idxs = idxs->conv2std<bool>();
		std::vector<off_t> offs;
		for (size_t i = 0; i < bool_idxs.size(); i++)
			if (bool_idxs[i])
				offs.push_back(i);
		return get_rows(offs);
	}

	std::vector<detail::matrix_store::const_ptr> in_mats(2);
	dense_matrix::ptr tmp = conv2(matrix_layout_t::L_ROW);
	in_mats[0] = tmp->get_raw_store();
	tmp = idxs->conv2(matrix_layout_t::L_COL);
	in_mats[1] = tmp->get_raw_store();
	std::vector<detail::matrix_store::ptr> out_mats;

	// TODO let's assume the filter results are always in memory.
	detail::matrix_store::ptr res = detail::matrix_store::create(
			get_num_rows(), get_num_cols(), matrix_layout_t::L_ROW,
			get_type(), store->get_num_nodes(), true);
	get_row_portion_op::const_ptr op(new get_row_portion_op(res,
				store->get_portion_size().first));
	bool ret = __mapply_portion(in_mats, op, out_mats);
	if (!ret)
		return dense_matrix::ptr();
	else
		return dense_matrix::create(op->get_result());
}

dense_matrix::ptr dense_matrix::get_rows(col_vec::ptr idxs) const
{
	if (idxs->get_length() == 0) {
		BOOST_LOG_TRIVIAL(error) << "cannot get 0 rows";
		return dense_matrix::ptr();
	}

	if (idxs->get_type() == get_scalar_type<bool>())
		return get_rows_bool(idxs);

	// In this case, we just read the rows from the current matrix physically
	// and outputs a materialized matrix.
	// If the output matrix is still wide
	if (idxs->get_length() < get_num_cols()
			// or the number of output rows is smaller than the number of
			// original rows.
			|| idxs->get_length() < get_num_rows()) {
		dense_matrix::ptr tmp = idxs->cast_ele_type(get_scalar_type<off_t>());
		tmp = tmp->conv2(matrix_layout_t::L_COL);
		tmp = tmp->conv_store(true, -1);
		detail::mem_col_matrix_store::const_ptr col_store
			= std::dynamic_pointer_cast<const detail::mem_col_matrix_store>(
					tmp->get_raw_store());
		assert(col_store);
		std::vector<off_t> std_idxs(idxs->get_length());
		memcpy(std_idxs.data(), col_store->get_col(0),
				idxs->get_length() * sizeof(off_t));
		return get_rows(std_idxs);
	}

	// In this case, the output matrix is a tall matrix and the number of
	// output rows is larger than the original matrix.
	// We use a virtual matrix to represent the output matrix.

	// If the matrix is very skinny or the output matrix is still small enough.
	if (get_num_cols() <= matrix_conf.get_block_size()
				|| idxs->get_length() <= detail::EM_matrix_store::CHUNK_SIZE) {
		dense_matrix::ptr tmp = conv2(matrix_layout_t::L_ROW);
		tmp = tmp->conv_store(true, -1);
		detail::mem_row_matrix_store::const_ptr row_store
			= std::dynamic_pointer_cast<const detail::mem_row_matrix_store>(
					tmp->get_raw_store());
		assert(row_store);

		std::vector<dense_matrix::const_ptr> mats(1,
				col_vec::create(idxs->cast_ele_type(get_scalar_type<off_t>())));
		detail::portion_mapply_op::const_ptr op(new repeat_rows_op(row_store,
					idxs->get_length()));
		return detail::mapply_portion(mats, op, matrix_layout_t::L_ROW);
	}
	// In this case, we output a wide and a large matrix.
	else {
		dense_matrix::ptr tmp = conv2(matrix_layout_t::L_ROW);
		tmp = tmp->conv_store(true, -1);
		std::vector<detail::mem_row_matrix_store::const_ptr> blocks
			= split_mat_vertial(tmp->get_raw_store(), matrix_conf.get_block_size());

		std::vector<detail::matrix_store::const_ptr> out_stores(blocks.size());
		std::vector<dense_matrix::const_ptr> mats(1,
				col_vec::create(idxs->cast_ele_type(get_scalar_type<off_t>())));
		for (size_t i = 0; i < blocks.size(); i++) {
			detail::portion_mapply_op::const_ptr op(new repeat_rows_op(blocks[i],
						idxs->get_length()));
			dense_matrix::ptr tmp = detail::mapply_portion(mats, op,
					matrix_layout_t::L_ROW);
			out_stores[i] = tmp->get_raw_store();
		}
		return block_matrix::create(detail::combined_matrix_store::create(
					out_stores, matrix_layout_t::L_ROW));
	}
}

col_vec::ptr dense_matrix::get_eles(dense_matrix::ptr idx) const
{
	if (idx->get_num_cols() != 2) {
		BOOST_LOG_TRIVIAL(error) << "The index matrix has to have two columns";
		return col_vec::ptr();
	}

	// This is useful to get a small number of elements from a matrix.
	// TODO we need to optimize for many other cases.
	detail::mem_matrix_store::const_ptr mem_store
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				get_raw_store());
	if (mem_store == NULL) {
		BOOST_LOG_TRIVIAL(error) << "Can't get elements from an EM matrix";
		return col_vec::ptr();
	}

	idx = idx->cast_ele_type(get_scalar_type<off_t>());
	idx = idx->conv2(matrix_layout_t::L_COL);
	idx = idx->conv_store(true, -1);
	detail::mem_matrix_store::const_ptr idx_store
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				idx->get_raw_store());
	assert(idx_store);

	detail::mem_matrix_store::ptr vec = detail::mem_matrix_store::create(
			idx->get_num_rows(), 1, matrix_layout_t::L_COL, get_type(), -1);
	size_t entry_size = mem_store->get_type().get_size();
	for (size_t i = 0; i < idx_store->get_num_rows(); i++) {
		const off_t *ridx_arr
			= reinterpret_cast<const off_t *>(idx_store->get_col(0));
		const off_t *cidx_arr
			= reinterpret_cast<const off_t *>(idx_store->get_col(1));
		assert(ridx_arr);
		assert(cidx_arr);
		off_t row_idx = ridx_arr[i];
		off_t col_idx = cidx_arr[i];
		memcpy(vec->get(i, 0), mem_store->get(row_idx, col_idx), entry_size);
	}
	return col_vec::create(vec);
}

////////////////////////////// Inner product //////////////////////////////////

dense_matrix::ptr dense_matrix::inner_prod(const dense_matrix &m,
		bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
		matrix_layout_t out_layout) const
{
	if (!verify_inner_prod(m, *left_op, *right_op))
		return dense_matrix::ptr();

	// If input matrices are sink matrices, we materialize them first.
	if (store->is_sink())
		materialize_self();
	if (m.store->is_sink())
		m.materialize_self();

	size_t long_dim1 = std::max(get_num_rows(), get_num_cols());
	size_t long_dim2 = std::max(m.get_num_rows(), m.get_num_cols());
	// We prefer to perform computation on the larger matrix.
	// If the matrix in the right operand is larger, we transpose
	// the entire computation.
	if (long_dim2 > long_dim1) {
		dense_matrix::ptr t_mat1 = this->transpose();
		dense_matrix::ptr t_mat2 = m.transpose();
		matrix_layout_t t_layout = out_layout;
		if (t_layout == matrix_layout_t::L_ROW)
			t_layout = matrix_layout_t::L_COL;
		else if (t_layout == matrix_layout_t::L_COL)
			t_layout = matrix_layout_t::L_ROW;
		dense_matrix::ptr t_res = t_mat2->inner_prod(*t_mat1,
				left_op, right_op, t_layout);
		return t_res->transpose();
	}

	if (is_wide())
		return inner_prod_wide(m, left_op, right_op, out_layout);
	else
		return inner_prod_tall(m, left_op, right_op, out_layout);
}

namespace
{

class inner_prod_tall_op: public detail::portion_mapply_op
{
	// I need to keep the right matrix to make sure the memory isn't deallocated.
	detail::matrix_store::const_ptr right;
	detail::local_matrix_store::const_ptr local_right;
	bulk_operate::const_ptr left_op;
	bulk_operate::const_ptr right_op;
public:
	inner_prod_tall_op(detail::matrix_store::const_ptr right,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			size_t out_num_rows, size_t out_num_cols): detail::portion_mapply_op(
				out_num_rows, out_num_cols, right_op->get_output_type()) {
		this->right = right;
		// We assume the right matrix is small, so we don't need to partition it.
		this->local_right = right->get_portion(0);
		assert(local_right->get_num_rows() == right->get_num_rows()
				&& local_right->get_num_cols() == right->get_num_cols());
		this->left_op = left_op;
		this->right_op = right_op;
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		// We can only resize the number of the rows in this operation.
		return local_num_cols == get_out_num_cols();
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;
	virtual portion_mapply_op::const_ptr transpose() const;
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("inner_prod(") + mats[0]->get_name() + ","
			+ local_right->get_name() + ")";
	}
};

class t_inner_prod_tall_op: public detail::portion_mapply_op
{
	inner_prod_tall_op op;
public:
	t_inner_prod_tall_op(const inner_prod_tall_op &_op): detail::portion_mapply_op(
			_op.get_out_num_cols(), _op.get_out_num_rows(),
			_op.get_output_type()), op(_op) {
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		return op.is_resizable(local_start_col, local_start_row,
				local_num_cols, local_num_rows);
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		assert(ins.size() == 1);
		std::vector<fm::detail::local_matrix_store::const_ptr> t_ins(ins.size());
		t_ins[0] = std::static_pointer_cast<const fm::detail::local_matrix_store>(
				ins[0]->transpose());
		fm::detail::local_matrix_store::ptr t_out
			= std::static_pointer_cast<fm::detail::local_matrix_store>(
					out.transpose());
		op.run(t_ins, *t_out);
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		return fm::detail::portion_mapply_op::const_ptr(
				new inner_prod_tall_op(op));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return op.to_string(mats);
	}
};

detail::portion_mapply_op::const_ptr inner_prod_tall_op::transpose() const
{
	return fm::detail::portion_mapply_op::const_ptr(
			new t_inner_prod_tall_op(*this));
}

void inner_prod_tall_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 1);
	assert(ins[0]->get_global_start_col() == out.get_global_start_col());
	assert(ins[0]->get_global_start_row() == out.get_global_start_row());
	out.reset_data();
	detail::inner_prod_tall(*ins[0], *local_right, *left_op, *right_op, out);
}

}

// A flag to indicate whether we need to convert the matrix layout.
// This flag is only used for testing.
bool inner_prod_conv = true;

dense_matrix::ptr dense_matrix::inner_prod_tall(
		const dense_matrix &m, bulk_operate::const_ptr left_op,
		bulk_operate::const_ptr right_op, matrix_layout_t out_layout) const
{
	detail::matrix_store::const_ptr right = m.get_raw_store();
	// If the left matrix is row-major, the right matrix should be
	// column-major. When the left matrix is tall, the right matrix should
	// be small. It makes sense to convert the right matrix to column major
	// before we break up the left matrix for parallel processing.
	if (!is_wide() && this->store_layout() == matrix_layout_t::L_ROW) {
		dense_matrix::ptr tmp = m.conv2(matrix_layout_t::L_COL);
		right = tmp->get_raw_store();
	}
	// TODO the right matrix might be a sparse matrix. conv_store
	// will turn it into a dense matrix.
	// We should optimize for the right sparse matrix later.
	if (right->is_virtual() || !right->is_in_mem() || right->get_num_nodes() > 0) {
		dense_matrix::ptr tmp = dense_matrix::create(right);
		tmp = tmp->conv_store(true, -1);
		right = tmp->get_raw_store();
	}
	assert(right->is_in_mem());
	assert(right->get_num_nodes() == -1);
	assert(!right->is_virtual());

	std::vector<detail::matrix_store::const_ptr> ins(1);
	// If the matrix is already stored in col-major order or this is
	// a relatively wide matrix, we don't need to convert its data layout.
	if (store_layout() == matrix_layout_t::L_COL || get_num_cols() > 32
			|| !inner_prod_conv)
		ins[0] = this->get_raw_store();
	else {
		// If this is a very tall and skinny matrix, col-major works better
		// even though we need to pay the overhead of converting the layout.
		dense_matrix::ptr tmp = conv2(matrix_layout_t::L_COL);
		ins[0] = tmp->get_raw_store();
	}

	if (out_layout == matrix_layout_t::L_NONE) {
		// If the left matrix is col-major, the output matrix should also
		// be col-major.
		if (ins[0]->store_layout() == matrix_layout_t::L_COL)
			out_layout = matrix_layout_t::L_COL;
		else
			// We don't care about the layout of the output matrix in this case.
			out_layout = matrix_layout_t::L_ROW;
	}

	inner_prod_tall_op::const_ptr mapply_op(new inner_prod_tall_op(right,
				left_op, right_op, get_num_rows(), m.get_num_cols()));
	detail::matrix_store::ptr res = __mapply_portion_virtual(ins, mapply_op, out_layout);
	assert(res);
	return dense_matrix::create(res);
}

dense_matrix::ptr dense_matrix::inner_prod_wide(
		const dense_matrix &m, bulk_operate::const_ptr left_op,
		bulk_operate::const_ptr right_op, matrix_layout_t out_layout) const
{
	// If the right matrix is stored in col major order,
	// the left matrix should be stored in row major order.
	detail::matrix_store::const_ptr left_mat;
	if (store_layout() == matrix_layout_t::L_COL
			&& m.store_layout() == matrix_layout_t::L_COL) {
		dense_matrix::ptr tmp = conv2(matrix_layout_t::L_ROW);
		left_mat = tmp->get_raw_store();
	}
	else
		left_mat = get_raw_store();
	assert(left_mat);
	// If the left matrix is stored in row major order,
	// the right matrix should be stored in col major order.
	detail::matrix_store::const_ptr right_mat;
	if (store_layout() == matrix_layout_t::L_ROW
			&& m.store_layout() == matrix_layout_t::L_ROW) {
		dense_matrix::ptr tmp = m.conv2(matrix_layout_t::L_COL);
		right_mat = tmp->get_raw_store();
	}
	else
		right_mat = m.get_raw_store();
	assert(right_mat);

	detail::matrix_store::ptr res = detail::IPW_matrix_store::create(
				left_mat, right_mat, left_op, right_op, out_layout);
	return dense_matrix::create(res);
}

////////////////////////////// Aggregation /////////////////////////////

namespace
{

/*
 * This allows us to aggregate on the shorter dimension.
 * It outputs a long vector, so the result doesn't need to be materialized
 * immediately.
 */
class matrix_short_agg_op: public detail::portion_mapply_op
{
	matrix_margin margin;
	agg_operate::const_ptr op;
public:
	matrix_short_agg_op(matrix_margin margin, agg_operate::const_ptr op,
			size_t out_num_rows, size_t out_num_cols): detail::portion_mapply_op(
				out_num_rows, out_num_cols, op->get_output_type()) {
		this->margin = margin;
		this->op = op;
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		assert(ins.size() == 1);
		// If we perform agg on rows, it implies that this matrix is a tall
		// matrix.
		detail::part_dim_t dim = margin == matrix_margin::MAR_ROW
			? detail::part_dim_t::PART_DIM1 : detail::part_dim_t::PART_DIM2;
		// The output matrix is actually a vector.
		if (out.get_num_rows() == 1) {
			assert(out.store_layout() == matrix_layout_t::L_ROW);
			char *arr = out.get_raw_arr();
			assert(arr);
			detail::local_ref_contig_col_matrix_store ref_out(arr, 0, 0,
					out.get_num_cols(), 1, out.get_type(), out.get_node_id());
			aggregate(*ins[0], *op, margin, dim, ref_out);
		}
		else {
			assert(out.store_layout() == matrix_layout_t::L_COL);
			assert(out.get_num_cols() == 1);
			aggregate(*ins[0], *op, margin, dim, out);
		}
	}

	virtual portion_mapply_op::const_ptr transpose() const {
		matrix_margin new_margin = (this->margin == matrix_margin::MAR_ROW ?
			matrix_margin::MAR_COL : matrix_margin::MAR_ROW);
		return portion_mapply_op::const_ptr(new matrix_short_agg_op(
					new_margin, op, get_out_num_cols(), get_out_num_rows()));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("agg(") + mats[0]->get_name() + ")";
	}
};

}

static detail::matrix_store::ptr aggregate(detail::matrix_store::const_ptr store,
		matrix_margin margin, agg_operate::const_ptr op)
{
	/*
	 * If we aggregate on the shorter dimension.
	 */
	if ((margin == matrix_margin::MAR_ROW && !store->is_wide())
			|| (margin == matrix_margin::MAR_COL && store->is_wide())) {
		// We can turn a wide matrix to a tall matrix to simplify implementation.
		if (margin == matrix_margin::MAR_COL && store->is_wide()) {
			margin = matrix_margin::MAR_ROW;
			store = store->transpose();
		}
		std::vector<detail::matrix_store::const_ptr> ins(1);
		ins[0] = store;
		size_t out_num_rows
			= margin == matrix_margin::MAR_ROW ? store->get_num_rows() : 1;
		size_t out_num_cols
			= margin == matrix_margin::MAR_COL ? store->get_num_cols() : 1;
		matrix_short_agg_op::const_ptr agg_op(new matrix_short_agg_op(
					margin, op, out_num_rows, out_num_cols));
		return __mapply_portion_virtual(ins, agg_op, matrix_layout_t::L_COL);
	}
	if (!op->has_combine()) {
		BOOST_LOG_TRIVIAL(error)
			<< "aggregation on the long dimension requires combine";
		return detail::matrix_store::ptr();
	}

	/*
	 * If we aggregate on the entire matrix or on the longer dimension.
	 */
	return detail::agg_matrix_store::create(store, margin, op);
}

dense_matrix::ptr dense_matrix::aggregate(matrix_margin margin,
			agg_operate::const_ptr op) const
{
	if (this->get_type() != op->get_input_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The matrix element type is different from the operator";
		return dense_matrix::ptr();
	}
	// If the input matrix is a sink matrix, we materialize it first.
	if (store->is_sink())
		materialize_self();
	auto res = fm::aggregate(store, margin, op);
	if (res)
		return dense_matrix::create(res);
	else
		return dense_matrix::ptr();
}

scalar_variable::ptr dense_matrix::aggregate(bulk_operate::const_ptr op) const
{
	return aggregate(agg_operate::create(op));
}

scalar_variable::ptr dense_matrix::aggregate(agg_operate::const_ptr op) const
{
	dense_matrix::ptr agg_mat = aggregate(matrix_margin::BOTH, op);
	if (agg_mat == NULL)
		return scalar_variable::ptr();

	detail::matrix_store::const_ptr _res = agg_mat->get_raw_store();
	assert(_res != NULL);
	// It's evaluated lazily
	detail::virtual_matrix_store::const_ptr virt_res
		= std::dynamic_pointer_cast<const detail::virtual_matrix_store>(_res);
	assert(virt_res);
	_res = virt_res->materialize(true, -1);
	assert(_res != NULL);
	assert(_res->get_num_rows() == 1 && _res->get_num_cols() == 1);
	assert(_res->is_in_mem());

	scalar_variable::ptr res = op->get_output_type().create_scalar();
	res->set_raw(dynamic_cast<const detail::mem_matrix_store &>(
				*_res).get_raw_arr(), res->get_size());
	return res;
}

namespace
{

class cum_short_dim_op: public detail::portion_mapply_op
{
	matrix_margin margin;
	agg_operate::const_ptr op;
public:
	cum_short_dim_op(matrix_margin margin, agg_operate::const_ptr op,
			size_t out_num_rows, size_t out_num_cols): detail::portion_mapply_op(
				out_num_rows, out_num_cols, op->get_output_type()) {
		this->margin = margin;
		this->op = op;
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		if (margin == matrix_margin::MAR_ROW)
			// We can only resize the number of the rows in this operation.
			return local_num_cols == get_out_num_cols();
		else
			// We can only resize the number of the cols in this operation.
			return local_num_rows == get_out_num_rows();
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		detail::part_dim_t dim = margin == matrix_margin::MAR_ROW
			? detail::part_dim_t::PART_DIM1 : detail::part_dim_t::PART_DIM2;
		detail::cum(*ins[0], NULL, *op, margin, dim, out);
	}
	virtual portion_mapply_op::const_ptr transpose() const {
		matrix_margin new_margin = this->margin == matrix_margin::MAR_ROW ?
			matrix_margin::MAR_COL : matrix_margin::MAR_ROW;
		return portion_mapply_op::const_ptr(new cum_short_dim_op(
					new_margin, op, get_out_num_cols(), get_out_num_rows()));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("cum(") + mats[0]->get_name() + ")";
	}
};

}

dense_matrix::ptr dense_matrix::cum(matrix_margin margin,
		agg_operate::const_ptr op) const
{
	if ((margin == matrix_margin::MAR_ROW && !is_wide())
			|| (margin == matrix_margin::MAR_COL && is_wide())) {
		std::vector<detail::matrix_store::const_ptr> ins(1);
		ins[0] = this->get_raw_store();
		cum_short_dim_op::const_ptr apply_op(new cum_short_dim_op(
					margin, op, get_num_rows(), get_num_cols()));
		detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
				apply_op, store_layout());
		return dense_matrix::create(ret);
	}
	else {
		std::vector<detail::matrix_store::const_ptr> ins(1);
		ins[0] = this->get_raw_store();
		auto portion_size = store->get_portion_size();
		detail::cum_long_dim_op::const_ptr apply_op(new detail::cum_long_dim_op(
					margin, op, portion_size.first * portion_size.second,
					get_num_rows(), get_num_cols()));
		detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
				apply_op, store_layout());
		return dense_matrix::create(ret);
	}
}

namespace
{

class matrix_margin_apply_op: public detail::portion_mapply_op
{
	matrix_margin margin;
	arr_apply_operate::const_ptr op;
public:
	matrix_margin_apply_op(matrix_margin margin, arr_apply_operate::const_ptr op,
			size_t out_num_rows, size_t out_num_cols): detail::portion_mapply_op(
				out_num_rows, out_num_cols, op->get_output_type()) {
		this->margin = margin;
		this->op = op;
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		if (margin == matrix_margin::MAR_ROW)
			// We can only resize the number of the rows in this operation.
			return local_num_cols == get_out_num_cols();
		else
			// We can only resize the number of the cols in this operation.
			return local_num_rows == get_out_num_rows();
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		detail::apply(margin, *op, *ins[0], out);
	}
	virtual portion_mapply_op::const_ptr transpose() const {
		matrix_margin new_margin = this->margin == matrix_margin::MAR_ROW ?
			matrix_margin::MAR_COL : matrix_margin::MAR_ROW;
		return portion_mapply_op::const_ptr(new matrix_margin_apply_op(
					new_margin, op, get_out_num_cols(), get_out_num_rows()));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("apply(") + mats[0]->get_name() + ")";
	}
};

}

dense_matrix::ptr dense_matrix::apply(matrix_margin margin,
		arr_apply_operate::const_ptr op) const
{
	// If input matrices are sink matrices, we materialize them first.
	if (store->is_sink())
		materialize_self();

	// In these two cases, we need to convert the matrix store layout
	// before we can apply the function to the matrix.
	detail::matrix_store::const_ptr this_mat;
	if (is_wide() && store_layout() == matrix_layout_t::L_COL
			&& margin == matrix_margin::MAR_ROW) {
		dense_matrix::ptr mat = conv2(matrix_layout_t::L_ROW);
		mat->materialize_self();
		this_mat = mat->get_raw_store();
	}
	else if (!is_wide() && store_layout() == matrix_layout_t::L_ROW
			&& margin == matrix_margin::MAR_COL) {
		dense_matrix::ptr mat = conv2(matrix_layout_t::L_COL);
		mat->materialize_self();
		this_mat = mat->get_raw_store();
	}
	else
		this_mat = get_raw_store();
	assert(this_mat);

	// In these two cases, apply the function on the rows/columns in the long
	// dimension. The previous two cases are handled as one of the two cases
	// after the matrix layout conversion from the previous cases.
	if (is_wide() && this_mat->store_layout() == matrix_layout_t::L_ROW
			&& margin == matrix_margin::MAR_ROW) {
#if 0
		// In this case, it's very difficult to make it work for NUMA matrix.
		assert(get_num_nodes() == -1);
		detail::mem_row_matrix_store::const_ptr row_mat
			= detail::mem_row_matrix_store::cast(this_mat);
		detail::mem_row_matrix_store::ptr ret_mat
			= detail::mem_row_matrix_store::create(get_num_rows(),
					op->get_num_out_eles(), op->get_output_type());
		for (size_t i = 0; i < get_num_rows(); i++) {
			local_cref_vec_store in_vec(row_mat->get_row(i), 0,
					this_mat->get_num_cols(), get_type(), -1);
			local_ref_vec_store out_vec(ret_mat->get_row(i), 0,
					ret_mat->get_num_cols(), ret_mat->get_type(), -1);
			op->run(in_vec, out_vec);
		}
		return mem_dense_matrix::create(ret_mat);
#endif
		BOOST_LOG_TRIVIAL(error)
			<< "it doesn't support to apply rows on a wide matrix";
		return dense_matrix::ptr();
	}
	else if (!is_wide() && this_mat->store_layout() == matrix_layout_t::L_COL
			&& margin == matrix_margin::MAR_COL) {
#if 0
		assert(get_num_nodes() == -1);
		detail::mem_col_matrix_store::const_ptr col_mat
			= detail::mem_col_matrix_store::cast(this_mat);
		detail::mem_col_matrix_store::ptr ret_mat
			= detail::mem_col_matrix_store::create(op->get_num_out_eles(),
					get_num_cols(), op->get_output_type());
		for (size_t i = 0; i < get_num_cols(); i++) {
			local_cref_vec_store in_vec(col_mat->get_col(i), 0,
					this_mat->get_num_rows(), get_type(), -1);
			local_ref_vec_store out_vec(ret_mat->get_col(i), 0,
					ret_mat->get_num_rows(), ret_mat->get_type(), -1);
			op->run(in_vec, out_vec);
		}
		return mem_dense_matrix::create(ret_mat);
#endif
		BOOST_LOG_TRIVIAL(error)
			<< "it doesn't support to apply columns on a tall matrix";
		return dense_matrix::ptr();
	}
	// There are four cases left. In these four cases, apply the function
	// on the rows/columns in the short dimension. We can use mapply to
	// parallelize the computation here.
	else {
		std::vector<detail::matrix_store::const_ptr> ins(1);
		ins[0] = this->get_raw_store();
		size_t out_num_rows;
		size_t out_num_cols;
		if (margin == matrix_margin::MAR_ROW) {
			out_num_rows = this->get_num_rows();
			out_num_cols = op->get_num_out_eles(this->get_num_cols());
		}
		else {
			out_num_rows = op->get_num_out_eles(this->get_num_rows());
			out_num_cols = this->get_num_cols();
		}
		matrix_margin_apply_op::const_ptr apply_op(new matrix_margin_apply_op(
					margin, op, out_num_rows, out_num_cols));
		matrix_layout_t output_layout = (margin == matrix_margin::MAR_ROW
				? matrix_layout_t::L_ROW : matrix_layout_t::L_COL);
		detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
				apply_op, output_layout);
		return dense_matrix::create(ret);
	}
}

////////////////////// Convert the data layout of a matrix ////////////////////

namespace
{

class conv_layout_op: public detail::portion_mapply_op
{
	matrix_layout_t layout;
public:
	conv_layout_op(matrix_layout_t layout, size_t num_rows, size_t num_cols,
			const scalar_type &type): detail::portion_mapply_op(num_rows,
				num_cols, type) {
		this->layout = layout;
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		assert(ins.size() == 1);
		assert(ins[0]->get_global_start_col() == out.get_global_start_col());
		assert(ins[0]->get_global_start_row() == out.get_global_start_row());
		out.copy_from(*ins[0]);
	}

	virtual portion_mapply_op::const_ptr transpose() const {
		matrix_layout_t new_layout;
		if (layout == matrix_layout_t::L_COL)
			new_layout = matrix_layout_t::L_ROW;
		else
			new_layout = matrix_layout_t::L_COL;
		return detail::portion_mapply_op::const_ptr(new conv_layout_op(new_layout,
					get_out_num_cols(), get_out_num_rows(), get_output_type()));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("conv_layout(") + mats[0]->get_name() + ")";
	}
};

}

dense_matrix::ptr dense_matrix::conv2(matrix_layout_t layout) const
{
	if (store_layout() == layout)
		return dense_matrix::create(get_raw_store());

#if 0
	// If the dense matrix has only one row or one column, it's very easy
	// to convert its layout. We don't need to copy data or run computation
	// at all. This only works for in-mem non-virtual matrices. If this is
	// a virtual matrix, getting a row/col from the matrix will trigger
	// materialization. If this is an EM matrix, it'll result in reading
	// the entire column to memory.
	if (get_num_cols() == 1 && !get_data().is_virtual()
			&& get_data().is_in_mem()) {
		detail::vec_store::const_ptr vec = get_data().get_col_vec(0);
		return dense_matrix::create(vec->conv2mat(get_num_rows(),
					get_num_cols(), layout == matrix_layout_t::L_ROW));
	}
	else if (get_num_rows() == 1 && !get_data().is_virtual()
			&& get_data().is_in_mem()) {
		detail::vec_store::const_ptr vec = get_data().get_row_vec(0);
		return dense_matrix::create(vec->conv2mat(get_num_rows(),
					get_num_cols(), layout == matrix_layout_t::L_ROW));
	}
#endif
	if (store->is_sink())
		materialize_self();

	std::vector<detail::matrix_store::const_ptr> ins(1);
	ins[0] = this->get_raw_store();
	conv_layout_op::const_ptr mapply_op(new conv_layout_op(layout,
				get_num_rows(), get_num_cols(), get_type()));
	detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
			mapply_op, layout);
	return dense_matrix::create(ret);
}

dense_matrix::ptr dense_matrix::row_sum() const
{
	bulk_operate::const_ptr add
		= bulk_operate::conv2ptr(get_type().get_basic_ops().get_add());
	return aggregate(matrix_margin::MAR_ROW, agg_operate::create(add));
}

dense_matrix::ptr dense_matrix::col_sum() const
{
	bulk_operate::const_ptr add
		= bulk_operate::conv2ptr(get_type().get_basic_ops().get_add());
	return aggregate(matrix_margin::MAR_COL, agg_operate::create(add));
}

dense_matrix::ptr dense_matrix::row_norm2() const
{
	detail::matrix_stats.inc_multiplies(get_num_rows() * get_num_cols());

	const bulk_uoperate *op = get_type().get_basic_uops().get_op(
			basic_uops::op_idx::SQ);
	dense_matrix::ptr sq_mat = this->sapply(bulk_uoperate::conv2ptr(*op));
	dense_matrix::ptr sums = sq_mat->row_sum();
	op = get_type().get_basic_uops().get_op(basic_uops::op_idx::SQRT);
	return sums->sapply(bulk_uoperate::conv2ptr(*op));
}

dense_matrix::ptr dense_matrix::col_norm2() const
{
	detail::matrix_stats.inc_multiplies(get_num_rows() * get_num_cols());

	const bulk_uoperate *op = get_type().get_basic_uops().get_op(
			basic_uops::op_idx::SQ);
	dense_matrix::ptr sq_mat = this->sapply(bulk_uoperate::conv2ptr(*op));
	dense_matrix::ptr sums = sq_mat->col_sum();
	op = get_type().get_basic_uops().get_op(basic_uops::op_idx::SQRT);
	return sums->sapply(bulk_uoperate::conv2ptr(*op));
}

namespace
{

class copy_op: public detail::portion_mapply_op
{
public:
	copy_op(size_t out_num_rows, size_t out_num_cols,
			const scalar_type &out_type): detail::portion_mapply_op(
			out_num_rows, out_num_cols, out_type) {
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		assert(ins.size() == 1);
#if 0
		// Using bypass-cache memory copy for the right matrix in tall matrix
		// multiplication significantly decreases the performance of tall
		// matrix multiplication.
		bool ret = out.large_copy_from(*ins[0]);
		// Large memory copy may fail, fall back to the default one.
		if (!ret)
#endif
			out.copy_from(*ins[0]);
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		assert(0);
		return detail::portion_mapply_op::const_ptr();
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return std::string();
	}
};

}

detail::matrix_store::const_ptr dense_matrix::_conv_store(bool in_mem,
		int num_nodes) const
{
	// If the current matrix is EM matrix and we want to convert it to
	// an EM matrix, don't do anything.
	if (!in_mem && !store->is_in_mem() && !store->is_virtual())
		return store;
	// If the current matrix is in-mem matrix and it stores in the same
	// number of NUMA nodes as requested, don't do anything.
	if (in_mem && store->is_in_mem() && store->get_num_nodes() == num_nodes
			&& !store->is_virtual())
		return store;

	if (store->is_virtual()) {
		auto ret = detail::virtual_matrix_store::cast(store)->materialize(
				in_mem, num_nodes);
		// Some virtual matrices may not materialize data in the storage
		// we want, we need to convert the storage of the materialized matrix
		// explicitly.
		if (ret->is_in_mem() != in_mem || ret->get_num_nodes() != num_nodes) {
			dense_matrix::ptr tmp = dense_matrix::create(ret);
			tmp = tmp->conv_store(in_mem, num_nodes);
			return tmp->get_raw_store();
		}
		else
			return ret;
	}
	else {
		std::vector<detail::matrix_store::const_ptr> in_mats(1);
		in_mats[0] = store;
		std::vector<detail::matrix_store::ptr> out_mats(1);
		out_mats[0] = detail::matrix_store::create(get_num_rows(), get_num_cols(),
				store_layout(), get_type(), num_nodes, in_mem);
		if (out_mats[0] == NULL)
			return detail::matrix_store::const_ptr();

		detail::portion_mapply_op::const_ptr op(new copy_op(get_num_rows(),
					get_num_cols(), get_type()));
		bool ret = detail::__mapply_portion(in_mats, op, out_mats);
		if (ret)
			return out_mats[0];
		else
			return detail::matrix_store::const_ptr();
	}
}

dense_matrix::ptr dense_matrix::conv_store(bool in_mem, int num_nodes) const
{
	detail::matrix_store::const_ptr store = _conv_store(in_mem, num_nodes);
	if (store)
		return dense_matrix::create(store);
	else
		return dense_matrix::ptr();
}

bool dense_matrix::move_store(bool in_mem, int num_nodes) const
{
	detail::matrix_store::const_ptr store = _conv_store(in_mem, num_nodes);
	if (store == NULL) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't move matrix store to another storage media";
		return false;
	}
	const_cast<dense_matrix *>(this)->store = store;
	return true;
}

bool dense_matrix::drop_cache()
{
	detail::cached_matrix_store::const_ptr cached
		= std::dynamic_pointer_cast<const detail::cached_matrix_store>(store);
	if (cached == NULL)
		return false;
	const_cast<detail::cached_matrix_store &>(*cached).drop_cache();
	this->store = cached->get_underlying();
	return true;
}

size_t dense_matrix::get_num_cached() const
{
	detail::cached_matrix_store::const_ptr cached
		= std::dynamic_pointer_cast<const detail::cached_matrix_store>(store);
	if (cached == NULL)
		return 0;
	return cached->get_num_cached_vecs();
}

dense_matrix::ptr dense_matrix::logic_not() const
{
	if (get_type() != get_scalar_type<bool>()) {
		BOOST_LOG_TRIVIAL(error) << "logic_not only works on boolean matrix";
		return dense_matrix::ptr();
	}

	bulk_uoperate::const_ptr op = bulk_uoperate::conv2ptr(
			*get_type().get_basic_uops().get_op(basic_uops::op_idx::NOT));
	return sapply(op);
}

dense_matrix::ptr dense_matrix::deep_copy() const
{
	std::vector<detail::matrix_store::const_ptr> ins(1);
	ins[0] = store;
	detail::portion_mapply_op::const_ptr op(new copy_op(get_num_rows(),
				get_num_cols(), get_type()));
	return dense_matrix::create(__mapply_portion(ins, op,
				store->store_layout()));
}

namespace
{

class groupby_long_row_mapply_op: public detail::portion_mapply_op
{
	factor_col_vector::const_ptr labels;
	detail::local_matrix_store::const_ptr llabels;
	agg_operate::const_ptr op;
	volatile bool success;
public:
	groupby_long_row_mapply_op(agg_operate::const_ptr op,
			factor_col_vector::const_ptr labels,
			size_t num_cols): detail::portion_mapply_op(
				labels->get_factor().get_num_levels(), num_cols,
				op->get_output_type()) {
		this->op = op;
		this->labels = labels;
		detail::matrix_store::const_ptr tmp = labels->get_raw_store();
		this->llabels = tmp->get_portion(0, 0, labels->get_num_rows(),
				labels->get_num_cols());
		assert(this->llabels);
		success = true;
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		// We can only resize the number of the cols in this operation.
		return local_num_rows == get_out_num_rows();
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const;

	virtual bool is_success() const {
		return success;
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("groupby_row(") + mats[0]->get_name() + ")";
	}
};

class groupby_long_col_mapply_op: public detail::portion_mapply_op
{
	factor_col_vector::const_ptr labels;
	detail::local_matrix_store::const_ptr llabels;
	agg_operate::const_ptr op;
	volatile bool success;
public:
	groupby_long_col_mapply_op(agg_operate::const_ptr op,
			factor_col_vector::const_ptr labels,
			size_t num_rows): detail::portion_mapply_op(
				num_rows, labels->get_factor().get_num_levels(),
				op->get_output_type()) {
		this->op = op;
		this->labels = labels;
		detail::matrix_store::const_ptr tmp = labels->get_raw_store();
		this->llabels = tmp->get_portion(0, 0, labels->get_num_rows(),
				labels->get_num_cols());
		assert(this->llabels);
		// We need to make sure the local matrix is wide.
		if (this->llabels->get_num_cols() == 1 && this->llabels->get_num_rows() > 1)
			this->llabels
				= std::static_pointer_cast<const detail::local_matrix_store>(
						this->llabels->transpose());
		success = true;
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		// We can only resize the number of the rows in this operation.
		return local_num_cols == get_out_num_cols();
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const;

	virtual bool is_success() const {
		return success;
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("groupby_col(") + mats[0]->get_name() + ")";
	}
};

detail::portion_mapply_op::const_ptr groupby_long_row_mapply_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new groupby_long_col_mapply_op(
				op, labels, get_out_num_cols()));
}

void groupby_long_row_mapply_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	// If we have failed on some portions, we don't need to run on
	// the remaining portions.
	if (!success)
		return;

	assert(ins.size() == 1);
	assert(ins[0]->store_layout() == matrix_layout_t::L_ROW);
	assert(out.store_layout() == matrix_layout_t::L_ROW);
	out.reset_data();
	std::vector<bool> agg_flags(out.get_num_rows());
	const detail::local_row_matrix_store &row_in
		= static_cast<const detail::local_row_matrix_store &>(*ins[0]);
	detail::local_row_matrix_store &row_out
		= static_cast<detail::local_row_matrix_store &>(out);
	// We are grouping very long rows. So this matrix has to be a wide matrix.
	bool ret = detail::groupby(*llabels, row_in, *op, matrix_margin::MAR_ROW,
			detail::part_dim_t::PART_DIM2, row_out, agg_flags);
	if (!ret)
		// This operation can only go one direction.
		// so it's fine if multiple threads want to set it concurrently.
		const_cast<groupby_long_row_mapply_op *>(this)->success = false;
}

detail::portion_mapply_op::const_ptr groupby_long_col_mapply_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new groupby_long_row_mapply_op(
				op, labels, get_out_num_rows()));
}

void groupby_long_col_mapply_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	// If we have failed on some portions, we don't need to run on
	// the remaining portions.
	if (!success)
		return;

	assert(ins.size() == 1);
	assert(ins[0]->store_layout() == matrix_layout_t::L_COL);
	assert(out.store_layout() == matrix_layout_t::L_COL);
	out.reset_data();
	std::vector<bool> agg_flags(out.get_num_cols());
	const detail::local_col_matrix_store &col_in
		= static_cast<const detail::local_col_matrix_store &>(*ins[0]);
	detail::local_col_matrix_store &col_out
		= static_cast<detail::local_col_matrix_store &>(out);
	// We are grouping very long cols. So this matrix has to be a tall matrix.
	bool ret = detail::groupby(*llabels, col_in, *op, matrix_margin::MAR_COL,
			detail::part_dim_t::PART_DIM1, col_out, agg_flags);
	if (!ret)
		// This operation can only go one direction.
		// so it's fine if multiple threads want to set it concurrently.
		const_cast<groupby_long_col_mapply_op *>(this)->success = false;
}

}

dense_matrix::ptr dense_matrix::groupby_row(factor_col_vector::const_ptr labels,
		agg_operate::const_ptr op) const
{
	if (labels->get_length() != get_num_rows()) {
		BOOST_LOG_TRIVIAL(error)
			<< "groupby_row: there should be the same #labels as #rows";
		return dense_matrix::ptr();
	}
	if (get_type() != op->get_input_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "groupby_row: the agg op requires diff element types";
		return dense_matrix::ptr();
	}
	if (!op->has_combine()) {
		BOOST_LOG_TRIVIAL(error) << "agg op needs to have combine";
		return dense_matrix::ptr();
	}

	// If input matrices are sink matrices, we materialize them first.
	if (store->is_sink())
		materialize_self();

	if (this->is_wide()) {
		std::vector<detail::matrix_store::const_ptr> mats(1);
		if (store_layout() == matrix_layout_t::L_ROW)
			mats[0] = store;
		else {
			dense_matrix::ptr tmp = conv2(matrix_layout_t::L_ROW);
			mats[0] = tmp->get_raw_store();
		}
		// I need to materialize the vector because this vector is shared by
		// all threads.
		labels->materialize_self();
		detail::portion_mapply_op::const_ptr groupby_op(
				new groupby_long_row_mapply_op(op, labels, get_num_cols()));
		// TODO I need to make this step virtual.
		detail::matrix_store::ptr ret = __mapply_portion_virtual(mats,
				groupby_op, matrix_layout_t::L_ROW);
		if (ret == NULL)
			return dense_matrix::ptr();
		else
			return dense_matrix::create(ret);
	}
	else {
		detail::matrix_store::const_ptr mat;
		if (store_layout() == matrix_layout_t::L_ROW)
			mat = store;
		else {
			dense_matrix::ptr tmp = conv2(matrix_layout_t::L_ROW);
			mat = tmp->get_raw_store();
		}
		return dense_matrix::create(detail::groupby_matrix_store::create(mat,
					labels, matrix_margin::MAR_ROW, op));
	}
}

dense_matrix::ptr dense_matrix::groupby_row(factor_col_vector::const_ptr labels,
		bulk_operate::const_ptr op) const
{
	agg_operate::const_ptr agg = agg_operate::create(op);
	if (agg == NULL)
		return dense_matrix::ptr();
	return groupby_row(labels, agg);
}


class combined_set_operate: public set_operate
{
	// Previous state for the local thread.
	struct prev_state {
		off_t store_idx;
		size_t out_row_idx;
		prev_state() {
			store_idx = -1;
			out_row_idx = 0;
		}
	};
	std::vector<size_t> accu_nrows;
	std::vector<detail::mem_matrix_store::const_ptr> stores;
	// This is per-thread state.
	std::vector<prev_state> prev_states;
public:
	combined_set_operate(
			const std::vector<detail::mem_matrix_store::const_ptr> &stores) {
		this->stores = stores;
		accu_nrows.resize(stores.size() + 1);
		accu_nrows[0] = 0;
		for (size_t i = 0; i < stores.size(); i++)
			accu_nrows[i + 1] = accu_nrows[i] + stores[i]->get_num_rows();
		prev_states.resize(detail::mem_thread_pool::get_global_num_threads());
	}

	virtual void set(void *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const;

	virtual const scalar_type &get_type() const {
		return stores[0]->get_type();
	}
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

void combined_set_operate::set(void *arr, size_t num_eles, off_t row_idx,
		off_t col_idx) const
{
	combined_set_operate *mutable_this = const_cast<combined_set_operate *>(this);
	int thread_id = detail::mem_thread_pool::get_curr_thread_id();
	// If we don't have previous matrix store.
	if (prev_states[thread_id].store_idx < 0
			// If we have jumped.
			|| prev_states[thread_id].out_row_idx + 1 != (size_t) row_idx) {
		// We have to find the right matrix store for this row.
		size_t i;
		for (i = 0; i < stores.size(); i++) {
			if (accu_nrows[i] <= (size_t) row_idx
					&& (size_t) row_idx < accu_nrows[i + 1]) {
				mutable_this->prev_states[thread_id].store_idx = i;
				break;
			}
		}
		// We have to be able to find a store.
		assert(i < stores.size());
	}
	off_t store_idx = prev_states[thread_id].store_idx;
	assert(store_idx >= 0);
	size_t local_row_idx = row_idx - accu_nrows[store_idx];
	const void *src_row = stores[store_idx]->get_row(local_row_idx);
	assert(num_eles == stores[store_idx]->get_num_cols());
	memcpy(arr, src_row, num_eles * stores[store_idx]->get_entry_size());
	mutable_this->prev_states[thread_id].out_row_idx = row_idx;
	if (local_row_idx == stores[store_idx]->get_num_rows() - 1)
		mutable_this->prev_states[thread_id].store_idx = -1;
}

/*
 * If all block matrices have the same block size, we can rbind individual
 * matrices from each block matrix first and create a new block matrix.
 * All individual matrices are tall matrices.
 */
static dense_matrix::ptr rbind_block(const std::vector<block_matrix::ptr> &mats)
{
	block_matrix::ptr mat = mats.front();
	std::vector<detail::combined_matrix_store::const_ptr> combined_stores(
			mats.size());
	for (size_t i = 0; i < combined_stores.size(); i++)
		combined_stores[i]
			= std::static_pointer_cast<const detail::combined_matrix_store>(
					mats[i]->get_raw_store());
	std::vector<detail::matrix_store::const_ptr> indiv_stores(mat->get_num_blocks());
	for (size_t i = 0; i < mat->get_num_blocks(); i++) {
		// First, get all individual matrices at location i.
		std::vector<detail::mem_matrix_store::const_ptr> mem_stores(mats.size());
		size_t nrow = 0;
		for (size_t j = 0; j < mem_stores.size(); j++) {
			dense_matrix::ptr mat = dense_matrix::create(combined_stores[j]->get_mat(i));
			// TODO it's better to copy data for column major layout.
			mat = mat->conv2(matrix_layout_t::L_ROW);
			mat->materialize_self();
			mem_stores[j] = std::dynamic_pointer_cast<const detail::mem_matrix_store>(
						mat->get_raw_store());
			assert(mem_stores[j]);
			nrow += mem_stores[j]->get_num_rows();
		}
		detail::matrix_store::ptr indiv_store = detail::matrix_store::create(
				nrow, mem_stores[0]->get_num_cols(), matrix_layout_t::L_ROW,
				mem_stores[0]->get_type(), matrix_conf.get_num_nodes(),
				mem_stores[0]->is_in_mem());
		if (indiv_store == NULL)
			return dense_matrix::ptr();

		// Copy the data over.
		indiv_store->set_data(combined_set_operate(mem_stores));
		indiv_stores[i] = indiv_store;
	}
	detail::combined_matrix_store::ptr combined
		= detail::combined_matrix_store::create(indiv_stores,
				mats[0]->store_layout());
	return block_matrix::create(combined);
}

static inline bool is_small(const dense_matrix &mat)
{
	return mat.get_num_rows() <= detail::mem_matrix_store::CHUNK_SIZE
		&& mat.get_num_cols() <= detail::mem_matrix_store::CHUNK_SIZE;
}

dense_matrix::ptr dense_matrix::rbind(const std::vector<dense_matrix::ptr> &mats)
{
	if (mats.empty())
		return dense_matrix::ptr();
	if (mats.size() == 1)
		return mats[0];

	std::vector<block_matrix::ptr> block_mats;
	std::vector<detail::matrix_store::const_ptr> stores(mats.size());
	size_t ncol = mats[0]->get_num_cols();
	size_t nrow = 0;
	bool in_mem = mats[0]->is_in_mem();
	bool small = true;
	const scalar_type &type = mats[0]->get_type();
	for (size_t i = 0; i < mats.size(); i++) {
		dense_matrix::ptr mat = mats[i];
		if (mat->get_num_cols() != ncol) {
			BOOST_LOG_TRIVIAL(error)
				<< "can't rbind two matrices with diff number of columns.";
			return dense_matrix::ptr();
		}
		if (mat->get_type() != type) {
			BOOST_LOG_TRIVIAL(error)
				<< "can't bind two matrices with diff element types";
			return dense_matrix::ptr();
		}
		block_matrix::ptr block_mat = std::dynamic_pointer_cast<block_matrix>(mat);
		if (block_mat)
			block_mats.push_back(block_mat);
		stores[i] = mat->get_raw_store();
		// We create a matrix in memory only if all matrices are in memory.
		in_mem = in_mem && mats[i]->is_in_mem();
		small = small && is_small(*mat);
		nrow += mats[i]->get_num_rows();
	}

	detail::matrix_store::ptr combined;
	// If all matrices are small, we can bind them physically.
	if (small) {
		detail::matrix_store::ptr res = detail::matrix_store::create(nrow,
				ncol, mats[0]->store_layout(), type, -1, true);
		if (res == NULL)
			return dense_matrix::ptr();

		off_t row_idx = 0;
		for (size_t i = 0; i < mats.size(); i++) {
			dense_matrix::ptr tmp = mats[i];
			tmp->materialize_self();
			detail::local_matrix_store::ptr res_part = res->get_portion(row_idx,
					0, tmp->get_num_rows(), tmp->get_num_cols());
			detail::local_matrix_store::const_ptr src_part
				= tmp->get_data().get_portion(0);
			res_part->copy_from(*src_part);
			row_idx += tmp->get_num_rows();
		}
		return dense_matrix::create(res);
	}
	else if (ncol > nrow) {
		// If the input matrices are tall and the combined matrix is wide,
		// We need to handle it differently. TODO right now, we don't handle
		// this case.
		if (mats[0]->get_num_cols() <= mats[0]->get_num_rows())
			return dense_matrix::ptr();

		std::vector<detail::matrix_store::const_ptr> indiv_stores;
		for (size_t i = 0; i < stores.size(); i++) {
			// If this is a combined matrix store, we get all its individual
			// matrices.
			detail::combined_matrix_store::const_ptr combined
				= std::dynamic_pointer_cast<const detail::combined_matrix_store>(
						stores[i]);
			if (combined) {
				for (size_t j = 0; j < combined->get_num_mats(); j++)
					indiv_stores.push_back(combined->get_mat(j));
			}
			else
				indiv_stores.push_back(stores[i]);
		}
		combined = detail::combined_matrix_store::create(indiv_stores,
				matrix_layout_t::L_ROW);
	}
	else if (in_mem) {
		// If the input matrices are wide and the combined matrix is tall,
		// We need to handle it differently. TODO right now, we don't handle
		// this case.
		if (mats[0]->get_num_cols() >= mats[0]->get_num_rows())
			return dense_matrix::ptr();

		// If all block matrices have the same block_size.
		bool same_block_size = true;
		for (size_t i = 1; i < block_mats.size(); i++) {
			if (block_mats[i]->get_block_size() != block_mats[0]->get_block_size()) {
				same_block_size = false;
				break;
			}
		}
		// If all matrices are block matrices with the same block size.
		if (block_mats.size() == mats.size() && same_block_size)
			return rbind_block(block_mats);

		std::vector<detail::mem_matrix_store::const_ptr> mem_stores(mats.size());
		for (size_t i = 0; i < stores.size(); i++) {
			dense_matrix::ptr mat = mats[i]->conv2(matrix_layout_t::L_ROW);
			mat->materialize_self();
			mem_stores[i] = std::dynamic_pointer_cast<const detail::mem_matrix_store>(
						mat->get_raw_store());
		}
		combined = detail::matrix_store::create(nrow, ncol,
				matrix_layout_t::L_ROW, type, matrix_conf.get_num_nodes(),
				in_mem);
		if (combined == NULL)
			return dense_matrix::ptr();

		combined->set_data(combined_set_operate(mem_stores));
	}
	else {
		BOOST_LOG_TRIVIAL(error) << "we don't support rbind on EM tall matrices";
		return dense_matrix::ptr();
	}
	return dense_matrix::create(combined);
}

dense_matrix::ptr dense_matrix::cbind(const std::vector<dense_matrix::ptr> &mats)
{
	if (mats.empty())
		return dense_matrix::ptr();
	if (mats.size() == 1)
		return mats[0];

	size_t nrow = mats[0]->get_num_rows();
	const scalar_type &type = mats[0]->get_type();
	for (size_t i = 0; i < mats.size(); i++) {
		dense_matrix::ptr mat = mats[i];
		if (mat->get_num_rows() != nrow) {
			BOOST_LOG_TRIVIAL(error)
				<< "can't rbind two matrices with diff number of rows.";
			return dense_matrix::ptr();
		}
		if (mat->get_type() != type) {
			BOOST_LOG_TRIVIAL(error)
				<< "can't bind two matrices with diff element types";
			return dense_matrix::ptr();
		}
	}
	std::vector<dense_matrix::ptr> tmats(mats.size());
	for (size_t i = 0; i < mats.size(); i++)
		tmats[i] = mats[i]->transpose();
	dense_matrix::ptr combined = rbind(tmats);
	if (combined)
		combined = combined->transpose();
	return combined;
}

dense_matrix::ptr dense_matrix::scale_cols(col_vec::const_ptr vals) const
{
	if (get_type() == vals->get_type()) {
		bulk_operate::const_ptr multiply
			= bulk_operate::conv2ptr(get_type().get_basic_ops().get_multiply());
		// When we scale columns, it's the same as applying the vector to
		// each row.
		return mapply_rows(vals, multiply);
	}
	else {
		auto match_res = match_type(*this, *vals);
		dense_matrix::ptr new_this = match_res.first;
		dense_matrix::ptr new_m = match_res.second;
		bulk_operate::const_ptr multiply = bulk_operate::conv2ptr(
				new_this->get_type().get_basic_ops().get_multiply());
		return new_this->mapply_rows(col_vec::create(new_m), multiply);
	}
}

dense_matrix::ptr dense_matrix::scale_rows(col_vec::const_ptr vals) const
{
	if (get_type() == vals->get_type()) {
		bulk_operate::const_ptr multiply
			= bulk_operate::conv2ptr(get_type().get_basic_ops().get_multiply());
		// When we scale rows, it's the same as applying the vector to
		// each column.
		return mapply_cols(vals, multiply);
	}
	else {
		auto match_res = match_type(*this, *vals);
		dense_matrix::ptr new_this = match_res.first;
		dense_matrix::ptr new_m = match_res.second;
		bulk_operate::const_ptr multiply = bulk_operate::conv2ptr(
				new_this->get_type().get_basic_ops().get_multiply());
		return new_this->mapply_cols(col_vec::create(new_m), multiply);
	}
}

namespace
{
struct comp_pair_second
{
	bool operator()(const std::pair<off_t, off_t> &p1,
			const std::pair<off_t, off_t> &p2) const {
		if (p1.second == p2.second)
			return p1.first < p2.first;
		else
			return p1.second < p2.second;
	}
};

}

dense_matrix::ptr dense_matrix::set_cols(const std::vector<off_t> &idxs,
			dense_matrix::ptr cols) const
{
	if (idxs.size() != cols->get_num_cols()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The number of new columns doesn't match the col index";
		return dense_matrix::ptr();
	}
	if (get_num_rows() != cols->get_num_rows()) {
		BOOST_LOG_TRIVIAL(error)
			<< "#rows in the new matrix doesn't match the one in this matrix";
		return dense_matrix::ptr();
	}
	if (cols->get_type() != get_type()) {
		BOOST_LOG_TRIVIAL(error) << "new cols has a different type";
		return dense_matrix::ptr();
	}
	for (size_t i = 0; i < idxs.size(); i++)
		if (idxs[i] < 0 || (size_t) idxs[i] >= get_num_cols()) {
			BOOST_LOG_TRIVIAL(error) << "The col index is out of range";
			return dense_matrix::ptr();
		}

	if (is_wide()) {
		cols = cols->conv2(matrix_layout_t::L_COL);
		cols = cols->conv_store(true, -1);

		std::shared_ptr<std::vector<off_t> > idx_ptr(new std::vector<off_t>());
		if (std::is_sorted(idxs.begin(), idxs.end()))
			*idx_ptr = idxs;
		else {
			std::vector<std::pair<off_t, off_t> > p(idxs.size());
			for (size_t i = 0; i < p.size(); i++) {
				p[i].first = i;
				p[i].second = idxs[i];
			}
			std::sort(p.begin(), p.end(), comp_pair_second());
			idx_ptr = std::shared_ptr<std::vector<off_t> >(
					new std::vector<off_t>(idxs.size()));
			std::vector<off_t> locs(idxs.size());
			for (size_t i = 0; i < p.size(); i++) {
				locs[i] = p[i].first;
				idx_ptr->at(i) = p[i].second;
			}
			cols = cols->get_cols(locs);
		}
		auto col_store
			= std::dynamic_pointer_cast<const detail::mem_col_matrix_store>(
					cols->get_raw_store());
		detail::portion_mapply_op::const_ptr op(new detail::set_col_mapply_op(
					idx_ptr, col_store, get_num_rows(), get_num_cols(),
					get_type()));

		std::vector<detail::matrix_store::const_ptr> ins(1);
		ins[0] = get_raw_store();
		return dense_matrix::create(detail::mapply_matrix_store::const_ptr(
					new detail::mapply_matrix_store(ins, op, store_layout())));
	}

	detail::matrix_store::const_ptr col_store;
	if (store_layout() == matrix_layout_t::L_COL)
		col_store = store;
	else {
		dense_matrix::ptr tmp = conv2(matrix_layout_t::L_COL);
		col_store = tmp->get_raw_store();
	}
	dense_matrix::ptr col_mat = dense_matrix::create(col_store);
	cols = cols->conv2(matrix_layout_t::L_COL);

	if (std::is_sorted(idxs.begin(), idxs.end())) {
		// We need to find out which ranges of columns need to be set.
		std::vector<std::pair<off_t, off_t> > set_ranges;
		set_ranges.push_back(std::pair<off_t, off_t>(idxs.front(),
					idxs.front() + 1));
		for (size_t i = 1; i < idxs.size(); i++) {
			if (set_ranges.back().second == idxs[i])
				set_ranges.back().second++;
			else
				set_ranges.push_back(std::pair<off_t, off_t>(idxs[i],
							idxs[i] + 1));
		}
		// Collect all sub-matrices.
		std::vector<dense_matrix::ptr> subs;
		if (set_ranges.front().first > 0)
			subs.push_back(col_mat->get_cols(0, set_ranges.front().first, 1));
		off_t curr_col = 0;
		for (size_t i = 0; i < set_ranges.size(); i++) {
			size_t num_cols = set_ranges[i].second - set_ranges[i].first;
			subs.push_back(cols->get_cols(curr_col, curr_col + num_cols, 1));
			curr_col += num_cols;
			if (i < set_ranges.size() - 1)
				subs.push_back(col_mat->get_cols(set_ranges[i].second,
							set_ranges[i + 1].first, 1));
			else if ((size_t) set_ranges[i].second < get_num_cols())
				subs.push_back(col_mat->get_cols(set_ranges[i].second,
							get_num_cols(), 1));
		}
		return dense_matrix::cbind(subs);
	}
	else {
		std::vector<dense_matrix::ptr> col_vec(get_num_cols());
		for (size_t i = 0; i < get_num_cols(); i++)
			col_vec[i] = col_mat->get_col(i);
		for (size_t i = 0; i < idxs.size(); i++)
			col_vec[idxs[i]] = cols->get_col(i);
		return dense_matrix::cbind(col_vec);
	}
}

dense_matrix::ptr dense_matrix::set_cols(size_t start, size_t stop, size_t step,
			dense_matrix::ptr cols) const
{
	std::vector<off_t> idxs;
	for (size_t i = start; i < stop; i += step)
		idxs.push_back(i);
	return set_cols(idxs, cols);
}

dense_matrix::ptr dense_matrix::set_rows(const std::vector<off_t> &idxs,
			dense_matrix::ptr rows) const
{
	if (idxs.size() != rows->get_num_rows()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The number of new rows doesn't match the row index";
		return dense_matrix::ptr();
	}

	dense_matrix::ptr tmp = transpose();
	if (tmp == NULL)
		return dense_matrix::ptr();
	dense_matrix::ptr cols = rows->transpose();
	if (cols == NULL)
		return dense_matrix::ptr();
	tmp = tmp->set_cols(idxs, cols);
	if (tmp == NULL)
		return dense_matrix::ptr();
	return tmp->transpose();
}

dense_matrix::ptr dense_matrix::set_rows(size_t start, size_t stop, size_t step,
			dense_matrix::ptr rows) const
{
	if ((stop - start) / step != rows->get_num_rows()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The number of new rows doesn't match the row index";
		return dense_matrix::ptr();
	}

	dense_matrix::ptr tmp = transpose();
	if (tmp == NULL)
		return dense_matrix::ptr();
	dense_matrix::ptr cols = rows->transpose();
	if (cols == NULL)
		return dense_matrix::ptr();
	tmp = tmp->set_cols(start, stop, step, cols);
	if (tmp == NULL)
		return dense_matrix::ptr();
	return tmp->transpose();
}

dense_matrix::ptr dense_matrix::set_eles(dense_matrix::ptr idx,
		col_vec::ptr vals) const
{
	if (idx->get_num_cols() != 2) {
		BOOST_LOG_TRIVIAL(error)
			<< "The index matrix should have two cols";
		return dense_matrix::ptr();
	}
	if (idx->get_type() != get_scalar_type<off_t>()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The index matrix should have type 'off_t'";
		return dense_matrix::ptr();
	}
	idx = idx->conv2(matrix_layout_t::L_COL);
	if (idx->get_num_rows() != get_num_rows()) {
		BOOST_LOG_TRIVIAL(error)
			<< "We only support replacing one element in every row";
		return dense_matrix::ptr();
	}
	col_vec::ptr row_idx = col_vec::create(idx->get_col(0));
	if (!row_idx->is_seq()) {
		BOOST_LOG_TRIVIAL(error)
			<< "We only support replacing one element in every row";
		return dense_matrix::ptr();
	}

	if (is_wide()) {
		BOOST_LOG_TRIVIAL(error)
			<< "We only support replacing elements on tall matrices";
		return dense_matrix::ptr();
	}

	std::vector<detail::matrix_store::const_ptr> ins(3);
	ins[0] = get_raw_store();
	ins[1] = idx->get_col(1)->get_raw_store();
	ins[2] = vals->get_raw_store();
	detail::portion_mapply_op::const_ptr op(new detail::set_ele_seq_mapply_op(
				true, get_num_rows(), get_num_cols(), get_type()));
	return dense_matrix::create(detail::mapply_matrix_store::const_ptr(
				new detail::mapply_matrix_store(ins, op, store_layout())));
}

dense_matrix::ptr mapply_ele(const std::vector<dense_matrix::const_ptr> &mats,
		detail::portion_mapply_op::const_ptr op, matrix_layout_t out_layout,
		bool par_access)
{
	// TODO I need to optimize this for block matrices.
	return detail::mapply_portion(mats, op, out_layout, par_access);
}

void dense_matrix::print(FILE *f) const
{
	dense_matrix::ptr tmp = conv2(matrix_layout_t::L_ROW);
	tmp = tmp->conv_store(true, -1);
	detail::mem_matrix_store::const_ptr mem_store
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				tmp->get_raw_store());
	assert(mem_store);

	for (size_t i = 0; i < get_num_rows(); i++) {
		std::string str = get_type().conv2str(mem_store->get_row(i),
				get_num_cols(), ",");
		fprintf(f, "%s\n", str.c_str());
	}
}

namespace
{

class ifelse_portion_op: public detail::portion_mapply_op
{
public:
	ifelse_portion_op(size_t out_num_rows, size_t out_num_cols,
			const scalar_type &type): detail::portion_mapply_op(out_num_rows,
				out_num_cols, type) {
	}

	virtual portion_mapply_op::const_ptr transpose() const {
		return portion_mapply_op::const_ptr(new ifelse_portion_op(
					get_out_num_cols(), get_out_num_rows(), get_output_type()));
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 3);
		return std::string("ifelse(") + mats[0]->get_name() + ","
			+ mats[1]->get_name() + "," + mats[2]->get_name() + ")";
	}
};

void ifelse_portion_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 3);
	assert(ins[0]->get_type() == get_scalar_type<bool>());
	assert(ins[1]->get_type() == ins[2]->get_type());
	assert(out.get_type() == ins[2]->get_type());
	assert(ins[0]->store_layout() == ins[1]->store_layout());
	assert(ins[0]->store_layout() == ins[2]->store_layout());
	assert(ins[0]->store_layout() == out.store_layout());

	const ifelse &ie = out.get_type().get_ifelse();
	if (ins[0]->get_raw_arr() && ins[1]->get_raw_arr()
			&& ins[2]->get_raw_arr() && out.get_raw_arr()) {
		const bool *test = reinterpret_cast<const bool *>(
				ins[0]->get_raw_arr());
		size_t num_eles = ins[0]->get_num_rows() * ins[0]->get_num_cols();
		ie.run(test, num_eles, ins[1]->get_raw_arr(), ins[2]->get_raw_arr(),
				out.get_raw_arr());
	}
	else if (ins[0]->store_layout() == matrix_layout_t::L_ROW) {
		detail::local_row_matrix_store::const_ptr row_in0
			= std::static_pointer_cast<const detail::local_row_matrix_store>(
					ins[0]);
		detail::local_row_matrix_store::const_ptr row_in1
			= std::static_pointer_cast<const detail::local_row_matrix_store>(
					ins[1]);
		detail::local_row_matrix_store::const_ptr row_in2
			= std::static_pointer_cast<const detail::local_row_matrix_store>(
					ins[2]);
		detail::local_row_matrix_store &row_out
			= static_cast<detail::local_row_matrix_store &>(out);
		for (size_t i = 0; i < ins[0]->get_num_rows(); i++) {
			const bool *test = reinterpret_cast<const bool *>(
					row_in0->get_row(i));
			ie.run(test, ins[0]->get_num_cols(), row_in1->get_row(i),
					row_in2->get_row(i), row_out.get_row(i));
		}
	}
	else {
		detail::local_col_matrix_store::const_ptr col_in0
			= std::static_pointer_cast<const detail::local_col_matrix_store>(
					ins[0]);
		detail::local_col_matrix_store::const_ptr col_in1
			= std::static_pointer_cast<const detail::local_col_matrix_store>(
					ins[1]);
		detail::local_col_matrix_store::const_ptr col_in2
			= std::static_pointer_cast<const detail::local_col_matrix_store>(
					ins[2]);
		detail::local_col_matrix_store &col_out
			= static_cast<detail::local_col_matrix_store &>(out);
		for (size_t i = 0; i < ins[0]->get_num_cols(); i++) {
			const bool *test = reinterpret_cast<const bool *>(
					col_in0->get_col(i));
			ie.run(test, ins[0]->get_num_rows(), col_in1->get_col(i),
					col_in2->get_col(i), col_out.get_col(i));
		}
	}
}

}

dense_matrix::ptr dense_matrix::ifelse(const dense_matrix &yes,
		const dense_matrix &no) const
{
	if (yes.get_type() != no.get_type()) {
		fprintf(stderr,
				"ifelse doesn't support yes and no of different types\n");
		return dense_matrix::ptr();
	}
	if (get_num_rows() != no.get_num_rows()
			|| get_num_cols() != no.get_num_cols()
			|| get_num_rows() != yes.get_num_rows()
			|| get_num_cols() != yes.get_num_cols()) {
		fprintf(stderr, "the size of test, yes and no has to be the same\n");
		return dense_matrix::ptr();
	}

	dense_matrix::ptr test_mat = cast_ele_type(get_scalar_type<bool>());
	dense_matrix::ptr yes_mat = yes.clone();
	dense_matrix::ptr no_mat = no.clone();

	// We need to make sure all matrices have the same layout.
	if (yes.store_layout() == no.store_layout())
		test_mat = test_mat->conv2(yes.store_layout());
	else if (test_mat->store_layout() == yes.store_layout())
		no_mat = no.conv2(test_mat->store_layout());
	else
		yes_mat = yes.conv2(test_mat->store_layout());

	detail::portion_mapply_op::const_ptr op
		= detail::portion_mapply_op::const_ptr(new ifelse_portion_op(
					test_mat->get_num_rows(), test_mat->get_num_cols(),
					yes.get_type()));

	std::vector<dense_matrix::const_ptr> mats(3);
	mats[0] = test_mat;
	mats[1] = yes_mat;
	mats[2] = no_mat;
	return mapply_ele(mats, op, test_mat->store_layout());
}

namespace
{

class groupby_ele_portion_op: public detail::portion_mapply_op
{
	agg_operate::const_ptr find_next;
	agg_operate::const_ptr agg_op;
	std::vector<generic_hashtable::ptr> tables;
public:
	groupby_ele_portion_op(agg_operate::const_ptr agg_op): detail::portion_mapply_op(
			0, 0, agg_op->get_output_type()) {
		this->find_next = agg_op->get_input_type().get_agg_ops().get_find_next();
		this->agg_op = agg_op;
		size_t nthreads = detail::mem_thread_pool::get_global_num_threads();
		tables.resize(nthreads);
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		return detail::portion_mapply_op::const_ptr();
	}

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return "";
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const;

	generic_hashtable::ptr get_agg() const;
};

void groupby_ele_portion_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins) const
{
	int thread_id = detail::mem_thread_pool::get_curr_thread_id();
	groupby_ele_portion_op *mutable_this = const_cast<groupby_ele_portion_op *>(this);
	if (tables[thread_id] == NULL)
		mutable_this->tables[thread_id] = ins[0]->get_type().create_hashtable(
				agg_op->get_output_type());
	generic_hashtable::ptr ltable = tables[thread_id];
	assert(ltable);

	size_t llength = ins[0]->get_num_rows() * ins[0]->get_num_cols();
	detail::local_matrix_store::ptr lmat;
	if (ins[0]->store_layout() == matrix_layout_t::L_COL)
		lmat = detail::local_matrix_store::ptr(
			new detail::local_buf_col_matrix_store(0, 0, ins[0]->get_num_rows(),
				ins[0]->get_num_cols(), ins[0]->get_type(), -1));
	else
		lmat = detail::local_matrix_store::ptr(
			new detail::local_buf_row_matrix_store(0, 0, ins[0]->get_num_rows(),
				ins[0]->get_num_cols(), ins[0]->get_type(), -1));
	lmat->copy_from(*ins[0]);
	lmat->get_type().get_sorter().serial_sort(lmat->get_raw_arr(), llength,
			false);

	local_vec_store::ptr lkeys(new local_buf_vec_store(0, llength,
				ins[0]->get_type(), -1));
	local_vec_store::ptr laggs(new local_buf_vec_store( 0, llength,
				agg_op->get_output_type(), -1));

	// Start to aggregate on the values.
	const char *arr = lmat->get_raw_arr();
	size_t key_idx = 0;
	while (llength > 0) {
		size_t num_same = 0;
		find_next->runAgg(llength, arr, &num_same);
		assert(num_same <= llength);
		memcpy(lkeys->get(key_idx), arr, lmat->get_entry_size());
		agg_op->runAgg(num_same, arr, laggs->get(key_idx));

		key_idx++;
		arr += num_same * lmat->get_entry_size();
		llength -= num_same;
	}
	ltable->insert(key_idx, lkeys->get_raw_arr(), laggs->get_raw_arr(), *agg_op);
}

generic_hashtable::ptr groupby_ele_portion_op::get_agg() const
{
	generic_hashtable::ptr ret;
	size_t i;
	// Find a local table.
	for (i = 0; i < tables.size(); i++) {
		if (tables[i]) {
			ret = tables[i];
			break;
		}
	}
	// We need to move to the next local table.
	i++;
	// Merge with other local tables if they exist.
	for (; i < tables.size(); i++)
		if (tables[i])
			ret->merge(*tables[i], *agg_op);
	return ret;
}

}

data_frame::ptr dense_matrix::groupby(agg_operate::const_ptr op, bool with_val,
		bool sorted) const
{
	std::vector<detail::matrix_store::const_ptr> stores(1);
	stores[0] = get_raw_store();
	groupby_ele_portion_op *_portion_op = new groupby_ele_portion_op(op);
	detail::portion_mapply_op::const_ptr portion_op(_portion_op);
	detail::__mapply_portion(stores, portion_op, matrix_layout_t::L_COL);
	generic_hashtable::ptr agg_res = _portion_op->get_agg();

	// The key-value pairs got from the hashtable aren't in any order.
	// We need to sort them before returning them.
	data_frame::ptr ret = data_frame::create();
	data_frame::ptr df = agg_res->conv2df();
	if (sorted) {
		vector::ptr keys = vector::create(df->get_vec(0));
		data_frame::ptr sorted_keys = keys->sort_with_index();
		detail::smp_vec_store::const_ptr idx_store
			= std::dynamic_pointer_cast<const detail::smp_vec_store>(
					sorted_keys->get_vec("idx"));
		detail::smp_vec_store::const_ptr agg_store
			= std::dynamic_pointer_cast<const detail::smp_vec_store>(df->get_vec(1));
		if (with_val)
			ret->add_vec("val", sorted_keys->get_vec("val"));
		ret->add_vec("agg", agg_store->get(*idx_store));
	}
	else {
		if (with_val)
			ret->add_vec("val", df->get_vec(0));
		ret->add_vec("agg", df->get_vec(1));
	}
	return ret;
}

}
