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
#include "rand_gen.h"
#include "one_val_matrix_store.h"
#include "local_matrix_store.h"
#include "virtual_matrix_store.h"
#include "mapply_matrix_store.h"
#include "vector.h"
#include "matrix_stats.h"

namespace fm
{

namespace detail
{

matrix_stats_t matrix_stats;

void matrix_stats_t::print_diff(const matrix_stats_t &orig) const
{
	if (this->mem_read_bytes != orig.mem_read_bytes)
		BOOST_LOG_TRIVIAL(info) << "in-mem read "
			<< (this->mem_read_bytes - orig.mem_read_bytes) << " bytes";
	if (this->mem_write_bytes != orig.mem_write_bytes)
		BOOST_LOG_TRIVIAL(info) << "in-mem write "
			<< (this->mem_write_bytes - orig.mem_write_bytes) << " bytes";
	if (this->EM_read_bytes != orig.EM_read_bytes)
		BOOST_LOG_TRIVIAL(info) << "ext-mem read "
			<< (this->EM_read_bytes - orig.EM_read_bytes) << " bytes";
	if (this->EM_write_bytes != orig.EM_write_bytes)
		BOOST_LOG_TRIVIAL(info) << "ext-mem write "
			<< (this->EM_write_bytes - orig.EM_write_bytes) << " bytes";
	if (this->double_multiplies != orig.double_multiplies)
		BOOST_LOG_TRIVIAL(info) << "multiply "
			<< (this->double_multiplies - orig.double_multiplies)
			<< " double float points";
}

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

bool dense_matrix::verify_aggregate(const bulk_operate &op) const
{
	if (op.left_entry_size() != op.right_entry_size()
			|| op.left_entry_size() != op.output_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The input and output type of the operator is different";
		return false;
	}

	if (this->get_entry_size() != op.left_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The matrix entry size is different from the operator";
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

	if (this->store_layout() != m.store_layout()) {
		BOOST_LOG_TRIVIAL(error)
			<< "two matrices in mapply2 don't have the same data layout";
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

bool dense_matrix::verify_apply(apply_margin margin, const arr_apply_operate &op) const
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
};

class sum_agg: public bulk_operate
{
public:
	virtual void runA(size_t num_eles, const void *left_arr1,
			void *output) const {
		const long double *t_input = (const long double *) left_arr1;
		long double *t_output = (long double *) output;
		if (num_eles == 0)
			return;
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

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<long double>();
	}

	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<long double>();
	}

	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<long double>();
	}
};

class double_multiply_operate: public bulk_operate
{
public:
	virtual void runAA(size_t num_eles, const void *left_arr,
			const void *right_arr, void *output_arr) const {
		const double *a = static_cast<const double *>(left_arr);
		const double *b = static_cast<const double *>(right_arr);
		long double *c = static_cast<long double *>(output_arr);
		for (size_t i = 0; i < num_eles; i++)
			c[i] = ((long double) a[i]) * ((long double) b[i]);
	}
	virtual void runAE(size_t num_eles, const void *left_arr,
			const void *right, void *output_arr) const {
		long double a = *static_cast<const double *>(right);
		const double *x = static_cast<const double *>(left_arr);
		long double *c = static_cast<long double *>(output_arr);
		for (size_t i = 0; i < num_eles; i++)
			c[i] = x[i] * a;
	}
	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr, void *output_arr) const {
		long double a = *static_cast<const double *>(left);
		const double *x = static_cast<const double *>(right_arr);
		long double *c = static_cast<long double *>(output_arr);
		for (size_t i = 0; i < num_eles; i++)
			c[i] = x[i] * a;
	}
	virtual void runA(size_t num_eles, const void *left_arr,
			void *output) const {
		assert(0);
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<double>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<double>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<long double>();
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
		scalar_variable::ptr res = sq_mat->aggregate(sum_agg());
		assert(res->get_type() == get_scalar_type<long double>());
		ret = sqrtl(*(long double *) res->get_raw());
	}
	else {
		const bulk_uoperate *op = get_type().get_basic_uops().get_op(
				basic_uops::op_idx::SQ);
		dense_matrix::ptr sq_mat = this->sapply(bulk_uoperate::conv2ptr(*op));
		scalar_variable::ptr res = sq_mat->aggregate(
				sq_mat->get_type().get_basic_ops().get_add());
		res->get_type().get_basic_uops().get_op(
				basic_uops::op_idx::SQRT)->runA(1, res->get_raw(), &ret);
	}
	return ret;
}

namespace
{

template<class T>
class multiply_tall_op: public detail::portion_mapply_op
{
	detail::local_matrix_store::const_ptr Bstore;
	std::vector<detail::local_matrix_store::ptr> Abufs;
	std::vector<detail::local_matrix_store::ptr> res_bufs;
public:
	multiply_tall_op(detail::local_matrix_store::const_ptr Bstore,
			size_t num_threads, size_t out_num_rows,
			size_t out_num_cols): detail::portion_mapply_op(
				out_num_rows, out_num_cols, get_scalar_type<T>()) {
		this->Bstore = Bstore;
		Abufs.resize(num_threads);
		res_bufs.resize(num_threads);
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	virtual detail::portion_mapply_op::const_ptr transpose() const;
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("(") + (mats[0]->get_name()
					+ "*") + Bstore->get_name() + std::string(")");
	}
};

template<class T>
class t_multiply_tall_op: public detail::portion_mapply_op
{
	multiply_tall_op<T> op;
public:
	t_multiply_tall_op(const multiply_tall_op<T> &_op): detail::portion_mapply_op(
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

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		return fm::detail::portion_mapply_op::const_ptr(
				new multiply_tall_op<T>(op));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return op.to_string(mats);
	}
};

template<class T>
detail::portion_mapply_op::const_ptr multiply_tall_op<T>::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new t_multiply_tall_op<T>(*this));
}

template<class T>
void multiply_tall_op<T>::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	detail::local_matrix_store::const_ptr Astore = ins[0];
	detail::matrix_stats.inc_multiplies(
			Astore->get_num_rows() * Astore->get_num_cols() * Bstore->get_num_cols());

	const T *Amat = (const T *) Astore->get_raw_arr();
	// Let's make sure all matrices have the same data layout as the result matrix.
	if (Amat == NULL || Astore->store_layout() != out.store_layout()) {
		detail::pool_task_thread *thread = dynamic_cast<detail::pool_task_thread *>(
				thread::get_curr_thread());
		int thread_id = thread->get_pool_thread_id();
		if (Abufs[thread_id] == NULL
				|| Astore->get_num_rows() != Abufs[thread_id]->get_num_rows()
				|| Astore->get_num_cols() != Abufs[thread_id]->get_num_cols()) {
			if (out.store_layout() == matrix_layout_t::L_COL)
				const_cast<multiply_tall_op<T> *>(this)->Abufs[thread_id]
					= detail::local_matrix_store::ptr(
							new fm::detail::local_buf_col_matrix_store(0, 0,
								Astore->get_num_rows(), Astore->get_num_cols(),
								Astore->get_type(), -1));
			else
				const_cast<multiply_tall_op<T> *>(this)->Abufs[thread_id]
					= detail::local_matrix_store::ptr(
							new fm::detail::local_buf_row_matrix_store(0, 0,
								Astore->get_num_rows(), Astore->get_num_cols(),
								Astore->get_type(), -1));
		}
		Abufs[thread_id]->copy_from(*Astore);
		Amat = (const T *) Abufs[thread_id]->get_raw_arr();
	}
	const T *Bmat = (const T *) Bstore->get_raw_arr();
	assert(Bstore->store_layout() == out.store_layout());
	assert(Amat);
	assert(Bmat);

	T *res_mat = (T *) out.get_raw_arr();
	detail::local_matrix_store::ptr res_buf;
	if (res_mat == NULL) {
		detail::pool_task_thread *thread = dynamic_cast<detail::pool_task_thread *>(
				thread::get_curr_thread());
		int thread_id = thread->get_pool_thread_id();
		if (res_bufs[thread_id] == NULL
				|| out.get_num_rows() != res_bufs[thread_id]->get_num_rows()
				|| out.get_num_cols() != res_bufs[thread_id]->get_num_cols())
			const_cast<multiply_tall_op<T> *>(this)->res_bufs[thread_id]
				= detail::local_matrix_store::ptr(
					new fm::detail::local_buf_col_matrix_store(0, 0,
						out.get_num_rows(), out.get_num_cols(),
						out.get_type(), -1));
		res_buf = res_bufs[thread_id];
		res_mat = (T *) res_buf->get_raw_arr();
	}

	if (out.store_layout() == matrix_layout_t::L_COL)
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
				Astore->get_num_rows(), Bstore->get_num_cols(),
				Astore->get_num_cols(), 1, Amat,
				Astore->get_num_rows(), Bmat, Bstore->get_num_rows(),
				0, res_mat, out.get_num_rows());
	else
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
				Astore->get_num_rows(), Bstore->get_num_cols(),
				Astore->get_num_cols(), 1, Amat,
				Astore->get_num_cols(), Bmat, Bstore->get_num_cols(),
				0, res_mat, out.get_num_cols());
	if (res_buf)
		out.copy_from(*res_buf);
}

template<class T>
class multiply_wide_op: public detail::portion_mapply_op
{
	std::vector<detail::local_matrix_store::ptr> Abufs;
	std::vector<detail::local_matrix_store::ptr> Bbufs;
	std::vector<detail::local_matrix_store::ptr> res_bufs;
	size_t out_num_rows;
	size_t out_num_cols;
	matrix_layout_t Alayout;
	matrix_layout_t Blayout;
public:
	multiply_wide_op(size_t num_threads, size_t out_num_rows, size_t out_num_cols,
			matrix_layout_t required_layout): detail::portion_mapply_op(
				0, 0, get_scalar_type<T>()) {
		Abufs.resize(num_threads);
		Bbufs.resize(num_threads);
		res_bufs.resize(num_threads);
		this->out_num_rows = out_num_rows;
		this->out_num_cols = out_num_cols;
		// We need to transpose the A matrix, so we want the data in the A matrix
		// to be organized in the opposite layout to the required.
		if (required_layout == matrix_layout_t::L_COL)
			Alayout = matrix_layout_t::L_ROW;
		else
			Alayout = matrix_layout_t::L_COL;
		Blayout = required_layout;
	}

	const std::vector<detail::local_matrix_store::ptr> &get_partial_results(
			) const {
		return res_bufs;
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &) const;

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
void multiply_wide_op<T>::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &) const
{
	assert(ins.size() == 2);
	detail::pool_task_thread *thread = dynamic_cast<detail::pool_task_thread *>(
			thread::get_curr_thread());
	int thread_id = thread->get_pool_thread_id();

	detail::local_matrix_store::const_ptr Astore = ins[0];
	const T *Amat = (const T *) Astore->get_raw_arr();
	if (Amat == NULL || Astore->store_layout() != Alayout) {
		if (Abufs[thread_id] == NULL
				|| Astore->get_num_rows() != Abufs[thread_id]->get_num_rows()
				|| Astore->get_num_cols() != Abufs[thread_id]->get_num_cols()) {
			if (Alayout == matrix_layout_t::L_ROW)
				const_cast<multiply_wide_op<T> *>(this)->Abufs[thread_id]
					= detail::local_matrix_store::ptr(
							new fm::detail::local_buf_row_matrix_store(0, 0,
								Astore->get_num_rows(), Astore->get_num_cols(),
								Astore->get_type(), -1));
			else
				const_cast<multiply_wide_op<T> *>(this)->Abufs[thread_id]
					= detail::local_matrix_store::ptr(
							new fm::detail::local_buf_col_matrix_store(0, 0,
								Astore->get_num_rows(), Astore->get_num_cols(),
								Astore->get_type(), -1));
		}
		Abufs[thread_id]->copy_from(*Astore);
		Amat = (const T *) Abufs[thread_id]->get_raw_arr();
	}
	assert(Amat);

	detail::local_matrix_store::const_ptr Bstore = ins[1];
	const T *Bmat = (const T *) Bstore->get_raw_arr();
	if (Bmat == NULL || Bstore->store_layout() != Blayout) {
		if (Bbufs[thread_id] == NULL
				|| Bstore->get_num_rows() != Bbufs[thread_id]->get_num_rows()
				|| Bstore->get_num_cols() != Bbufs[thread_id]->get_num_cols()) {
			if (Blayout == matrix_layout_t::L_COL)
				const_cast<multiply_wide_op<T> *>(this)->Bbufs[thread_id]
					= detail::local_matrix_store::ptr(
							new fm::detail::local_buf_col_matrix_store(0, 0,
								Bstore->get_num_rows(), Bstore->get_num_cols(),
								Bstore->get_type(), -1));
			else
				const_cast<multiply_wide_op<T> *>(this)->Bbufs[thread_id]
					= detail::local_matrix_store::ptr(
							new fm::detail::local_buf_row_matrix_store(0, 0,
								Bstore->get_num_rows(), Bstore->get_num_cols(),
								Bstore->get_type(), -1));
		}
		Bbufs[thread_id]->copy_from(*Bstore);
		Bmat = (const T *) Bbufs[thread_id]->get_raw_arr();
	}
	assert(Bmat);

	if (res_bufs[thread_id] == NULL) {
		if (Blayout == matrix_layout_t::L_COL)
			const_cast<multiply_wide_op<T> *>(this)->res_bufs[thread_id]
				= detail::local_matrix_store::ptr(
						new fm::detail::local_buf_col_matrix_store(0, 0,
							out_num_rows, out_num_cols, get_scalar_type<T>(), -1));
		else
			const_cast<multiply_wide_op<T> *>(this)->res_bufs[thread_id]
				= detail::local_matrix_store::ptr(
						new fm::detail::local_buf_row_matrix_store(0, 0,
							out_num_rows, out_num_cols, get_scalar_type<T>(), -1));
		res_bufs[thread_id]->reset_data();
	}
	assert(res_bufs[thread_id]->store_layout() == Blayout);
	T *res_mat = (T *) res_bufs[thread_id]->get_raw_arr();
	// The A matrix is the transpose of the matrix we need. Since the A matrix
	// is stored in contiguous memory and is organized in row major, we can
	// easily interpret it as its transpose by switching its #rows and #cols.
	if (Blayout == matrix_layout_t::L_COL)
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
				Astore->get_num_cols(), Bstore->get_num_cols(),
				Astore->get_num_rows(), 1, Amat,
				Astore->get_num_cols(), Bmat, Bstore->get_num_rows(),
				1, res_mat, out_num_rows);
	else
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
				Astore->get_num_cols(), Bstore->get_num_cols(),
				Astore->get_num_rows(), 1, Amat,
				Astore->get_num_rows(), Bmat, Bstore->get_num_cols(),
				1, res_mat, out_num_cols);
}

}

static dense_matrix::ptr blas_multiply_tall(const dense_matrix &m1,
		const dense_matrix &m2, matrix_layout_t out_layout)
{
	if (out_layout == matrix_layout_t::L_NONE)
		out_layout = m1.store_layout();

	assert(m1.get_type() == get_scalar_type<double>());
	assert(m2.get_type() == get_scalar_type<double>());
	// We assume the right matrix is small, so we don't need to partition it.
	detail::local_matrix_store::const_ptr local_right
		= m2.get_data().get_portion(0);
	assert(local_right->get_num_rows() == m2.get_num_rows()
			&& local_right->get_num_cols() == m2.get_num_cols());
	// Although all data in local_right is stored in contiguous memory,
	// get_raw_arr() may still return NULL. I have to make sure get_raw_arr()
	// works in any way.
	detail::local_matrix_store::ptr tmp;
	if (out_layout == matrix_layout_t::L_COL)
		tmp = detail::local_matrix_store::ptr(
				new detail::local_buf_col_matrix_store(0, 0,
					local_right->get_num_rows(), local_right->get_num_cols(),
					local_right->get_type(), -1));
	else
		tmp = detail::local_matrix_store::ptr(
				new detail::local_buf_row_matrix_store(0, 0,
					local_right->get_num_rows(), local_right->get_num_cols(),
					local_right->get_type(), -1));
	tmp->copy_from(*local_right);
	local_right = tmp;

	std::vector<detail::matrix_store::const_ptr> ins(1);
	ins[0] = m1.get_raw_store();
	detail::mem_thread_pool::ptr threads
		= detail::mem_thread_pool::get_global_mem_threads();
	multiply_tall_op<double>::const_ptr mapply_op(new multiply_tall_op<double>(
				local_right, threads->get_num_threads(), m1.get_num_rows(),
				m2.get_num_cols()));
	return dense_matrix::create(__mapply_portion_virtual(ins, mapply_op,
				out_layout));
}

static dense_matrix::ptr blas_multiply_wide(const dense_matrix &m1,
		const dense_matrix &m2, matrix_layout_t out_layout)
{
	detail::matrix_stats.inc_multiplies(
			m1.get_num_rows() * m1.get_num_cols() * m2.get_num_cols());

	matrix_layout_t required_layout;
	// If both input matrices have the same data layout, it's easy.
	if (m1.store_layout() == m2.store_layout())
		required_layout = m1.store_layout();
	// If they are different, we convert the smaller matrix.
	else if (m1.get_num_rows() * m1.get_num_cols()
			> m2.get_num_rows() * m2.get_num_cols())
		required_layout = m1.store_layout();
	else
		required_layout = m2.store_layout();
	if (out_layout == matrix_layout_t::L_NONE)
		out_layout = required_layout;

	assert(m1.get_type() == get_scalar_type<double>());
	assert(m2.get_type() == get_scalar_type<double>());

	detail::mem_thread_pool::ptr threads
		= detail::mem_thread_pool::get_global_mem_threads();
	size_t nthreads = threads->get_num_threads();

	std::vector<detail::matrix_store::const_ptr> mats(2);
	mats[0] = m1.get_data().transpose();
	assert(mats[0]);
	mats[1] = m2.get_raw_store();
	assert(mats[1]);
	size_t out_num_rows = m1.get_num_rows();
	size_t out_num_cols = m2.get_num_cols();
	std::shared_ptr<multiply_wide_op<double> > op(new multiply_wide_op<double> (
				nthreads, out_num_rows, out_num_cols, required_layout));
	__mapply_portion(mats, op, required_layout);
	std::vector<detail::local_matrix_store::ptr> local_ms
		= op->get_partial_results();
	assert(local_ms.size() == nthreads);

	// Aggregate the results from omp threads.
	// This matrix is small. We can always keep it in memory.
	detail::local_matrix_store::ptr local_res;
	if (required_layout == matrix_layout_t::L_ROW)
		local_res = detail::local_matrix_store::ptr(
				new detail::local_buf_row_matrix_store(0, 0,
					out_num_rows, out_num_cols, m1.get_type(), -1));
	else
		local_res = detail::local_matrix_store::ptr(
				new detail::local_buf_col_matrix_store(0, 0,
					out_num_rows, out_num_cols, m1.get_type(), -1));
	local_res->reset_data();
	const bulk_operate &add = get_scalar_type<double>().get_basic_ops().get_add();
	for (size_t j = 0; j < local_ms.size(); j++) {
		// It's possible that the local matrix store doesn't exist
		// because the input matrix is very small.
		if (local_ms[j])
			detail::mapply2(*local_res, *local_ms[j], add, *local_res);
	}

	detail::matrix_store::ptr res = detail::matrix_store::create(out_num_rows,
			out_num_cols, out_layout, m1.get_type(), -1, true);
	detail::local_matrix_store::ptr tmp = res->get_portion(0);
	assert(tmp->get_num_rows() == res->get_num_rows()
			&& tmp->get_num_cols() == res->get_num_cols());
	// This works for in-mem matrix. TODO Maybe it's not the best way of
	// copying data to the output matrix.
	tmp->copy_from(*local_res);
	return dense_matrix::create(res);
}

dense_matrix::ptr dense_matrix::multiply(const dense_matrix &mat,
		matrix_layout_t out_layout, bool use_blas) const
{
	if (get_type() == get_scalar_type<double>() && use_blas) {
		if (is_wide())
			return blas_multiply_wide(*this, mat, out_layout);
		else
			return blas_multiply_tall(*this, mat, out_layout);
	}
	else if (get_type() == get_scalar_type<double>()) {
		bulk_operate::const_ptr add = bulk_operate::conv2ptr(
				get_scalar_type<long double>().get_basic_ops().get_add());
		bulk_operate::const_ptr multiply(new double_multiply_operate());
		dense_matrix::ptr res;
		if (is_wide())
			res = inner_prod(mat, multiply, add, out_layout);
		else
			res = inner_prod(mat, multiply, add, out_layout);
		assert(res->get_type() == get_scalar_type<long double>());
		dense_matrix::ptr ret = res->cast_ele_type(get_scalar_type<double>());
		return ret;
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

class multiply_scalar_op: public detail::portion_mapply_op
{
	scalar_variable::const_ptr var;
	const bulk_operate &op;
public:
	multiply_scalar_op(scalar_variable::const_ptr var, size_t out_num_rows,
			size_t out_num_cols): detail::portion_mapply_op(out_num_rows,
				out_num_cols, var->get_type()),
			op(var->get_type().get_basic_ops().get_multiply()) {
		this->var = var;
	}
	void run(const std::vector<std::shared_ptr<const detail::local_matrix_store> > &ins,
			detail::local_matrix_store &out) const;
	detail::portion_mapply_op::const_ptr transpose() const {
		return detail::portion_mapply_op::const_ptr(new multiply_scalar_op(
					var, get_out_num_cols(), get_out_num_rows()));
	}

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("(") + mats[0]->get_name()
			+ " * scalar" + std::string(")");
	}
};

void multiply_scalar_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 1);
	detail::matrix_stats.inc_multiplies(
			ins[0]->get_num_rows() * ins[0]->get_num_cols());

	assert(ins[0]->store_layout() == out.store_layout());
	assert(ins[0]->get_num_rows() == out.get_num_rows());
	assert(ins[0]->get_num_cols() == out.get_num_cols());
	if (out.store_layout() == matrix_layout_t::L_COL) {
		const detail::local_col_matrix_store &col_in
			= static_cast<const detail::local_col_matrix_store &>(*ins[0]);
		detail::local_col_matrix_store &col_out
			= static_cast<detail::local_col_matrix_store &>(out);
		for (size_t i = 0; i < out.get_num_cols(); i++)
			op.runAE(out.get_num_rows(), col_in.get_col(i), var->get_raw(),
					col_out.get_col(i));
	}
	else {
		const detail::local_row_matrix_store &row_in
			= static_cast<const detail::local_row_matrix_store &>(*ins[0]);
		detail::local_row_matrix_store &row_out
			= static_cast<detail::local_row_matrix_store &>(out);
		for (size_t i = 0; i < out.get_num_rows(); i++)
			op.runAE(out.get_num_cols(), row_in.get_row(i), var->get_raw(),
					row_out.get_row(i));
	}
}

}

dense_matrix::ptr dense_matrix::_multiply_scalar(
		scalar_variable::const_ptr var) const
{
	if (get_type() != var->get_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't multiply a scalar of incompatible type";
		return dense_matrix::ptr();
	}

	std::vector<detail::matrix_store::const_ptr> stores(1);
	stores[0] = store;
	detail::portion_mapply_op::const_ptr op(new multiply_scalar_op(
				var, get_num_rows(), get_num_cols()));
	detail::matrix_store::ptr ret = __mapply_portion_virtual(stores, op,
			store_layout());
	return dense_matrix::create(ret);
}

namespace
{

/*
 * This class set elements in a container randomly.
 * set_operate can't change its own state and has to be thread-safe when
 * running on multiple threads. However, random generators aren't
 * thread-safe, so we have to create a random generator for each thread.
 */
class rand_init: public set_operate
{
public:
	enum rand_dist_type {
		NORM,
		UNIF,
		MAX_NUM,
	};
private:
	class rand_gen_wrapper {
		rand_gen::ptr gen;
	public:
		rand_gen_wrapper(rand_gen::ptr gen) {
			this->gen = gen;
		}

		rand_gen &get_gen() {
			return *gen;
		}
	};

	pthread_key_t gen_key;
	const scalar_type &type;
	const scalar_variable &var1;
	const scalar_variable &var2;
	rand_dist_type rand_dist;

	rand_gen &get_rand_gen() const {
		void *addr = pthread_getspecific(gen_key);
		if (addr == NULL) {
			if (rand_dist == rand_dist_type::NORM)
				addr = new rand_gen_wrapper(type.create_randn_gen(var1, var2));
			else if (rand_dist == rand_dist_type::UNIF)
				addr = new rand_gen_wrapper(type.create_randu_gen(var1, var2));
			else
				assert(0);
			int ret = pthread_setspecific(gen_key, addr);
			assert(ret == 0);
		}
		rand_gen_wrapper *wrapper = (rand_gen_wrapper *) addr;
		return wrapper->get_gen();
	}

	static void destroy_rand_gen(void *gen) {
		rand_gen_wrapper *wrapper = (rand_gen_wrapper *) gen;
		delete wrapper;
		printf("destroy rand gen\n");
	}
public:
	rand_init(const scalar_variable &_var1, const scalar_variable &_var2,
			rand_dist_type rand_dist): type(_var1.get_type()), var1(
				_var1), var2(_var2) {
		int ret = pthread_key_create(&gen_key, destroy_rand_gen);
		this->rand_dist = rand_dist;
		assert(ret == 0);
	}

	~rand_init() {
		pthread_key_delete(gen_key);
	}

	virtual void set(void *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		get_rand_gen().gen(arr, num_eles);
	}
	virtual const scalar_type &get_type() const {
		return get_rand_gen().get_type();
	}
};

}

dense_matrix::ptr dense_matrix::_create_randu(const scalar_variable &min,
		const scalar_variable &max, size_t nrow, size_t ncol,
		matrix_layout_t layout, int num_nodes, bool in_mem)
{
	assert(min.get_type() == max.get_type());
	detail::matrix_store::ptr store = detail::matrix_store::create(
			nrow, ncol, layout, min.get_type(), num_nodes, in_mem);
	store->set_data(rand_init(min, max, rand_init::rand_dist_type::UNIF));
	return dense_matrix::ptr(new dense_matrix(store));
}

dense_matrix::ptr dense_matrix::_create_randn(const scalar_variable &mean,
		const scalar_variable &var, size_t nrow, size_t ncol,
		matrix_layout_t layout, int num_nodes, bool in_mem)
{
	assert(mean.get_type() == var.get_type());
	detail::matrix_store::ptr store = detail::matrix_store::create(
			nrow, ncol, layout, mean.get_type(), num_nodes, in_mem);
	store->set_data(rand_init(mean, var, rand_init::rand_dist_type::NORM));
	return dense_matrix::ptr(new dense_matrix(store));
}

dense_matrix::ptr dense_matrix::_create_const(scalar_variable::ptr val,
		size_t nrow, size_t ncol, matrix_layout_t layout, int num_nodes,
		bool in_mem)
{
	detail::matrix_store::ptr store(new detail::one_val_matrix_store(
				val, nrow, ncol, layout, num_nodes));
	return dense_matrix::ptr(new dense_matrix(store));
}

void dense_matrix::materialize_self() const
{
	if (!store->is_virtual())
		return;
	const_cast<dense_matrix *>(this)->store
		= detail::virtual_matrix_store::cast(store)->materialize();
}

dense_matrix::ptr dense_matrix::append_cols(
		const std::vector<dense_matrix::ptr> &mats)
{
	std::vector<detail::matrix_store::const_ptr> stores(mats.size());
	for (size_t i = 0; i < mats.size(); i++)
		stores[i] = mats[i]->store;
	detail::matrix_store::const_ptr ret = store->append_cols(stores);
	return dense_matrix::create(ret);
}

/********************************* mapply ************************************/

namespace detail
{

namespace
{

class local_empty_matrix_store: public local_buf_col_matrix_store
{
public:
	local_empty_matrix_store(): local_buf_col_matrix_store(0, 0, 0, 0,
			get_scalar_type<int>(), -1) {
	}
} empty_store;

class mapply_task: public thread_task
{
	std::vector<detail::local_matrix_store::const_ptr> local_stores;
	detail::local_matrix_store::ptr local_res;
	const portion_mapply_op &op;
public:
	mapply_task(
			const std::vector<detail::local_matrix_store::const_ptr> &local_stores,
			const portion_mapply_op &_op,
			detail::local_matrix_store::ptr local_res): op(_op) {
		this->local_stores = local_stores;
		this->local_res = local_res;
	}

	void run() {
		if (local_res)
			op.run(local_stores, *local_res);
		else
			op.run(local_stores, empty_store);
	}
};

class EM_mat_mapply_dispatcher: public detail::EM_portion_dispatcher
{
	std::vector<matrix_store::const_ptr> mats;
	size_t num_EM_mats;
	EM_matrix_store::ptr res_mat;
	portion_mapply_op::const_ptr op;
public:
	EM_mat_mapply_dispatcher(const std::vector<matrix_store::const_ptr> &mats,
			EM_matrix_store::ptr res_mat, portion_mapply_op::const_ptr op,
			size_t tot_len, size_t portion_size): detail::EM_portion_dispatcher(
				tot_len, portion_size) {
		this->mats = mats;
		this->res_mat = res_mat;
		this->op = op;
		num_EM_mats = 0;
		for (size_t i = 0; i < mats.size(); i++)
			if (!mats[i]->is_in_mem())
				num_EM_mats++;
	}

	virtual void create_task(off_t global_start, size_t length);
};

class mapply_portion_compute: public portion_compute
{
	std::vector<detail::local_matrix_store::const_ptr> local_stores;
	size_t num_required_reads;
	size_t num_reads;
	EM_matrix_store::ptr to_mat;
	const portion_mapply_op &op;
public:
	mapply_portion_compute(size_t num_required_reads, EM_matrix_store::ptr mat,
			const portion_mapply_op &_op): op(_op) {
		this->to_mat = mat;
		this->num_required_reads = num_required_reads;
		this->num_reads = 0;
	}

	void set_buf(
			const std::vector<detail::local_matrix_store::const_ptr> &stores) {
		this->local_stores = stores;
	}

	virtual void run(char *buf, size_t size);
};

void mapply_portion_compute::run(char *buf, size_t size)
{
	assert(!local_stores.empty());
	num_reads++;
	if (num_required_reads == num_reads) {
		const detail::local_matrix_store &first_mat = *local_stores.front();
		if (to_mat) {
			size_t global_start_row = first_mat.get_global_start_row();
			size_t global_start_col = first_mat.get_global_start_col();
			size_t res_num_rows;
			size_t res_num_cols;
			if (to_mat->is_wide()) {
				res_num_rows = to_mat->get_num_rows();
				res_num_cols = first_mat.get_num_cols();
			}
			else {
				res_num_rows = first_mat.get_num_rows();
				res_num_cols = to_mat->get_num_cols();
			}
			detail::local_matrix_store::ptr local_res;
			if (to_mat->store_layout() == matrix_layout_t::L_ROW)
				local_res = detail::local_matrix_store::ptr(
						new detail::local_buf_row_matrix_store(global_start_row,
							global_start_col, res_num_rows, res_num_cols,
							to_mat->get_type(), -1));
			else
				local_res = detail::local_matrix_store::ptr(
						new detail::local_buf_col_matrix_store(global_start_row,
							global_start_col, res_num_rows, res_num_cols,
							to_mat->get_type(), -1));
			op.run(local_stores, *local_res);
			to_mat->write_portion_async(local_res, global_start_row,
					global_start_col);
		}
		else
			op.run(local_stores, empty_store);
	}
}

void EM_mat_mapply_dispatcher::create_task(off_t global_start, size_t length)
{
	std::vector<detail::local_matrix_store::const_ptr> local_stores(
			mats.size());
	mapply_portion_compute *mapply_compute = new mapply_portion_compute(
			num_EM_mats, res_mat, *op);
	mapply_portion_compute::ptr compute(mapply_compute);
	for (size_t j = 0; j < local_stores.size(); j++) {
		size_t global_start_row;
		size_t global_start_col;
		size_t num_rows;
		size_t num_cols;
		if (mats[j]->is_wide()) {
			global_start_row = 0;
			global_start_col = global_start;
			num_rows = mats[j]->get_num_rows();
			num_cols = length;
		}
		else {
			global_start_row = global_start;
			global_start_col = 0;
			num_rows = length;
			num_cols = mats[j]->get_num_cols();
		}
		local_stores[j] = mats[j]->get_portion_async(global_start_row,
				global_start_col, num_rows, num_cols, compute);
	}
	mapply_compute->set_buf(local_stores);
}

}

// This function computes the result of mapply.
// We allow to not output a matrix.
matrix_store::ptr __mapply_portion(
		const std::vector<matrix_store::const_ptr> &mats,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout)
{
	assert(mats.size() >= 1);
	bool in_mem = mats.front()->is_in_mem();
	detail::matrix_stats.inc_write_bytes(
			op->get_out_num_rows() * op->get_out_num_cols()
			* op->get_output_type().get_size(), in_mem);

	size_t num_chunks = mats.front()->get_num_portions();
	std::pair<size_t, size_t> first_size = mats.front()->get_portion_size();
	size_t tot_len;
	size_t portion_size;
	if (mats.front()->is_wide()) {
		tot_len = mats.front()->get_num_cols();
		portion_size = first_size.second;
		if (op->get_out_num_cols() > 0)
			assert(op->get_out_num_cols() == mats.front()->get_num_cols());
		for (size_t i = 1; i < mats.size(); i++) {
			assert(portion_size == mats[i]->get_portion_size().second);
			assert(mats[i]->get_num_cols() == tot_len);
			in_mem = in_mem && mats[i]->is_in_mem();
		}
	}
	else {
		tot_len = mats.front()->get_num_rows();
		portion_size = first_size.first;
		if (op->get_out_num_rows() > 0)
			assert(op->get_out_num_rows() == mats.front()->get_num_rows());
		for (size_t i = 1; i < mats.size(); i++) {
			assert(portion_size == mats[i]->get_portion_size().first);
			assert(mats[i]->get_num_rows() == tot_len);
			in_mem = in_mem && mats[i]->is_in_mem();
		}
	}

	if (in_mem) {
		int num_nodes = mats[0]->get_num_nodes();
		for (size_t i = 1; i < mats.size(); i++)
			assert(num_nodes == mats[i]->get_num_nodes());
		detail::mem_matrix_store::ptr res;
		if (op->get_out_num_rows() > 0 && op->get_out_num_cols() > 0)
			res = detail::mem_matrix_store::create(
					op->get_out_num_rows(), op->get_out_num_cols(),
					out_layout, op->get_output_type(), num_nodes);

		std::vector<detail::local_matrix_store::const_ptr> local_stores(
				mats.size());
		detail::mem_thread_pool::ptr mem_threads
			= detail::mem_thread_pool::get_global_mem_threads();
		for (size_t i = 0; i < num_chunks; i++) {
			detail::local_matrix_store::ptr local_res;
			if (res)
				local_res = res->get_portion(i);
			int node_id = -1;
			for (size_t j = 0; j < local_stores.size(); j++) {
				local_stores[j] = mats[j]->get_portion(i);
				if (node_id < 0)
					node_id = local_stores[j]->get_node_id();
				else
					assert(node_id == local_stores[j]->get_node_id());
			}

			// If the local matrix portion is not assigned to any node, 
			// assign the tasks in round robin fashion.
			if (node_id < 0)
				node_id = i % mem_threads->get_num_nodes();
			mem_threads->process_task(node_id,
					new mapply_task(local_stores, *op, local_res));
		}
		mem_threads->wait4complete();
		return res;
	}
	else {
		detail::EM_matrix_store::ptr res;
		if (op->get_out_num_rows() > 0 && op->get_out_num_cols() > 0)
			res = detail::EM_matrix_store::create(
					op->get_out_num_rows(), op->get_out_num_cols(),
					out_layout, op->get_output_type());
		mem_thread_pool::ptr threads = mem_thread_pool::get_global_mem_threads();
		EM_mat_mapply_dispatcher::ptr dispatcher(
				new EM_mat_mapply_dispatcher(mats, res, op, tot_len,
					portion_size));
		for (size_t i = 0; i < threads->get_num_threads(); i++) {
			io_worker_task *task = new io_worker_task(dispatcher, 1);
			for (size_t j = 0; j < mats.size(); j++) {
				if (!mats[j]->is_in_mem()) {
					const EM_object *obj
						= dynamic_cast<const EM_object *>(mats[j].get());
					task->register_EM_obj(const_cast<EM_object *>(obj));
				}
			}
			if (res)
				task->register_EM_obj(res.get());
			threads->process_task(i % threads->get_num_nodes(), task);
		}
		threads->wait4complete();
		return res;
	}
}

matrix_store::ptr __mapply_portion_virtual(
		const std::vector<matrix_store::const_ptr> &stores,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout)
{
	return matrix_store::ptr(new mapply_matrix_store(stores, op,
				out_layout, op->get_out_num_rows(), op->get_out_num_cols()));
}

dense_matrix::ptr mapply_portion(
		const std::vector<dense_matrix::const_ptr> &mats,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout)
{
	std::vector<matrix_store::const_ptr> stores(mats.size());
	for (size_t i = 0; i < stores.size(); i++)
		stores[i] = mats[i]->get_raw_store();
	matrix_store::const_ptr ret(new mapply_matrix_store(stores, op,
				out_layout, op->get_out_num_rows(), op->get_out_num_cols()));
	return dense_matrix::create(ret);
}

}

/*************************** mapply applications *****************************/

///////////////////////// Scale rows and columns //////////////////////////////

namespace
{

class scale_col_op: public detail::portion_mapply_op
{
	detail::mem_vec_store::const_ptr vals;
public:
	scale_col_op(detail::mem_vec_store::const_ptr vals, size_t out_num_rows,
			size_t out_num_cols, const scalar_type &type): detail::portion_mapply_op(
				out_num_rows, out_num_cols, type) {
		this->vals = vals;
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;
	virtual portion_mapply_op::const_ptr transpose() const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("(") + mats[0]->get_name() + "* vec" + std::string(")");
	}
};

class scale_row_op: public detail::portion_mapply_op
{
	detail::mem_vec_store::const_ptr vals;
public:
	scale_row_op(detail::mem_vec_store::const_ptr vals, size_t out_num_rows,
			size_t out_num_cols, const scalar_type &type): detail::portion_mapply_op(
				out_num_rows, out_num_cols, type) {
		this->vals = vals;
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;
	virtual portion_mapply_op::const_ptr transpose() const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("(") + std::string("vec *")
			+ mats[0]->get_name() + std::string(")");
	}
};

detail::portion_mapply_op::const_ptr scale_col_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new scale_row_op(vals,
				get_out_num_cols(), get_out_num_rows(), get_output_type()));
}

void scale_col_op::run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 1);
	detail::matrix_stats.inc_multiplies(
			ins[0]->get_num_rows() * ins[0]->get_num_cols());

	assert(ins[0]->get_global_start_col() == out.get_global_start_col());
	assert(ins[0]->get_global_start_row() == out.get_global_start_row());
	// This is a tall matrix. We divide the matrix horizontally.
	if (ins[0]->get_num_cols() == get_out_num_cols()) {
		// If we use get_raw_arr, it may not work with NUMA vector.
		const char *arr = vals->get_sub_arr(0, vals->get_length());
		assert(arr);
		local_cref_vec_store lvals(arr, 0, vals->get_length(),
				vals->get_type(), -1);
		detail::scale_cols(*ins[0], lvals, out);
	}
	// We divide the matrix vertically.
	else {
		assert(vals->get_length() == get_out_num_cols());
		off_t global_start = ins[0]->get_global_start_col();
		size_t len = ins[0]->get_num_cols();
		local_vec_store::const_ptr portion = vals->get_portion(global_start,
				len);
		assert(portion);
		detail::scale_cols(*ins[0], *portion, out);
	}
}

void scale_row_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 1);
	detail::matrix_stats.inc_multiplies(
			ins[0]->get_num_rows() * ins[0]->get_num_cols());

	assert(ins[0]->get_global_start_col() == out.get_global_start_col());
	assert(ins[0]->get_global_start_row() == out.get_global_start_row());
	// This is a wide matrix. We divide the matrix vertically.
	if (ins[0]->get_num_rows() == get_out_num_rows()) {
		// If we use get_raw_arr, it may not work with NUMA vector.
		const char *arr = vals->get_sub_arr(0, vals->get_length());
		assert(arr);
		local_cref_vec_store lvals(arr, 0, vals->get_length(),
				vals->get_type(), -1);
		detail::scale_rows(*ins[0], lvals, out);
	}
	// We divide the tall matrix horizontally.
	else {
		assert(vals->get_length() == get_out_num_rows());
		off_t global_start = ins[0]->get_global_start_row();
		size_t len = ins[0]->get_num_rows();
		local_vec_store::const_ptr portion = vals->get_portion(global_start,
				len);
		assert(portion);
		detail::scale_rows(*ins[0], *portion, out);
	}
}

detail::portion_mapply_op::const_ptr scale_row_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new scale_col_op(vals,
				get_out_num_cols(), get_out_num_rows(), get_output_type()));
}

}

dense_matrix::ptr dense_matrix::scale_cols(vector::const_ptr vals) const
{
	if (!vals->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error) << "Can't scale columns with an EM vector";
		return dense_matrix::ptr();
	}

	assert(get_num_cols() == vals->get_length());
	assert(get_type() == vals->get_type());
	std::vector<detail::matrix_store::const_ptr> ins(1);
	ins[0] = this->get_raw_store();
	scale_col_op::const_ptr mapply_op(new scale_col_op(
				detail::mem_vec_store::cast(vals->get_raw_store()),
				get_num_rows(), get_num_cols(), get_type()));
	detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
			mapply_op, this->store_layout());
	return dense_matrix::create(ret);
}

dense_matrix::ptr dense_matrix::scale_rows(vector::const_ptr vals) const
{
	if (!vals->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error) << "Can't scale rows with an EM vector";
		return dense_matrix::ptr();
	}

	assert(get_num_rows() == vals->get_length());
	assert(get_type() == vals->get_type());
	std::vector<detail::matrix_store::const_ptr> ins(1);
	ins[0] = this->get_raw_store();
	scale_row_op::const_ptr mapply_op(new scale_row_op(
				detail::mem_vec_store::cast(vals->get_raw_store()),
				get_num_rows(), get_num_cols(), get_type()));
	detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
			mapply_op, this->store_layout());
	return dense_matrix::create(ret);
}

//////////////////////////// Cast the element types ///////////////////////////

namespace
{

class cast_type_op: public detail::portion_mapply_op
{
	vector::const_ptr vals;
public:
	cast_type_op(size_t out_num_rows, size_t out_num_cols,
			const scalar_type &type): detail::portion_mapply_op(
				out_num_rows, out_num_cols, type) {
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;
	virtual portion_mapply_op::const_ptr transpose() const {
		return portion_mapply_op::const_ptr(new cast_type_op(get_out_num_rows(),
					get_out_num_cols(), get_output_type()));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("cast(") + mats[0]->get_name() + ")";
	}
};

void cast_type_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	const type_cast &cast = ins[0]->get_type().get_type_cast(out.get_type());
	assert(ins[0]->store_layout() == out.store_layout());
	if (ins[0]->get_raw_arr() && out.get_raw_arr())
		cast.cast(ins[0]->get_num_rows() * ins[0]->get_num_cols(),
				ins[0]->get_raw_arr(), out.get_raw_arr());
	else if (ins[0]->store_layout() == matrix_layout_t::L_ROW) {
		const detail::local_row_matrix_store &row_in
			= static_cast<const detail::local_row_matrix_store &>(*ins[0]);
		detail::local_row_matrix_store &row_out
			= static_cast<detail::local_row_matrix_store &>(out);
		for (size_t i = 0; i < row_in.get_num_rows(); i++)
			cast.cast(row_in.get_num_cols(), row_in.get_row(i),
					row_out.get_row(i));
	}
	else {
		const detail::local_col_matrix_store &col_in
			= static_cast<const detail::local_col_matrix_store &>(*ins[0]);
		detail::local_col_matrix_store &col_out
			= static_cast<detail::local_col_matrix_store &>(out);
		for (size_t i = 0; i < col_in.get_num_cols(); i++)
			cast.cast(col_in.get_num_rows(), col_in.get_col(i),
					col_out.get_col(i));
	}
}

}

dense_matrix::ptr dense_matrix::cast_ele_type(const scalar_type &type) const
{
	if (!type_cast::require_cast(get_type(), type))
		return dense_matrix::create(get_raw_store());
	else {
		std::vector<detail::matrix_store::const_ptr> ins(1);
		ins[0] = this->get_raw_store();
		cast_type_op::const_ptr mapply_op(new cast_type_op(get_num_rows(),
					get_num_cols(), type));
		detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
				mapply_op, this->store_layout());
		return dense_matrix::create(ret);
	}
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
		return std::string("op(") + mats[0]->get_name() + ", " + mats[1]->get_name() + ")";
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
	detail::mapply2(*ins[0], *ins[1], *op, out);
}

}

dense_matrix::ptr dense_matrix::mapply2(const dense_matrix &m,
		bulk_operate::const_ptr op) const
{
	// The same shape and the same data layout.
	if (!verify_mapply2(m, *op))
		return dense_matrix::ptr();

	std::vector<detail::matrix_store::const_ptr> ins(2);
	ins[0] = this->get_raw_store();
	ins[1] = m.get_raw_store();
	mapply2_op::const_ptr mapply_op(new mapply2_op(op, get_num_rows(),
				get_num_cols()));
	return dense_matrix::create(__mapply_portion_virtual(ins, mapply_op,
				this->store_layout()));
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
		return std::string("uop(") + mats[0]->get_name() + ")";
	}
};

void sapply_op::run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 1);
	assert(ins[0]->get_global_start_col() == out.get_global_start_col());
	assert(ins[0]->get_global_start_row() == out.get_global_start_row());
	detail::sapply(*ins[0], *op, out);
}

}

dense_matrix::ptr dense_matrix::sapply(bulk_uoperate::const_ptr op) const
{
	std::vector<detail::matrix_store::const_ptr> ins(1);
	ins[0] = this->get_raw_store();
	sapply_op::const_ptr mapply_op(new sapply_op(op, get_num_rows(),
				get_num_cols()));
	detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
			mapply_op, this->store_layout());
	return dense_matrix::create(ret);
}

dense_matrix::dense_matrix(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, int num_nodes, bool in_mem)
{
	store = detail::matrix_store::ptr(new detail::one_val_matrix_store(
				type.create_scalar(), nrow, ncol, layout, num_nodes));
}

dense_matrix::ptr dense_matrix::create(size_t nrow, size_t ncol,
		matrix_layout_t layout, const scalar_type &type, int num_nodes,
		bool in_mem)
{
	// If nothing is specified, it creates a zero matrix.
	detail::matrix_store::ptr store(new detail::one_val_matrix_store(
				type.create_scalar(), nrow, ncol, layout, num_nodes));
	return dense_matrix::ptr(new dense_matrix(store));
}

dense_matrix::ptr dense_matrix::create(size_t nrow, size_t ncol,
		matrix_layout_t layout, const scalar_type &type, const set_operate &op,
		int num_nodes, bool in_mem)
{
	detail::matrix_store::ptr store = detail::matrix_store::create(
			nrow, ncol, layout, type, num_nodes, in_mem);
	store->set_data(op);
	return dense_matrix::ptr(new dense_matrix(store));
}

vector::ptr dense_matrix::get_col(off_t idx) const
{
	detail::vec_store::const_ptr vec = get_data().get_col_vec(idx);
	if (vec)
		return vector::create(vec);
	else
		return vector::ptr();
}

vector::ptr dense_matrix::get_row(off_t idx) const
{
	detail::vec_store::const_ptr vec = get_data().get_row_vec(idx);
	if (vec)
		return vector::create(vec);
	else
		return vector::ptr();
}

dense_matrix::ptr dense_matrix::get_cols(const std::vector<off_t> &idxs) const
{
	if (store_layout() == matrix_layout_t::L_COL) {
		const detail::matrix_store::const_ptr ret = get_data().get_cols(idxs);
		if (ret)
			return dense_matrix::ptr(new dense_matrix(ret));
		else
			return dense_matrix::ptr();
	}
	else
		return dense_matrix::ptr();
}

dense_matrix::ptr dense_matrix::get_rows(const std::vector<off_t> &idxs) const
{
	if (store_layout() == matrix_layout_t::L_ROW) {
		const detail::matrix_store::const_ptr ret = get_data().get_rows(idxs);
		if (ret)
			return dense_matrix::ptr(new dense_matrix(ret));
		else
			return dense_matrix::ptr();
	}
	else
		return dense_matrix::ptr();
}

dense_matrix::ptr dense_matrix::transpose() const
{
	return dense_matrix::ptr(new dense_matrix(get_data().transpose()));
}

////////////////////////////// Inner product //////////////////////////////////

dense_matrix::ptr dense_matrix::inner_prod(const dense_matrix &m,
		bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
		matrix_layout_t out_layout) const
{
	if (!verify_inner_prod(m, *left_op, *right_op))
		return dense_matrix::ptr();

	if (out_layout == matrix_layout_t::L_NONE) {
		if (this->store_layout() == matrix_layout_t::L_ROW)
			out_layout = matrix_layout_t::L_ROW;
		else if (this->is_wide())
			out_layout = matrix_layout_t::L_ROW;
		else
			out_layout = matrix_layout_t::L_COL;
	}

	detail::matrix_store::ptr res;
	if (is_wide())
		res = inner_prod_wide(m, left_op, right_op, out_layout);
	else
		res = inner_prod_tall(m, left_op, right_op, out_layout);

	return dense_matrix::ptr(new dense_matrix(res));
}

namespace
{

class inner_prod_tall_op: public detail::portion_mapply_op
{
	detail::local_matrix_store::const_ptr local_right;
	bulk_operate::const_ptr left_op;
	bulk_operate::const_ptr right_op;
public:
	inner_prod_tall_op(detail::local_matrix_store::const_ptr local_right,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			size_t out_num_rows, size_t out_num_cols): detail::portion_mapply_op(
				out_num_rows, out_num_cols, right_op->get_output_type()) {
		this->left_op = left_op;
		this->right_op = right_op;
		this->local_right = local_right;
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
	detail::inner_prod(*ins[0], *local_right, *left_op, *right_op, out);
}

}

detail::matrix_store::ptr dense_matrix::inner_prod_tall(
		const dense_matrix &m, bulk_operate::const_ptr left_op,
		bulk_operate::const_ptr right_op, matrix_layout_t out_layout) const
{
	// We assume the right matrix is small, so we don't need to partition it.
	detail::local_matrix_store::const_ptr local_right = m.get_data().get_portion(0);
	assert(local_right->get_num_rows() == m.get_num_rows()
			&& local_right->get_num_cols() == m.get_num_cols());
	// If the left matrix is row-major, the right matrix should be
	// column-major. When the left matrix is tall, the right matrix should
	// be small. It makes sense to convert the right matrix to column major
	// before we break up the left matrix for parallel processing.
	if (!is_wide() && this->store_layout() == matrix_layout_t::L_ROW)
		local_right = local_right->conv2(matrix_layout_t::L_COL);

	std::vector<detail::matrix_store::const_ptr> ins(1);
	ins[0] = this->get_raw_store();
	inner_prod_tall_op::const_ptr mapply_op(new inner_prod_tall_op(local_right,
				left_op, right_op, get_num_rows(), m.get_num_cols()));
	return __mapply_portion_virtual(ins, mapply_op, out_layout);
}

namespace
{

class inner_prod_wide_op: public detail::portion_mapply_op
{
	const bulk_operate &left_op;
	const bulk_operate &right_op;
	detail::matrix_store::const_ptr res;
	std::vector<detail::local_matrix_store::ptr> local_ms;
public:
	inner_prod_wide_op(const bulk_operate &_left_op, const bulk_operate &_right_op,
			detail::matrix_store::const_ptr res,
			size_t num_threads): detail::portion_mapply_op( 0, 0,
				res->get_type()), left_op(_left_op), right_op(_right_op) {
		local_ms.resize(num_threads);
		this->res = res;
	}

	const std::vector<detail::local_matrix_store::ptr> &get_partial_results() const {
		return local_ms;
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

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
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	detail::pool_task_thread *curr
		= dynamic_cast<detail::pool_task_thread *>(thread::get_curr_thread());
	int thread_id = curr->get_pool_thread_id();
	detail::local_matrix_store::ptr local_m = local_ms[thread_id];
	if (local_m == NULL) {
		int node_id = curr->get_node_id();
		if (res->store_layout() == matrix_layout_t::L_COL)
			local_m = detail::local_matrix_store::ptr(
					new detail::local_buf_col_matrix_store(0, 0,
						res->get_num_rows(), res->get_num_cols(),
						right_op.get_output_type(), node_id));
		else
			local_m = detail::local_matrix_store::ptr(
					new detail::local_buf_row_matrix_store(0, 0,
						res->get_num_rows(), res->get_num_cols(),
						right_op.get_output_type(), node_id));
		local_m->reset_data();
		assert((size_t) thread_id < local_ms.size());
		const_cast<inner_prod_wide_op *>(this)->local_ms[thread_id] = local_m;
	}
	detail::local_matrix_store::const_ptr store
		= std::static_pointer_cast<const detail::local_matrix_store>(
				ins[0]->transpose());
	detail::inner_prod(*store, *ins[1], left_op, right_op, *local_m);
}

}

detail::matrix_store::ptr dense_matrix::inner_prod_wide(
		const dense_matrix &m, bulk_operate::const_ptr left_op,
		bulk_operate::const_ptr right_op, matrix_layout_t out_layout) const
{
	// This matrix is small. We can always keep it in memory.
	detail::matrix_store::ptr res = detail::matrix_store::create(
			get_num_rows(), m.get_num_cols(), out_layout,
			right_op->get_output_type(), -1, true);

	detail::mem_thread_pool::ptr threads
		= detail::mem_thread_pool::get_global_mem_threads();
	size_t nthreads = threads->get_num_threads();

	std::vector<detail::matrix_store::const_ptr> mats(2);
	mats[0] = get_data().transpose();
	assert(mats[0]);
	mats[1] = m.get_raw_store();
	assert(mats[1]);
	std::shared_ptr<inner_prod_wide_op> op(new inner_prod_wide_op(
				*left_op, *right_op, res, nthreads));
	__mapply_portion(mats, op, out_layout);
	std::vector<detail::local_matrix_store::ptr> local_ms
		= op->get_partial_results();
	assert(local_ms.size() == nthreads);

	// Aggregate the results from omp threads.
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
	return res;
}

////////////////////////////// Aggregation /////////////////////////////

namespace
{

class aggregate_task: public thread_task
{
	detail::local_matrix_store::const_ptr local_store;
	const bulk_operate &op;
	detail::agg_margin margin;
	local_ref_vec_store local_res;
public:
	aggregate_task(detail::local_matrix_store::const_ptr local_store,
			const bulk_operate &_op, detail::agg_margin margin,
			char *_local_res, size_t res_len): op(_op), local_res(_local_res,
				0, res_len, _op.get_output_type(), -1) {
		this->local_store = local_store;
		this->margin = margin;
	}

	void run() {
		detail::aggregate(*local_store, op, margin, local_res);
	}
};

class EM_mat_agg_dispatcher: public detail::EM_portion_dispatcher
{
	detail::matrix_store::const_ptr mat;
	const bulk_operate &op;
	// This is a row-major matrix. Each row contains the partial aggregation
	// result of a portion in the original matrix.
	detail::mem_matrix_store::ptr partial_res;
	detail::agg_margin margin;
public:
	EM_mat_agg_dispatcher(detail::matrix_store::const_ptr mat,
			detail::mem_matrix_store::ptr partial_res, const bulk_operate &_op,
			detail::agg_margin margin, size_t tot_len,
			size_t portion_size): detail::EM_portion_dispatcher(tot_len,
				portion_size), op(_op) {
		this->mat = mat;
		this->partial_res = partial_res;
		this->margin = margin;
	}

	virtual void create_task(off_t global_start, size_t length);
};

class agg_portion_compute: public detail::portion_compute
{
	detail::local_matrix_store::const_ptr local_store;
	const bulk_operate &op;
	detail::agg_margin margin;
	local_ref_vec_store local_res;
public:
	agg_portion_compute(const bulk_operate &_op, detail::agg_margin margin,
			char *_local_res, size_t res_len): op(_op), local_res(_local_res,
				0, res_len, _op.get_output_type(), -1) {
		this->margin = margin;
	}

	void set_buf(detail::local_matrix_store::const_ptr local_store) {
		this->local_store = local_store;
	}

	virtual void run(char *buf, size_t size) {
		detail::aggregate(*local_store, op, margin, local_res);
	}
};

void EM_mat_agg_dispatcher::create_task(off_t global_start, size_t length)
{
	assert(global_start % get_portion_size() == 0);
	off_t res_idx = global_start / get_portion_size();
	char *local_res = partial_res->get_row(res_idx);
	agg_portion_compute *agg_compute = new agg_portion_compute(op, margin,
			local_res, partial_res->get_num_cols());
	agg_portion_compute::ptr compute(agg_compute);
	size_t global_start_row;
	size_t global_start_col;
	size_t num_rows;
	size_t num_cols;
	if (mat->is_wide()) {
		global_start_row = 0;
		global_start_col = global_start;
		num_rows = mat->get_num_rows();
		num_cols = length;
	}
	else {
		global_start_row = global_start;
		global_start_col = 0;
		num_rows = length;
		num_cols = mat->get_num_cols();
	}
	detail::local_matrix_store::const_ptr local_store = mat->get_portion_async(
			global_start_row, global_start_col, num_rows, num_cols, compute);
	agg_compute->set_buf(local_store);
}

/*
 * This allows us to aggregate on the shorter dimension.
 */
class matrix_margin_agg_op: public detail::portion_mapply_op
{
	detail::agg_margin margin;
	const bulk_operate &op;
public:
	matrix_margin_agg_op(detail::agg_margin margin, const bulk_operate &_op,
			size_t out_num_rows, size_t out_num_cols): detail::portion_mapply_op(
				out_num_rows, out_num_cols, _op.get_output_type()), op(_op) {
		this->margin = margin;
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		assert(ins.size() == 1);
		// The output matrix is actually a vector.
		if (out.get_num_rows() == 1) {
			assert(out.store_layout() == matrix_layout_t::L_ROW);
			local_ref_vec_store res(
					static_cast<detail::local_row_matrix_store &>(out).get_row(0),
					0, out.get_num_cols(), out.get_type(), -1);
			aggregate(*ins[0], op, margin, res);
		}
		else {
			assert(out.store_layout() == matrix_layout_t::L_COL);
			assert(out.get_num_cols() == 1);
			local_ref_vec_store res(
					static_cast<detail::local_col_matrix_store &>(out).get_col(0),
					0, out.get_num_rows(), out.get_type(), -1);
			aggregate(*ins[0], op, margin, res);
		}
	}

	virtual portion_mapply_op::const_ptr transpose() const {
		detail::agg_margin new_margin
			= (this->margin == detail::agg_margin::MAR_ROW ?
			detail::agg_margin::MAR_COL : detail::agg_margin::MAR_ROW);
		return portion_mapply_op::const_ptr(new matrix_margin_agg_op(
					new_margin, op, get_out_num_cols(), get_out_num_rows()));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("agg(") + mats[0]->get_name() + ")";
	}
};

}

vector::ptr aggregate(detail::matrix_store::const_ptr store,
		detail::agg_margin margin, const bulk_operate &op)
{
	/*
	 * If we aggregate on the shorter dimension.
	 */
	if ((margin == detail::agg_margin::MAR_ROW && !store->is_wide())
			|| (margin == detail::agg_margin::MAR_COL && store->is_wide())) {
		std::vector<detail::matrix_store::const_ptr> ins(1);
		ins[0] = store;
		size_t out_num_rows;
		size_t out_num_cols;
		if (margin == detail::agg_margin::MAR_ROW) {
			out_num_rows = store->get_num_rows();
			out_num_cols = 1;
		}
		else {
			out_num_rows = 1;
			out_num_cols = store->get_num_cols();
		}
		matrix_margin_agg_op::const_ptr agg_op(new matrix_margin_agg_op(
					margin, op, out_num_rows, out_num_cols));
		matrix_layout_t output_layout = (margin == detail::agg_margin::MAR_ROW
				? matrix_layout_t::L_COL : matrix_layout_t::L_ROW);
		detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
				agg_op, output_layout);
		ret->materialize_self();
		// TODO if the result matrix is in external memory, getting
		// the row/column will read the entire row/column into memory.
		if (ret->get_num_cols() == 1)
			return vector::create(ret->get_col_vec(0));
		else
			return vector::create(ret->get_row_vec(0));
	}

	/*
	 * If we aggregate on the entire matrix or on the longer dimension.
	 */

	size_t num_chunks = store->get_num_portions();
	detail::mem_row_matrix_store::ptr partial_res;
	if (margin == detail::agg_margin::BOTH)
		partial_res = detail::mem_row_matrix_store::create(num_chunks,
				1, op.get_output_type());
	// For the next two cases, I assume the partial result is small enough
	// to be kept in memory.
	else if (margin == detail::agg_margin::MAR_ROW)
		partial_res = detail::mem_row_matrix_store::create(num_chunks,
				store->get_num_rows(), op.get_output_type());
	else if (margin == detail::agg_margin::MAR_COL)
		partial_res = detail::mem_row_matrix_store::create(num_chunks,
				store->get_num_cols(), op.get_output_type());
	else
		// This shouldn't happen.
		assert(0);

	detail::mem_thread_pool::ptr threads
		= detail::mem_thread_pool::get_global_mem_threads();
	if (store->is_in_mem()) {
		for (size_t i = 0; i < num_chunks; i++) {
			detail::local_matrix_store::const_ptr local_store
				= store->get_portion(i);

			int node_id = local_store->get_node_id();
			// If the local matrix portion is not assigned to any node, 
			// assign the tasks in round robin fashion.
			if (node_id < 0)
				node_id = i % threads->get_num_nodes();
			threads->process_task(node_id,
					new aggregate_task(local_store, op, margin,
						partial_res->get_row(i), partial_res->get_num_cols()));
		}
	}
	else {
		size_t tot_len;
		size_t portion_size;
		if (store->is_wide()) {
			tot_len = store->get_num_cols();
			portion_size = store->get_portion_size().second;
		}
		else {
			tot_len = store->get_num_rows();
			portion_size = store->get_portion_size().first;
		}
		EM_mat_agg_dispatcher::ptr dispatcher(
				new EM_mat_agg_dispatcher(store, partial_res, op,
					margin, tot_len, portion_size));
		for (size_t i = 0; i < threads->get_num_threads(); i++) {
			detail::io_worker_task *task = new detail::io_worker_task(dispatcher);
			const detail::EM_object *obj
				= dynamic_cast<const detail::EM_object *>(store.get());
			task->register_EM_obj(const_cast<detail::EM_object *>(obj));
			threads->process_task(i % threads->get_num_nodes(), task);
		}
	}
	threads->wait4complete();

	// The last step is to aggregate the partial results from all portions.
	// It runs in serial. I hope it's not a bottleneck.
	detail::local_matrix_store::const_ptr local_res = partial_res->get_portion(
			0, 0, partial_res->get_num_rows(), partial_res->get_num_cols());
	detail::smp_vec_store::ptr res = detail::smp_vec_store::create(
			partial_res->get_num_cols(), partial_res->get_type());
	local_ref_vec_store local_vec(res->get_raw_arr(), 0, res->get_length(),
			res->get_type(), -1);
	detail::aggregate(*local_res, op, detail::agg_margin::MAR_COL, local_vec);
	return vector::create(res);
}

scalar_variable::ptr dense_matrix::aggregate(const bulk_operate &op) const
{
	if (!verify_aggregate(op))
		return scalar_variable::ptr();
	vector::ptr res_vec = fm::aggregate(store, detail::agg_margin::BOTH, op);
	assert(res_vec->get_length() == 1);
	assert(res_vec->is_in_mem());

	scalar_variable::ptr res = op.get_output_type().create_scalar();
	res->set_raw(dynamic_cast<const detail::mem_vec_store &>(
				res_vec->get_data()).get_raw_arr(), res->get_size());
	return res;
}

namespace
{

class matrix_margin_apply_op: public detail::portion_mapply_op
{
	apply_margin margin;
	arr_apply_operate::const_ptr op;
public:
	matrix_margin_apply_op(apply_margin margin, arr_apply_operate::const_ptr op,
			size_t out_num_rows, size_t out_num_cols): detail::portion_mapply_op(
				out_num_rows, out_num_cols, op->get_output_type()) {
		this->margin = margin;
		this->op = op;
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		detail::apply(margin, *op, *ins[0], out);
	}
	virtual portion_mapply_op::const_ptr transpose() const {
		apply_margin new_margin = this->margin == apply_margin::MAR_ROW ?
			apply_margin::MAR_COL : apply_margin::MAR_ROW;
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

dense_matrix::ptr dense_matrix::apply(apply_margin margin,
		arr_apply_operate::const_ptr op) const
{
	assert(op->get_num_out_eles() > 0);
	// In these two cases, we need to convert the matrix store layout
	// before we can apply the function to the matrix.
	detail::matrix_store::const_ptr this_mat;
	if (is_wide() && store_layout() == matrix_layout_t::L_COL
			&& margin == apply_margin::MAR_ROW) {
		dense_matrix::ptr mat = conv2(matrix_layout_t::L_ROW);
		mat->materialize_self();
		this_mat = mat->get_raw_store();
	}
	else if (!is_wide() && store_layout() == matrix_layout_t::L_ROW
			&& margin == apply_margin::MAR_COL) {
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
			&& margin == apply_margin::MAR_ROW) {
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
			&& margin == apply_margin::MAR_COL) {
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
		if (margin == apply_margin::MAR_ROW) {
			out_num_rows = this->get_num_rows();
			out_num_cols = op->get_num_out_eles();
		}
		else {
			out_num_rows = op->get_num_out_eles();
			out_num_cols = this->get_num_cols();
		}
		matrix_margin_apply_op::const_ptr apply_op(new matrix_margin_apply_op(
					margin, op, out_num_rows, out_num_cols));
		matrix_layout_t output_layout = (margin == apply_margin::MAR_ROW
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

	std::vector<detail::matrix_store::const_ptr> ins(1);
	ins[0] = this->get_raw_store();
	conv_layout_op::const_ptr mapply_op(new conv_layout_op(layout,
				get_num_rows(), get_num_cols(), get_type()));
	detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
			mapply_op, layout);
	return dense_matrix::create(ret);
}

vector::ptr dense_matrix::row_sum() const
{
	return fm::aggregate(store, detail::agg_margin::MAR_ROW,
			get_type().get_basic_ops().get_add());
}

vector::ptr dense_matrix::col_sum() const
{
	return fm::aggregate(store, detail::agg_margin::MAR_COL,
			get_type().get_basic_ops().get_add());
}

vector::ptr dense_matrix::row_norm2() const
{
	detail::matrix_stats.inc_multiplies(get_num_rows() * get_num_cols());

	const bulk_uoperate *op = get_type().get_basic_uops().get_op(
			basic_uops::op_idx::SQ);
	dense_matrix::ptr sq_mat = this->sapply(bulk_uoperate::conv2ptr(*op));
	vector::ptr sums = sq_mat->row_sum();
	op = get_type().get_basic_uops().get_op(basic_uops::op_idx::SQRT);
	dense_matrix::ptr sqrt_mat = sums->conv2mat(sums->get_length(), 1,
			false)->sapply(bulk_uoperate::conv2ptr(*op));
	return sqrt_mat->get_col(0);
}

vector::ptr dense_matrix::col_norm2() const
{
	detail::matrix_stats.inc_multiplies(get_num_rows() * get_num_cols());

	const bulk_uoperate *op = get_type().get_basic_uops().get_op(
			basic_uops::op_idx::SQ);
	dense_matrix::ptr sq_mat = this->sapply(bulk_uoperate::conv2ptr(*op));
	vector::ptr sums = sq_mat->col_sum();
	op = get_type().get_basic_uops().get_op(basic_uops::op_idx::SQRT);
	dense_matrix::ptr sqrt_mat = sums->conv2mat(sums->get_length(), 1,
			false)->sapply(bulk_uoperate::conv2ptr(*op));
	return sqrt_mat->get_col(0);
}

}
