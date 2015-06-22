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

#include "log.h"

#include "dense_matrix.h"
#include "bulk_operate.h"
#include "mem_dense_matrix.h"
#include "NUMA_dense_matrix.h"
#include "EM_dense_matrix.h"
#include "generic_type.h"
#include "rand_gen.h"
#include "one_val_matrix_store.h"
#include "local_matrix_store.h"
#include "virtual_matrix_store.h"

namespace fm
{

void dense_matrix::multiply_scalar_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 1);
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

dense_matrix::ptr dense_matrix::create(detail::matrix_store::const_ptr store)
{
	assert(store->is_in_mem());
	return mem_dense_matrix::create(detail::mem_matrix_store::cast(store));
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

dense_matrix::ptr dense_matrix::create(size_t nrow, size_t ncol,
		const scalar_type &type, matrix_layout_t layout, bool in_mem)
{
	if (in_mem)
		return mem_dense_matrix::create(nrow, ncol, layout, type);
	else {
		return dense_matrix::ptr();
#if 0
		if (layout == matrix_layout_t::L_ROW) {
			fprintf(stderr, "EM row dense matrix isn't supported\n");
			return dense_matrix::ptr();
		}
		else
			return EM_col_dense_matrix::create(nrow, ncol, type);
#endif
	}
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

dense_matrix::ptr dense_matrix::multiply(const dense_matrix &mat,
		matrix_layout_t out_layout) const
{
	if (get_type() == get_scalar_type<double>()) {
		const bulk_operate &add
			= get_scalar_type<long double>().get_basic_ops().get_add();
		dense_matrix::ptr res;
		if (is_wide())
			res = inner_prod(mat, double_multiply_operate(), add);
		else
			res = inner_prod(mat, double_multiply_operate(), add);
		assert(res->get_type() == get_scalar_type<long double>());
		dense_matrix::ptr ret = res->cast_ele_type(get_scalar_type<double>());
		ret->materialize_self();
		return ret;
	}
	else
		return inner_prod(mat, get_type().get_basic_ops().get_multiply(),
				get_type().get_basic_ops().get_add(), out_layout);
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
	const scalar_variable &min;
	const scalar_variable &max;

	rand_gen &get_rand_gen() const {
		void *addr = pthread_getspecific(gen_key);
		if (addr == NULL)
			addr = new rand_gen_wrapper(type.create_rand_gen(min, max));
		rand_gen_wrapper *wrapper = (rand_gen_wrapper *) addr;
		return wrapper->get_gen();
	}

	static void destroy_rand_gen(void *gen) {
		rand_gen_wrapper *wrapper = (rand_gen_wrapper *) gen;
		delete wrapper;
		printf("destroy rand gen\n");
	}
public:
	rand_init(const scalar_variable &_min, const scalar_variable &_max): type(
			_min.get_type()), min(_min), max(_max) {
		int ret = pthread_key_create(&gen_key, destroy_rand_gen);
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

mem_dense_matrix::ptr mem_dense_matrix::_create_rand(const scalar_variable &min,
		const scalar_variable &max, size_t nrow, size_t ncol,
		matrix_layout_t layout, int num_nodes)
{
	assert(min.get_type() == max.get_type());
	detail::mem_matrix_store::ptr store = detail::mem_matrix_store::create(
			nrow, ncol, layout, min.get_type(), num_nodes);
	store->set_data(rand_init(min, max));
	return mem_dense_matrix::ptr(new mem_dense_matrix(store));
}

mem_dense_matrix::ptr mem_dense_matrix::_create_const(scalar_variable::ptr val,
		size_t nrow, size_t ncol, matrix_layout_t layout, int num_nodes)
{
	detail::mem_matrix_store::ptr store(new detail::one_val_matrix_store(
				val, nrow, ncol, layout, num_nodes));
	return mem_dense_matrix::ptr(new mem_dense_matrix(store));
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
	return mem_dense_matrix::create(detail::mem_matrix_store::cast(ret));
}

}
