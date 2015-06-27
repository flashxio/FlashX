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
#include "NUMA_dense_matrix.h"
#include "EM_dense_matrix.h"
#include "generic_type.h"
#include "rand_gen.h"
#include "one_val_matrix_store.h"
#include "local_matrix_store.h"
#include "virtual_matrix_store.h"
#include "mapply_matrix_store.h"
#include "vector.h"

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
		if (addr == NULL) {
			addr = new rand_gen_wrapper(type.create_rand_gen(min, max));
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

dense_matrix::ptr dense_matrix::_create_rand(const scalar_variable &min,
		const scalar_variable &max, size_t nrow, size_t ncol,
		matrix_layout_t layout, int num_nodes, bool in_mem)
{
	assert(min.get_type() == max.get_type());
	detail::matrix_store::ptr store = detail::matrix_store::create(
			nrow, ncol, layout, min.get_type(), num_nodes, in_mem);
	store->set_data(rand_init(min, max));
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
		op.run(local_stores, *local_res);
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
	EM_matrix_store &to_mat;
	const portion_mapply_op &op;
public:
	mapply_portion_compute(size_t num_required_reads, EM_matrix_store &mat,
			const portion_mapply_op &_op): to_mat(mat), op(_op) {
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
		size_t global_start_row = first_mat.get_global_start_row();
		size_t global_start_col = first_mat.get_global_start_col();
		size_t res_num_rows;
		size_t res_num_cols;
		if (to_mat.is_wide()) {
			res_num_rows = to_mat.get_num_rows();
			res_num_cols = first_mat.get_num_cols();
		}
		else {
			res_num_rows = first_mat.get_num_rows();
			res_num_cols = to_mat.get_num_cols();
		}
		detail::local_matrix_store::ptr local_res;
		if (to_mat.store_layout() == matrix_layout_t::L_ROW)
			local_res = detail::local_matrix_store::ptr(
					new detail::local_buf_row_matrix_store(global_start_row,
						global_start_col, res_num_rows, res_num_cols,
						to_mat.get_type(), -1));
		else
			local_res = detail::local_matrix_store::ptr(
					new detail::local_buf_col_matrix_store(global_start_row,
						global_start_col, res_num_rows, res_num_cols,
						to_mat.get_type(), -1));
		op.run(local_stores, *local_res);
		to_mat.write_portion_async(local_res, global_start_row, global_start_col);
	}
}

void EM_mat_mapply_dispatcher::create_task(off_t global_start, size_t length)
{
	std::vector<detail::local_matrix_store::const_ptr> local_stores(
			mats.size());
	mapply_portion_compute *mapply_compute = new mapply_portion_compute(
			num_EM_mats, *res_mat, *op);
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

matrix_store::ptr __mapply_portion(
		const std::vector<matrix_store::const_ptr> &mats,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout)
{
	assert(mats.size() >= 1);
	size_t num_chunks = mats.front()->get_num_portions();
	std::pair<size_t, size_t> first_size = mats.front()->get_portion_size();
	size_t tot_len;
	size_t portion_size;
	bool in_mem = mats.front()->is_in_mem();
	if (mats.front()->is_wide()) {
		tot_len = mats.front()->get_num_cols();
		portion_size = first_size.second;
		assert(op->get_out_num_cols() == mats.front()->get_num_cols());
		for (size_t i = 1; i < mats.size(); i++) {
			assert(portion_size == mats[i]->get_portion_size().second);
			assert(mats[i]->store_layout() == mats.front()->store_layout());
			assert(mats[i]->get_num_cols() == tot_len);
			in_mem = in_mem && mats[i]->is_in_mem();
		}
	}
	else {
		tot_len = mats.front()->get_num_rows();
		portion_size = first_size.first;
		assert(op->get_out_num_rows() == mats.front()->get_num_rows());
		for (size_t i = 1; i < mats.size(); i++) {
			assert(portion_size == mats[i]->get_portion_size().first);
			assert(mats[i]->store_layout() == mats.front()->store_layout());
			assert(mats[i]->get_num_rows() == tot_len);
			in_mem = in_mem && mats[i]->is_in_mem();
		}
	}

	if (in_mem) {
		int num_nodes = mats[0]->get_num_nodes();
		for (size_t i = 1; i < mats.size(); i++)
			assert(num_nodes == mats[i]->get_num_nodes());
		detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
				op->get_out_num_rows(), op->get_out_num_cols(),
				out_layout, op->get_output_type(), num_nodes);

		std::vector<detail::local_matrix_store::const_ptr> local_stores(
				mats.size());
		detail::mem_thread_pool::ptr mem_threads
			= detail::mem_thread_pool::get_global_mem_threads();
		for (size_t i = 0; i < num_chunks; i++) {
			detail::local_matrix_store::ptr local_res = res->get_portion(i);
			for (size_t j = 0; j < local_stores.size(); j++) {
				local_stores[j] = mats[j]->get_portion(i);
				assert(local_res->get_node_id() == local_stores[j]->get_node_id());
			}

			int node_id = local_res->get_node_id();
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
		detail::EM_matrix_store::ptr res = detail::EM_matrix_store::create(
				op->get_out_num_rows(), op->get_out_num_cols(),
				out_layout, op->get_output_type());
		mem_thread_pool::ptr threads = mem_thread_pool::get_global_mem_threads();
		EM_mat_mapply_dispatcher::ptr dispatcher(
				new EM_mat_mapply_dispatcher(mats, res, op, tot_len,
					portion_size));
		for (size_t i = 0; i < threads->get_num_threads(); i++) {
			io_worker_task *task = new io_worker_task(dispatcher);
			for (size_t j = 0; j < mats.size(); j++) {
				if (!mats[j]->is_in_mem()) {
					const EM_object *obj
						= dynamic_cast<const EM_object *>(mats[j].get());
					task->register_EM_obj(const_cast<EM_object *>(obj));
				}
			}
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

dense_matrix::ptr dense_matrix::inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op,
		matrix_layout_t out_layout) const
{
	if (!verify_inner_prod(m, left_op, right_op))
		return dense_matrix::ptr();

	if (out_layout == matrix_layout_t::L_NONE) {
		if (this->store_layout() == matrix_layout_t::L_ROW)
			out_layout = matrix_layout_t::L_ROW;
		else if (this->is_wide())
			out_layout = matrix_layout_t::L_ROW;
		else
			out_layout = matrix_layout_t::L_COL;
	}

	detail::matrix_store::ptr res = detail::matrix_store::create(
			get_num_rows(), m.get_num_cols(), out_layout,
			right_op.get_output_type(), get_data().get_num_nodes(),
			is_in_mem());
	if (is_wide())
		inner_prod_wide(m.get_data(), left_op, right_op, *res);
	else
		inner_prod_tall(m.get_data(), left_op, right_op, *res);

	return dense_matrix::ptr(new dense_matrix(res));
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

void dense_matrix::inner_prod_tall(const detail::matrix_store &m,
		const bulk_operate &left_op, const bulk_operate &right_op,
		detail::matrix_store &res) const
{
	// We assume the right matrix is small, so we don't need to partition it.
	detail::local_matrix_store::const_ptr local_right = m.get_portion(0);
	assert(local_right->get_num_rows() == m.get_num_rows()
			&& local_right->get_num_cols() == m.get_num_cols());
	// If the left matrix is row-major, the right matrix should be
	// column-major. When the left matrix is tall, the right matrix should
	// be small. It makes sense to convert the right matrix to column major
	// before we break up the left matrix for parallel processing.
	if (!is_wide() && this->store_layout() == matrix_layout_t::L_ROW)
		local_right = local_right->conv2(matrix_layout_t::L_COL);
	const detail::matrix_store &this_store = get_data();
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
	const detail::matrix_store &res;
	std::vector<detail::local_matrix_store::ptr> &local_ms;
public:
	inner_prod_wide_task(detail::local_matrix_store::const_ptr local_store,
			detail::local_matrix_store::const_ptr local_store2,
			const bulk_operate &_left_op, const bulk_operate &_right_op,
			const detail::matrix_store &_res,
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

void dense_matrix::inner_prod_wide(const detail::matrix_store &m,
		const bulk_operate &left_op, const bulk_operate &right_op,
		detail::matrix_store &res) const
{
	assert(this->get_num_rows() == res.get_num_rows());
	assert(m.get_num_cols() == res.get_num_cols());

	const detail::matrix_store &this_store = get_data();
	size_t num_chunks = this_store.get_num_portions();
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	int nthreads = mem_threads->get_num_threads();
	std::vector<detail::local_matrix_store::ptr> local_ms(nthreads);
	for (size_t i = 0; i < num_chunks; i++) {
		detail::local_matrix_store::const_ptr local_store
			= this_store.get_portion(i);
		detail::local_matrix_store::const_ptr local_store2
			= m.get_portion(i);
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
	for (int j = 0; j < nthreads; j++) {
		// It's possible that the local matrix store doesn't exist
		// because the input matrix is very small.
		if (local_ms[j])
			detail::mapply2(*local_res, *local_ms[j], right_op, *local_res);
	}
}

////////////////////////////// Aggregation /////////////////////////////

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

class EM_mat_agg_dispatcher: public detail::EM_portion_dispatcher
{
	detail::matrix_store::const_ptr mat;
	const bulk_operate &op;
	char *res;
public:
	EM_mat_agg_dispatcher(detail::matrix_store::const_ptr mat, char *res,
			const bulk_operate &_op, size_t tot_len,
			size_t portion_size): detail::EM_portion_dispatcher(
				tot_len, portion_size), op(_op) {
		this->mat = mat;
		this->res = res;
	}

	virtual void create_task(off_t global_start, size_t length);
};

class agg_portion_compute: public detail::portion_compute
{
	detail::local_matrix_store::const_ptr local_store;
	const bulk_operate &op;
	char *local_res;
public:
	agg_portion_compute(char *local_res, const bulk_operate &_op): op(_op) {
		this->local_res = local_res;
	}

	void set_buf(detail::local_matrix_store::const_ptr local_store) {
		this->local_store = local_store;
	}

	virtual void run(char *buf, size_t size) {
		detail::aggregate(*local_store, op, local_res);
	}
};

void EM_mat_agg_dispatcher::create_task(off_t global_start, size_t length)
{
	assert(global_start % get_portion_size() == 0);
	off_t res_idx = global_start / get_portion_size();
	char *local_res = res + res_idx * op.output_entry_size();
	agg_portion_compute *agg_compute = new agg_portion_compute(local_res, op);
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

}

scalar_variable::ptr dense_matrix::aggregate(const bulk_operate &op) const
{
	if (!verify_aggregate(op))
		return scalar_variable::ptr();
	scalar_variable::ptr res = op.get_output_type().create_scalar();

	const detail::matrix_store &this_store = get_data();
	size_t num_chunks = this_store.get_num_portions();
	std::unique_ptr<char[]> raw_arr(new char[res->get_size() * num_chunks]);

	detail::mem_thread_pool::ptr threads
		= detail::mem_thread_pool::get_global_mem_threads();
	if (is_in_mem()) {
		for (size_t i = 0; i < num_chunks; i++) {
			detail::local_matrix_store::const_ptr local_store
				= this_store.get_portion(i);

			int node_id = local_store->get_node_id();
			// If the local matrix portion is not assigned to any node, 
			// assign the tasks in round robin fashion.
			if (node_id < 0)
				node_id = i % threads->get_num_nodes();
			threads->process_task(node_id,
					new aggregate_task(local_store, op,
						raw_arr.get() + i * op.output_entry_size()));
		}
	}
	else {
		size_t tot_len;
		size_t portion_size;
		if (is_wide()) {
			tot_len = get_num_cols();
			portion_size = store->get_portion_size().second;
		}
		else {
			tot_len = get_num_rows();
			portion_size = store->get_portion_size().first;
		}
		EM_mat_agg_dispatcher::ptr dispatcher(
				new EM_mat_agg_dispatcher(store, raw_arr.get(), op,
					tot_len, portion_size));
		for (size_t i = 0; i < threads->get_num_threads(); i++) {
			detail::io_worker_task *task = new detail::io_worker_task(dispatcher);
			const detail::EM_object *obj
				= dynamic_cast<const detail::EM_object *>(store.get());
			task->register_EM_obj(const_cast<detail::EM_object *>(obj));
			threads->process_task(i % threads->get_num_nodes(), task);
		}
	}
	threads->wait4complete();

	char raw_res[res->get_size()];
	op.runA(num_chunks, raw_arr.get(), raw_res);
	res->set_raw(raw_res, res->get_size());
	return res;
}

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
};

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
};

}

dense_matrix::ptr dense_matrix::conv2(matrix_layout_t layout) const
{
	std::vector<detail::matrix_store::const_ptr> ins(1);
	ins[0] = this->get_raw_store();
	conv_layout_op::const_ptr mapply_op(new conv_layout_op(layout,
				get_num_rows(), get_num_cols(), get_type()));
	detail::matrix_store::ptr ret = __mapply_portion_virtual(ins,
			mapply_op, layout);
	return dense_matrix::create(ret);
}

}
