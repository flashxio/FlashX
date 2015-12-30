#ifndef __DENSE_MATRIX_H__
#define __DENSE_MATRIX_H__

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

#include <stdlib.h>

#include <vector>
#include <memory>

#include "safs_file.h"

#include "generic_type.h"
#include "matrix_header.h"
#include "bulk_operate.h"
#include "bulk_operate_ext.h"
#include "matrix_store.h"
#include "mem_matrix_store.h"
#include "virtual_matrix_store.h"
#include "factor.h"

namespace fm
{

class bulk_operate;
class bulk_uoperate;
class set_operate;
class arr_apply_operate;
class vector;

enum matrix_margin
{
	MAR_ROW = 1,
	MAR_COL = 2,
	BOTH,
};

class dense_matrix;

namespace detail
{

class local_matrix_store;

class portion_mapply_op
{
	size_t out_num_rows;
	size_t out_num_cols;
	const scalar_type &type;
public:
	typedef std::shared_ptr<const portion_mapply_op> const_ptr;

	portion_mapply_op(size_t out_num_rows, size_t out_num_cols,
			const scalar_type &_type): type(_type) {
		this->out_num_rows = out_num_rows;
		this->out_num_cols = out_num_cols;
	}

	virtual portion_mapply_op::const_ptr transpose() const = 0;

	/*
	 * There are three versions of performing computation on the portions.
	 * The first version performs computation only on input portions;
	 * the second version performs computation on input portions and
	 * outputs only one matrix;
	 * the third version performs computation on input portions and outputs
	 * multiple matrices.
	 */

	virtual void run(
			const std::vector<std::shared_ptr<const local_matrix_store> > &ins) const;
	virtual void run(
			const std::vector<std::shared_ptr<const local_matrix_store> > &ins,
			local_matrix_store &out) const;
	virtual void run(
			const std::vector<std::shared_ptr<const local_matrix_store> > &ins,
			const std::vector<std::shared_ptr<local_matrix_store> > &outs) const;

	virtual std::string to_string(
			const std::vector<matrix_store::const_ptr> &mats) const = 0;

	/*
	 * Give a hint if this operation is aggregation, so we can optimize
	 * the backend accordingly. When this is an aggregation operation,
	 * the second `run' method has to be implemented.
	 */
	virtual bool is_agg() const {
		return false;
	}

	size_t get_out_num_rows() const {
		return out_num_rows;
	}
	size_t get_out_num_cols() const {
		return out_num_cols;
	}
	const scalar_type &get_output_type() const {
		return type;
	}
};

/*
 * These two functions return a virtual matrix that records the computation.
 */

std::shared_ptr<dense_matrix> mapply_portion(
		const std::vector<std::shared_ptr<const dense_matrix> > &mats,
		// A user can specify the layout of the output dense matrix.
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout,
		bool par_access = true);
matrix_store::ptr __mapply_portion_virtual(
		const std::vector<matrix_store::const_ptr> &store,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout,
		bool par_access = true);

/*
 * These three functions return a materialized matrix.
 * The first version determines the storage of the output matrix automatically.
 * The second version allows users to specify the storage for the output matrix.
 * The third version not only allows users to specify the storage for the output
 * matrix, but also allows multiple output matrices.
 */

matrix_store::ptr __mapply_portion(
		const std::vector<matrix_store::const_ptr> &mats,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout,
		bool par_access = true);
matrix_store::ptr __mapply_portion(
		const std::vector<matrix_store::const_ptr> &mats,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout,
		bool out_in_mem, int out_num_nodes, bool par_access = true);
bool __mapply_portion(
		const std::vector<matrix_store::const_ptr> &mats,
		portion_mapply_op::const_ptr op,
		const std::vector<matrix_store::ptr> &out_mats, bool par_access = true);
}

/*
 * This class represents a dense matrix and is able to perform computation
 * on the matrix. However, this class can't modify the matrix data. The only
 * way to modify the matrix is to have the pointer point to another matrix.
 */
class dense_matrix
{
public:
	typedef std::shared_ptr<dense_matrix> ptr;
	typedef std::shared_ptr<const dense_matrix> const_ptr;
private:
	detail::matrix_store::const_ptr store;

	static ptr _create_randu(const scalar_variable &min, const scalar_variable &max,
			size_t nrow, size_t ncol, matrix_layout_t layout, int num_nodes,
			bool in_mem, safs::safs_file_group::ptr group);
	static ptr _create_randn(const scalar_variable &mean, const scalar_variable &var,
			size_t nrow, size_t ncol, matrix_layout_t layout, int num_nodes,
			bool in_mem, safs::safs_file_group::ptr group);
	static ptr _create_const(scalar_variable::ptr val, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes, bool in_mem,
			safs::safs_file_group::ptr group);

	detail::matrix_store::ptr inner_prod_tall(const dense_matrix &m,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t out_layout) const;
	detail::matrix_store::ptr inner_prod_wide(const dense_matrix &m,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t out_layout) const;

	detail::matrix_store::const_ptr _conv_store(bool in_mem, int num_nodes) const;
protected:
	dense_matrix(detail::matrix_store::const_ptr store) {
		this->store = store;
	}
	bool verify_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const;
	bool verify_mapply2(const dense_matrix &m,
			const bulk_operate &op) const;
	bool verify_apply(matrix_margin margin, const arr_apply_operate &op) const;
public:
	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, int num_nodes = -1, bool in_mem = true,
			safs::safs_file_group::ptr group = NULL);
	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, const set_operate &op, int num_nodes = -1,
			bool in_mem = true, safs::safs_file_group::ptr group = NULL);
	static ptr create(detail::matrix_store::const_ptr store) {
		return dense_matrix::ptr(new dense_matrix(store));
	}

	template<class T>
	static ptr create_randu(T _min, T _max, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes = -1, bool in_mem = true,
			safs::safs_file_group::ptr group = NULL) {
		scalar_variable_impl<T> min(_min);
		scalar_variable_impl<T> max(_max);
		return _create_randu(min, max, nrow, ncol, layout, num_nodes, in_mem,
				group);
	}
	template<class T>
	static ptr create_randn(T _mean, T _var, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes = -1, bool in_mem = true,
			safs::safs_file_group::ptr group = NULL) {
		scalar_variable_impl<T> mean(_mean);
		scalar_variable_impl<T> var(_var);
		return _create_randn(mean, var, nrow, ncol, layout, num_nodes, in_mem,
				group);
	}

	template<class T>
	static ptr create_const(T _val, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes = -1, bool in_mem = true,
			safs::safs_file_group::ptr group = NULL) {
		scalar_variable::ptr val(new scalar_variable_impl<T>(_val));
		return _create_const(val, nrow, ncol, layout, num_nodes, in_mem, group);
	}

	dense_matrix() {
	}
	dense_matrix(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, int num_nodes = -1, bool in_mem = true,
			safs::safs_file_group::ptr group = NULL);

	const detail::matrix_store &get_data() const {
		return *store;
	}

	detail::matrix_store::const_ptr get_raw_store() const {
		return store;
	}

	size_t get_entry_size() const {
		return store->get_entry_size();
	}

	size_t get_num_rows() const {
		return store->get_num_rows();
	}

	size_t get_num_cols() const {
		return store->get_num_cols();
	}

	const scalar_type &get_type() const {
		return store->get_type();
	}

	matrix_layout_t store_layout() const {
		return store->store_layout();
	}

	bool is_in_mem() const {
		return store->is_in_mem();
	}

	bool is_wide() const {
		return store->is_wide();
	}

	bool is_virtual() const {
		return store->is_virtual();
	}

	void materialize_self() const;
	void set_materialize_level(materialize_level level);

	template<class T>
	bool is_type() const {
		return get_type() == get_scalar_type<T>();
	}

	/*
	 * We can't change the matrix data that it points to, but we can change
	 * the pointer in the class so that it can point to another matrix data.
	 */
	void assign(const dense_matrix &mat) {
		store = mat.store;
	}

	std::shared_ptr<vector> get_col(off_t idx) const;
	std::shared_ptr<vector> get_row(off_t idx) const;
	dense_matrix::ptr get_cols(const std::vector<off_t> &idxs) const;
	dense_matrix::ptr get_rows(const std::vector<off_t> &idxs) const;
	/*
	 * Clone the matrix.
	 * The class can't modify the matrix data that it points to, but it
	 * can modify the pointer. If someone changes in the pointer in the cloned
	 * matrix, it doesn't affect the current matrix.
	 */
	dense_matrix::ptr clone() const {
		return ptr(new dense_matrix(get_raw_store()));
	}
	dense_matrix::ptr deep_copy() const;

	dense_matrix::ptr transpose() const;
	/*
	 * This converts the data layout of the dense matrix.
	 * It actually generates a virtual matrix that represents the matrix
	 * with required data layout.
	 */
	dense_matrix::ptr conv2(matrix_layout_t layout) const;
	/*
	 * This method converts the storage media of the matrix.
	 * It can convert an in-memory matrix to an EM matrix, or vice versa.
	 * The output matrix is always materialized.
	 */
	dense_matrix::ptr conv_store(bool in_mem, int num_nodes) const;
	bool move_store(bool in_mem, int num_nodes) const;

	dense_matrix::ptr inner_prod(const dense_matrix &m,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t out_layout = matrix_layout_t::L_NONE) const;
	vector::ptr aggregate(matrix_margin margin, agg_operate::const_ptr op) const;
	std::shared_ptr<scalar_variable> aggregate(agg_operate::const_ptr op) const;
	std::shared_ptr<scalar_variable> aggregate(bulk_operate::const_ptr op) const;

	/*
	 * This operator groups rows based on the labels in the factor vector
	 * and aggregate the elements of each column.
	 * It outputs a dense matrix whose #cols == this->#cols and #rows == #levels.
	 * Each row of the output dense matrix is the aggregation of all rows in
	 * the input dense matrix that have the same factor.
	 */
	dense_matrix::ptr groupby_row(factor_vector::const_ptr labels,
			agg_operate::const_ptr) const;
	dense_matrix::ptr groupby_row(factor_vector::const_ptr labels,
			bulk_operate::const_ptr) const;

	dense_matrix::ptr mapply_cols(std::shared_ptr<const vector> vals,
			bulk_operate::const_ptr op) const;
	dense_matrix::ptr mapply_rows(std::shared_ptr<const vector> vals,
			bulk_operate::const_ptr op) const;
	dense_matrix::ptr mapply2(const dense_matrix &m,
			bulk_operate::const_ptr op) const;
	dense_matrix::ptr sapply(bulk_uoperate::const_ptr op) const;
	dense_matrix::ptr apply(matrix_margin margin,
			arr_apply_operate::const_ptr op) const;
	dense_matrix::ptr apply_scalar(scalar_variable::const_ptr var,
			bulk_operate::const_ptr) const;

	dense_matrix::ptr cast_ele_type(const scalar_type &type) const;

	dense_matrix::ptr scale_cols(std::shared_ptr<const vector> vals) const {
		bulk_operate::const_ptr multiply
			= bulk_operate::conv2ptr(get_type().get_basic_ops().get_multiply());
		// When we scale columns, it's the same as applying the vector to
		// each row.
		return mapply_rows(vals, multiply);
	}
	dense_matrix::ptr scale_rows(std::shared_ptr<const vector> vals) const {
		bulk_operate::const_ptr multiply
			= bulk_operate::conv2ptr(get_type().get_basic_ops().get_multiply());
		// When we scale rows, it's the same as applying the vector to
		// each column.
		return mapply_cols(vals, multiply);
	}
	dense_matrix::ptr multiply(const dense_matrix &mat,
			matrix_layout_t out_layout = matrix_layout_t::L_NONE,
			bool use_blas = false) const;

	dense_matrix::ptr add(const dense_matrix &mat) const {
		const bulk_operate &op = get_type().get_basic_ops().get_add();
		return this->mapply2(mat, bulk_operate::conv2ptr(op));
	}
	dense_matrix::ptr minus(const dense_matrix &mat) const {
		const bulk_operate &op = get_type().get_basic_ops().get_sub();
		return this->mapply2(mat, bulk_operate::conv2ptr(op));
	}
	/*
	 * This performs element-wise multiplication between two matrices.
	 */
	dense_matrix::ptr multiply_ele(const dense_matrix &mat) const {
		const bulk_operate &op = get_type().get_basic_ops().get_multiply();
		return this->mapply2(mat, bulk_operate::conv2ptr(op));
	}
	dense_matrix::ptr div(const dense_matrix &mat) const {
		const bulk_operate &op = get_type().get_basic_ops().get_divide();
		return this->mapply2(mat, bulk_operate::conv2ptr(op));
	}
	dense_matrix::ptr pmax(const dense_matrix &mat) const {
		const bulk_operate &op = *get_type().get_basic_ops().get_op(
				basic_ops::op_idx::MAX);
		return this->mapply2(mat, bulk_operate::conv2ptr(op));
	}

	dense_matrix::ptr abs() const {
		bulk_uoperate::const_ptr op = bulk_uoperate::conv2ptr(
				*get_type().get_basic_uops().get_op(basic_uops::op_idx::ABS));
		return sapply(op);
	}

	dense_matrix::ptr logic_not() const;

	std::shared_ptr<vector> row_sum() const;
	std::shared_ptr<vector> col_sum() const;
	std::shared_ptr<vector> row_norm2() const;
	std::shared_ptr<vector> col_norm2() const;

	std::shared_ptr<scalar_variable> sum() const {
		if (get_type() == get_scalar_type<bool>()) {
			dense_matrix::ptr tmp = cast_ele_type(get_scalar_type<size_t>());
			return tmp->sum();
		}
		else
			return aggregate(bulk_operate::conv2ptr(
						get_type().get_basic_ops().get_add()));
	}

	std::shared_ptr<scalar_variable> max() const {
		return aggregate(bulk_operate::conv2ptr(
					*get_type().get_basic_ops().get_op(basic_ops::op_idx::MAX)));
	}

	template<class T>
	dense_matrix::ptr multiply_scalar(T val) const {
		scalar_variable::ptr var(new scalar_variable_impl<T>(val));
		bulk_operate::const_ptr op = bulk_operate::conv2ptr(
				var->get_type().get_basic_ops().get_multiply());
		return apply_scalar(var, op);
	}

	template<class T>
	dense_matrix::ptr add_scalar(T val) const {
		scalar_variable::ptr var(new scalar_variable_impl<T>(val));
		bulk_operate::const_ptr op = bulk_operate::conv2ptr(
				var->get_type().get_basic_ops().get_add());
		return apply_scalar(var, op);
	}

	template<class T>
	dense_matrix::ptr minus_scalar(T val) const {
		scalar_variable::ptr var(new scalar_variable_impl<T>(val));
		bulk_operate::const_ptr op = bulk_operate::conv2ptr(
				var->get_type().get_basic_ops().get_sub());
		return apply_scalar(var, op);
	}

	template<class T>
	dense_matrix::ptr lt_scalar(T val) const {
		scalar_variable::ptr var(new scalar_variable_impl<T>(val));
		bulk_operate::const_ptr op = bulk_operate::conv2ptr(
				*var->get_type().get_basic_ops().get_op(basic_ops::op_idx::LT));
		return apply_scalar(var, op);
	}

	template<class T>
	dense_matrix::ptr pmax_scalar(T val) const {
		scalar_variable::ptr var(new scalar_variable_impl<T>(val));
		bulk_operate::const_ptr op = bulk_operate::conv2ptr(
				*var->get_type().get_basic_ops().get_op(basic_ops::op_idx::MAX));
		return apply_scalar(var, op);
	}

	double norm2() const;
};

class col_vec: public dense_matrix
{
	col_vec(detail::matrix_store::const_ptr mat): dense_matrix(mat) {
		assert(mat->get_num_cols() == 1);
	}
public:
	template<class T>
	static ptr create_randn(size_t len) {
		dense_matrix::ptr mat = dense_matrix::create_randn<T>(0, 1, len, 1,
				matrix_layout_t::L_COL);
		return ptr(new col_vec(mat->get_raw_store()));
	}
	template<class T>
	static ptr create_randu(size_t len) {
		dense_matrix::ptr mat = dense_matrix::create_randu<T>(0, 1, len, 1,
				matrix_layout_t::L_COL);
		return ptr(new col_vec(mat->get_raw_store()));
	}

	col_vec(): dense_matrix(NULL) {
	}

	col_vec(size_t len, const scalar_type &type): dense_matrix(len, 1,
			matrix_layout_t::L_COL, type) {
	}

	size_t get_length() const {
		return get_num_rows();
	}

	col_vec operator=(const dense_matrix &mat) {
		assert(mat.get_num_cols() == 1);
		assign(mat);
		return *this;
	}
};

template<class T>
dense_matrix operator*(const dense_matrix &m, T val)
{
	dense_matrix::ptr ret = m.multiply_scalar<T>(val);
	assert(ret);
	// TODO I shouldn't materialize immediately.
	ret->materialize_self();
	return *ret;
}

template<class T>
dense_matrix operator*(T val, const dense_matrix &m)
{
	dense_matrix::ptr ret = m.multiply_scalar<T>(val);
	assert(ret);
	// TODO I shouldn't materialize immediately.
	ret->materialize_self();
	return *ret;
}

inline dense_matrix operator*(const dense_matrix &m1, const dense_matrix &m2)
{
	dense_matrix::ptr ret = m1.multiply(m2);
	assert(ret);
	// TODO I shouldn't materialize immediately.
	ret->materialize_self();
	return *ret;
}

inline dense_matrix operator*(const dense_matrix &m1, const col_vec &m2)
{
	dense_matrix::ptr ret = m1.multiply(m2);
	assert(ret);
	// TODO I shouldn't materialize immediately.
	ret->materialize_self();
	return *ret;
}

inline dense_matrix operator+(const dense_matrix &m1, const dense_matrix &m2)
{
	dense_matrix::ptr ret = m1.add(m2);
	assert(ret);
	// TODO I shouldn't materialize immediately.
	ret->materialize_self();
	return *ret;
}

inline dense_matrix operator-(const dense_matrix &m1, const dense_matrix &m2)
{
	dense_matrix::ptr ret = m1.minus(m2);
	assert(ret);
	// TODO I shouldn't materialize immediately.
	ret->materialize_self();
	return *ret;
}

template<class T>
inline T as_scalar(const dense_matrix &m)
{
	assert(m.get_type() == get_scalar_type<T>());
	m.materialize_self();
	assert(m.is_in_mem());
	detail::mem_matrix_store::const_ptr mem_m
		= detail::mem_matrix_store::cast(m.get_raw_store());
	return mem_m->get<T>(0, 0);
}

inline dense_matrix t(const dense_matrix &m)
{
	dense_matrix::ptr ret = m.transpose();
	assert(ret);
	return *ret;
}

}

#endif
