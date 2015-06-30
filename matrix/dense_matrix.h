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

#include "generic_type.h"
#include "matrix_header.h"
#include "bulk_operate.h"
#include "matrix_store.h"
#include "mem_matrix_store.h"

namespace fm
{

class bulk_operate;
class bulk_uoperate;
class set_operate;
class arr_apply_operate;
class vector;

enum apply_margin
{
	MAR_ROW = 1,
	MAR_COL = 2,
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

	virtual void run(
			const std::vector<std::shared_ptr<const local_matrix_store> > &ins,
			local_matrix_store &out) const = 0;

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

std::shared_ptr<dense_matrix> mapply_portion(
		const std::vector<std::shared_ptr<const dense_matrix> > &mats,
		// A user can specify the layout of the output dense matrix.
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout);

matrix_store::ptr __mapply_portion(
		const std::vector<matrix_store::const_ptr> &mats,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout);

matrix_store::ptr __mapply_portion_virtual(
		const std::vector<matrix_store::const_ptr> &store,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout);

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
			bool in_mem);
	static ptr _create_randn(const scalar_variable &mean, const scalar_variable &var,
			size_t nrow, size_t ncol, matrix_layout_t layout, int num_nodes,
			bool in_mem);
	static ptr _create_const(scalar_variable::ptr val, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes, bool in_mem);

	detail::matrix_store::ptr inner_prod_tall(const dense_matrix &m,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t out_layout) const;
	detail::matrix_store::ptr inner_prod_wide(const dense_matrix &m,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t out_layout) const;
	dense_matrix::ptr _multiply_scalar(scalar_variable::const_ptr var) const;
protected:
	dense_matrix(detail::matrix_store::const_ptr store) {
		this->store = store;
	}
	virtual bool verify_aggregate(const bulk_operate &op) const;
	virtual bool verify_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const;
	virtual bool verify_mapply2(const dense_matrix &m,
			const bulk_operate &op) const;
	virtual bool verify_apply(apply_margin margin, const arr_apply_operate &op) const;
public:
	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, int num_nodes = -1, bool in_mem = true);
	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, const set_operate &op, int num_nodes = -1,
			bool in_mem = true);
	static ptr create(detail::matrix_store::const_ptr store) {
		return dense_matrix::ptr(new dense_matrix(store));
	}

	template<class T>
	static ptr create_randu(T _min, T _max, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes = -1, bool in_mem = true) {
		scalar_variable_impl<T> min(_min);
		scalar_variable_impl<T> max(_max);
		return _create_randu(min, max, nrow, ncol, layout, num_nodes, in_mem);
	}
	template<class T>
	static ptr create_randn(T _mean, T _var, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes = -1, bool in_mem = true) {
		scalar_variable_impl<T> mean(_mean);
		scalar_variable_impl<T> var(_var);
		return _create_randn(mean, var, nrow, ncol, layout, num_nodes, in_mem);
	}

	template<class T>
	static ptr create_const(T _val, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes = -1, bool in_mem = true) {
		scalar_variable::ptr val(new scalar_variable_impl<T>(_val));
		return _create_const(val, nrow, ncol, layout, num_nodes, in_mem);
	}

	dense_matrix() {
	}

	virtual ~dense_matrix() {
	}

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
	virtual dense_matrix::ptr clone() const {
		return ptr(new dense_matrix(get_raw_store()));
	}

	dense_matrix::ptr transpose() const;
	dense_matrix::ptr conv2(matrix_layout_t layout) const;

	dense_matrix::ptr inner_prod(const dense_matrix &m,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t out_layout = matrix_layout_t::L_NONE) const;
	std::shared_ptr<scalar_variable> aggregate(const bulk_operate &op) const;

	dense_matrix::ptr mapply2(const dense_matrix &m,
			bulk_operate::const_ptr op) const;
	dense_matrix::ptr sapply(bulk_uoperate::const_ptr op) const;
	dense_matrix::ptr apply(apply_margin margin,
			arr_apply_operate::const_ptr op) const;

	dense_matrix::ptr scale_cols(std::shared_ptr<const vector> vals) const;
	dense_matrix::ptr scale_rows(std::shared_ptr<const vector> vals) const;
	dense_matrix::ptr cast_ele_type(const scalar_type &type) const;

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

	std::shared_ptr<vector> row_sum() const;
	std::shared_ptr<vector> col_sum() const;
	std::shared_ptr<vector> row_norm2() const;
	std::shared_ptr<vector> col_norm2() const;

	dense_matrix::ptr append_cols(const std::vector<dense_matrix::ptr> &mats);

	template<class T>
	dense_matrix::ptr multiply_scalar(T val) const {
		scalar_variable::ptr var(new scalar_variable_impl<T>(val));
		return _multiply_scalar(var);
	}

	double norm2() const;
};

template<class T>
static dense_matrix operator*(const dense_matrix &m, T val)
{
	dense_matrix::ptr ret = m.multiply_scalar<T>(val);
	assert(ret);
	return *ret;
}

template<class T>
static dense_matrix operator*(T val, const dense_matrix &m)
{
	dense_matrix::ptr ret = m.multiply_scalar<T>(val);
	assert(ret);
	return *ret;
}

static dense_matrix operator*(const dense_matrix &m1, const dense_matrix &m2)
{
	dense_matrix::ptr ret = m1.multiply(m2);
	assert(ret);
	return *ret;
}

static dense_matrix operator+(const dense_matrix &m1, const dense_matrix &m2)
{
	dense_matrix::ptr ret = m1.add(m2);
	assert(ret);
	return *ret;
}

static dense_matrix operator-(const dense_matrix &m1, const dense_matrix &m2)
{
	dense_matrix::ptr ret = m1.minus(m2);
	assert(ret);
	return *ret;
}

template<class T>
static T as_scalar(const dense_matrix &m)
{
	assert(m.get_type() == get_scalar_type<T>());
	m.materialize_self();
	assert(m.is_in_mem());
	detail::mem_matrix_store::const_ptr mem_m
		= detail::mem_matrix_store::cast(m.get_raw_store());
	return mem_m->get<T>(0, 0);
}

static inline dense_matrix t(const dense_matrix &m)
{
	dense_matrix::ptr ret = m.transpose();
	assert(ret);
	return *ret;
}

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
		return ptr(new col_vec(mat->store));
	}
	template<class T>
	static ptr create_randu(size_t len) {
		dense_matrix::ptr mat = dense_matrix::create_randu<T>(0, 1, len, 1,
				matrix_layout_t::L_COL);
		return ptr(new col_vec(mat->store));
	}

	col_vec(): dense_matrix(NULL) {
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

}

#endif
