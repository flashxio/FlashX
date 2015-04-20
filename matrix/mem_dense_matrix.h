#ifndef __MEM_DENSE_MATRIX_H__
#define __MEM_DENSE_MATRIX_H__

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
#include <assert.h>
#include <string.h>
#include <malloc.h>

#include "common.h"
#include "dense_matrix.h"
#include "bulk_operate.h"
#include "matrix_config.h"
#include "matrix_header.h"
#include "raw_data_array.h"
#include "mem_matrix_store.h"

namespace fm
{

class mem_vector;

class mem_dense_matrix: public dense_matrix
{
public:
	typedef std::shared_ptr<mem_dense_matrix> ptr;
private:
	mem_dense_matrix(detail::matrix_store::const_ptr store): dense_matrix(store) {
	}

	static ptr _create_rand(const scalar_variable &min, const scalar_variable &max,
			size_t nrow, size_t ncol, matrix_layout_t layout);
	static ptr _create_const(const scalar_variable &val, size_t nrow, size_t ncol,
			matrix_layout_t layout);

	void inner_prod_tall(const detail::mem_matrix_store &m,
			const bulk_operate &left_op, const bulk_operate &right_op,
			detail::mem_matrix_store &res) const;
	void inner_prod_wide(const detail::mem_matrix_store &m,
			const bulk_operate &left_op, const bulk_operate &right_op,
			detail::mem_matrix_store &res) const;
	bool verify_inner_prod(const dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op) const;
public:
	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type);
	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, const set_operate &op);

	template<class T>
	static ptr create_rand(T _min, T _max, size_t nrow, size_t ncol,
			matrix_layout_t layout) {
		scalar_variable_impl<T> min(_min);
		scalar_variable_impl<T> max(_max);
		return _create_rand(min, max, nrow, ncol, layout);
	}

	template<class T>
	static ptr create_const(T _val, size_t nrow, size_t ncol,
			matrix_layout_t layout) {
		scalar_variable_impl<T> val(_val);
		return _create_const(val, nrow, ncol, layout);
	}

	static ptr cast(dense_matrix::ptr m);

	virtual dense_matrix::ptr get_cols(const std::vector<off_t> &idxs) const;
	virtual dense_matrix::ptr get_rows(const std::vector<off_t> &idxs) const;

	virtual dense_matrix::ptr transpose() const;

	virtual dense_matrix::ptr inner_prod(const dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op) const;
	virtual std::shared_ptr<scalar_variable> aggregate(
			const bulk_operate &op) const;
	/*
	 * A subclass should define this method for element-wise operations.
	 */
	virtual dense_matrix::ptr mapply2(const dense_matrix &m,
			const bulk_operate &op) const;
	virtual dense_matrix::ptr sapply(const bulk_uoperate &op) const;
	virtual dense_matrix::ptr apply(apply_margin margin,
			const arr_apply_operate &op) const;

	template<class T>
	T get(size_t row, size_t col) const {
		const detail::mem_matrix_store &store
			= (const detail::mem_matrix_store &) get_data();
		return store.get<T>(row, col);
	}
};

template<class EntryType>
class type_mem_dense_matrix
{
	mem_dense_matrix::ptr m;

	type_mem_dense_matrix(size_t nrow, size_t ncol, matrix_layout_t layout) {
		m = mem_dense_matrix::create(nrow, ncol, layout,
				get_scalar_type<EntryType>());
	}

	type_mem_dense_matrix(size_t nrow, size_t ncol, matrix_layout_t layout,
			const type_set_operate<EntryType> &op) {
		m = mem_dense_matrix::create(nrow, ncol, layout,
				get_scalar_type<EntryType>(), op);
	}

	type_mem_dense_matrix(mem_dense_matrix::ptr m) {
		this->m = m;
	}
public:
	typedef std::shared_ptr<type_mem_dense_matrix<EntryType> > ptr;

	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout) {
		return ptr(new type_mem_dense_matrix<EntryType>(nrow, ncol, layout));
	}

	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const type_set_operate<EntryType> &op) {
		return ptr(new type_mem_dense_matrix<EntryType>(nrow, ncol, layout, op));
	}

	static ptr create(mem_dense_matrix::ptr m) {
		assert(m->get_type() == get_scalar_type<EntryType>());
		return ptr(new type_mem_dense_matrix<EntryType>(m));
	}

	size_t get_num_rows() const {
		return m->get_num_rows();
	}

	size_t get_num_cols() const {
		return m->get_num_cols();
	}

	EntryType get(size_t row, size_t col) const {
		const detail::mem_matrix_store &mem_store
			= (const detail::mem_matrix_store &) m->get_data();
		return mem_store.get<EntryType>(row, col);
	}

	const mem_dense_matrix::ptr get_matrix() const {
		return m;
	}
};

typedef type_mem_dense_matrix<int> I_mem_dense_matrix;
typedef type_mem_dense_matrix<double> D_mem_dense_matrix;

}

#endif
