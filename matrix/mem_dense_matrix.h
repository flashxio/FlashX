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
class mem_dense_matrix;

class mem_dense_matrix: public dense_matrix
{
public:
	typedef std::shared_ptr<mem_dense_matrix> ptr;
	typedef std::shared_ptr<const mem_dense_matrix> const_ptr;
private:
	mem_dense_matrix(detail::matrix_store::const_ptr store): dense_matrix(store) {
	}

	static ptr _create_rand(const scalar_variable &min, const scalar_variable &max,
			size_t nrow, size_t ncol, matrix_layout_t layout, int num_nodes);
	static ptr _create_const(scalar_variable::ptr val, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes);

	void inner_prod_tall(const detail::mem_matrix_store &m,
			const bulk_operate &left_op, const bulk_operate &right_op,
			detail::mem_matrix_store &res) const;
	void inner_prod_wide(const detail::mem_matrix_store &m,
			const bulk_operate &left_op, const bulk_operate &right_op,
			detail::mem_matrix_store &res) const;
	bool verify_inner_prod(const dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op) const;
public:
	/*
	 * We may create a matrix optimized for SMP, and we may also create a matrix
	 * optimized for NUMA. By default, it creates a SMP matrix.
	 */
	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, int num_nodes = -1);
	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, const set_operate &op, int num_nodes = -1);
	static ptr create(detail::mem_matrix_store::const_ptr store) {
		return ptr(new mem_dense_matrix(store));
	}

	template<class T>
	static ptr create_rand(T _min, T _max, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes = -1) {
		scalar_variable_impl<T> min(_min);
		scalar_variable_impl<T> max(_max);
		return _create_rand(min, max, nrow, ncol, layout, num_nodes);
	}

	template<class T>
	static ptr create_const(T _val, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes = -1) {
		scalar_variable::ptr val(new scalar_variable_impl<T>(_val));
		return _create_const(val, nrow, ncol, layout, num_nodes);
	}

	static ptr cast(dense_matrix::ptr m);
	static const_ptr cast(dense_matrix::const_ptr m);

	int get_num_nodes() const {
		return ((const detail::mem_matrix_store &) get_data()).get_num_nodes();
	}

	virtual dense_matrix::ptr clone() const {
		return ptr(new mem_dense_matrix(get_raw_store()));
	}

	virtual std::shared_ptr<vector> get_col(off_t idx) const;
	virtual std::shared_ptr<vector> get_row(off_t idx) const;
	virtual dense_matrix::ptr get_cols(const std::vector<off_t> &idxs) const;
	virtual dense_matrix::ptr get_rows(const std::vector<off_t> &idxs) const;

	virtual dense_matrix::ptr transpose() const;

	virtual dense_matrix::ptr inner_prod(const dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op,
			matrix_layout_t out_layout) const;
	virtual std::shared_ptr<scalar_variable> aggregate(
			const bulk_operate &op) const;
	/*
	 * A subclass should define this method for element-wise operations.
	 */
	virtual dense_matrix::ptr mapply2(const dense_matrix &m,
			const bulk_operate &op) const;
	virtual dense_matrix::ptr sapply(const bulk_uoperate &op) const;
	virtual dense_matrix::ptr apply(apply_margin margin,
			arr_apply_operate::const_ptr op) const;
	virtual dense_matrix::ptr scale_cols(
			std::shared_ptr<const mem_vector> vals) const;
	virtual dense_matrix::ptr scale_rows(
			std::shared_ptr<const mem_vector> vals) const;

	template<class T>
	T get(size_t row, size_t col) const {
		const detail::mem_matrix_store &store
			= (const detail::mem_matrix_store &) get_data();
		return store.get<T>(row, col);
	}
};

}

#endif
