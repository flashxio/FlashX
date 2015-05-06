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

#include <memory>

#include "generic_type.h"
#include "matrix_header.h"
#include "bulk_operate.h"
#include "matrix_store.h"

namespace fm
{

class bulk_operate;
class bulk_uoperate;
class set_operate;
class arr_apply_operate;

enum apply_margin
{
	MAR_ROW = 1,
	MAR_COL = 2,
};

/*
 * This class represents a dense matrix and is able to perform computation
 * on the matrix. However, this class can't modify the matrix data. The only
 * way to modify the matrix is to have the pointer point to another matrix.
 */
class dense_matrix
{
	detail::matrix_store::const_ptr store;

protected:
	dense_matrix(detail::matrix_store::const_ptr store) {
		this->store = store;
	}
	detail::matrix_store::const_ptr get_raw_store() const {
		return store;
	}

	virtual bool verify_aggregate(const bulk_operate &op) const;
	virtual bool verify_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const;
	virtual bool verify_mapply2(const dense_matrix &m,
			const bulk_operate &op) const;
	virtual bool verify_apply(apply_margin margin, const arr_apply_operate &op) const;
public:
	typedef std::shared_ptr<dense_matrix> ptr;
	typedef std::shared_ptr<const dense_matrix> const_ptr;

	static ptr create(size_t nrow, size_t ncol, const scalar_type &type,
			matrix_layout_t layout, bool in_mem);

	virtual ~dense_matrix() {
	}

	const detail::matrix_store &get_data() const {
		return *store;
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

	virtual dense_matrix::ptr get_cols(
			const std::vector<off_t> &idxs) const = 0;
	virtual dense_matrix::ptr get_rows(
			const std::vector<off_t> &idxs) const = 0;
	/*
	 * Clone the matrix.
	 * The class can't modify the matrix data that it points to, but it
	 * can modify the pointer. If someone changes in the pointer in the cloned
	 * matrix, it doesn't affect the current matrix.
	 */
	virtual dense_matrix::ptr clone() const = 0;

	virtual dense_matrix::ptr transpose() const = 0;

	virtual dense_matrix::ptr inner_prod(const dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op) const = 0;
	virtual std::shared_ptr<scalar_variable> aggregate(
			const bulk_operate &op) const = 0;
	/*
	 * A subclass should define this method for element-wise operations.
	 */
	virtual dense_matrix::ptr mapply2(const dense_matrix &m,
			const bulk_operate &op) const = 0;
	virtual dense_matrix::ptr sapply(const bulk_uoperate &op) const = 0;
	virtual dense_matrix::ptr apply(apply_margin margin,
			const arr_apply_operate &op) const = 0;
	virtual dense_matrix::ptr scale_cols(const mem_vector &vals) const = 0;

	dense_matrix::ptr multiply(const dense_matrix &mat) const {
		return inner_prod(mat, get_type().get_basic_ops().get_multiply(),
				get_type().get_basic_ops().get_add());
	}

	dense_matrix::ptr add(const dense_matrix &mat) const {
		return this->mapply2(mat, get_type().get_basic_ops().get_add());
	}

	template<class T>
	dense_matrix::ptr multiply_scalar(T val) const {
		assert(get_type() == get_scalar_type<T>());
		return this->sapply(bulk_uoperate_impl<multiply_uop<T>, T, T>(
					multiply_uop<T>(val)));
	}

	double norm2() const;
};

}

#endif
