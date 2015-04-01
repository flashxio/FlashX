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

class dense_matrix
{
	size_t nrow;
	size_t ncol;
	bool in_mem;
	// This is kind of redundant because we can always get the entry size
	// from the type. However, getting the entry size of the type requires
	// to call a virtual method. Storing the entry size here can avoid
	// the function call. It doesn't increase the size of the data structure
	// due to the data alignment by the compiler.
	int entry_size;
	// The type is a reference. It makes the dense matrix object uncopiable.
	// Maybe this is what we want.
	const scalar_type &type;
protected:
	dense_matrix(size_t nrow, size_t ncol, const scalar_type &_type,
			bool in_mem): type(_type) {
		this->nrow = nrow;
		this->ncol = ncol;
		this->in_mem = in_mem;
		this->entry_size = type.get_size();
	}
	virtual bool verify_aggregate(const bulk_operate &op) const;
	virtual bool verify_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const;
	virtual bool verify_mapply2(const dense_matrix &m,
			const bulk_operate &op) const;
	virtual bool verify_apply(apply_margin margin, const arr_apply_operate &op) const;
public:
	typedef std::shared_ptr<dense_matrix> ptr;

	static ptr create(size_t nrow, size_t ncol, const scalar_type &type,
			matrix_layout_t layout, bool in_mem);
	static ptr load(const std::string &file_name);

	virtual ~dense_matrix() {
	}

	bool write_header(FILE *f) const;
	virtual bool write2file(const std::string &file_name) const = 0;

	size_t get_entry_size() const {
		return entry_size;
	}

	size_t get_num_rows() const {
		return nrow;
	}

	size_t get_num_cols() const {
		return ncol;
	}

	virtual matrix_layout_t store_layout() const = 0;

	/**
	 * Resize the matrix to a new dimension.
	 * The current matrix must have enough elements. Otherwise, resizing
	 * fails.
	 */
	bool resize(size_t new_nrow, size_t new_ncol) {
		if (nrow * ncol < new_nrow * new_ncol)
			return false;
		else {
			nrow = new_nrow;
			ncol = new_ncol;
			return true;
		}
	}

	bool is_in_mem() const {
		return in_mem;
	}

	template<class T>
	bool is_type() const {
		return type.get_type() == fm::get_type<T>();
	}

	const scalar_type &get_type() const {
		return type;
	}

	/**
	 * The shape of a matrix: a tall matrix or a wide matrix.
	 * We care about the shape of a large matrix. We deal with a tall matrix
	 * different from a wide matrix.
	 */
	bool is_wide() const {
		return get_num_cols() > get_num_rows();
	}

	virtual dense_matrix::ptr shallow_copy() const = 0;
	virtual dense_matrix::ptr deep_copy() const = 0;
	/**
	 * This method will create a new matrix with different dimensions and
	 * data layout. But the new matrix shares the same data as the original
	 * matrix.
	 */
	virtual dense_matrix::ptr conv2(size_t nrow, size_t ncol, bool byrow) const = 0;
	virtual dense_matrix::ptr transpose() const = 0;

	virtual void reset_data() = 0;
	virtual void set_data(const set_operate &op) = 0;
	virtual bool copy_from(const dense_matrix &mat) = 0;

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
