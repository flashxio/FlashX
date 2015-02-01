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

namespace fm
{

class bulk_operate;
class bulk_uoperate;
class set_operate;
class scalar_type;

enum matrix_type
{
	VECTOR,
	DENSE,
	SPARSE,
};

enum matrix_layout_t
{
	L_COL,
	L_ROW,
};

/*
 * The matrix header contains the basic information about the matrix
 * when the matrix is stored on disks.
 */
struct matrix_header
{
	union {
		struct info {
			matrix_type type;
			size_t entry_size;
			size_t nrows;
			size_t ncols;
			matrix_layout_t layout;
			// It doesn't include the entry type, which should be determined by
			// users at runtime.
		} d;

		char page[4096];
	} u;
};

class dense_matrix
{
	size_t nrow;
	size_t ncol;
	size_t entry_size;
	bool in_mem;
protected:
	dense_matrix(size_t nrow, size_t ncol, size_t entry_size, bool in_mem) {
		this->nrow = nrow;
		this->ncol = ncol;
		this->entry_size = entry_size;
		this->in_mem = in_mem;
	}
	virtual bool verify_aggregate(const bulk_operate &op, scalar_type &res) const;
	virtual bool verify_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const;
	virtual bool verify_mapply2(const dense_matrix &m,
			const bulk_operate &op) const;
public:
	typedef std::shared_ptr<dense_matrix> ptr;

	static ptr create(size_t nrow, size_t ncol, size_t entry_size,
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

	/*
	 * Right now we can only verify the type by looking at the entry size.
	 * TODO we might need to verify it more thoroughly.
	 */
	template<class T>
	bool is_type() const {
		return sizeof(T) == get_entry_size();
	}

	prim_type get_type() const {
		// TODO this is a temporary solution.
		switch(get_entry_size()) {
			case sizeof(bool): return prim_type::P_BOOL;
			case sizeof(int): return prim_type::P_INTEGER;
			case sizeof(double): return prim_type::P_DOUBLE;
			default:
					fprintf(stderr, "wrong type");
					return prim_type::NUM_TYPES;
		}
	}

	/**
	 * The shape of a matrix: a tall matrix or a wide matrix.
	 * We care about the shape of a large matrix. We deal with a tall matrix
	 * different from a wide matrix.
	 */
	bool is_wide() const {
		return get_num_cols() > get_num_rows();
	}

	virtual dense_matrix::ptr clone() const = 0;
	/**
	 * This method will create a new matrix with different dimensions and
	 * data layout. But the new matrix shares the same data as the original
	 * matrix.
	 */
	virtual dense_matrix::ptr conv2(size_t nrow, size_t ncol, bool byrow) const = 0;
	virtual dense_matrix::ptr transpose() const = 0;

	virtual void reset_data() = 0;
	virtual void set_data(const set_operate &op) = 0;

	virtual dense_matrix::ptr inner_prod(const dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op) const = 0;
	virtual bool aggregate(const bulk_operate &op, scalar_type &res) const = 0;
	/*
	 * A subclass should define this method for element-wise operations.
	 */
	virtual dense_matrix::ptr mapply2(const dense_matrix &m,
			const bulk_operate &op) const = 0;
	virtual dense_matrix::ptr sapply(const bulk_uoperate &op) const = 0;
};

}

#endif
