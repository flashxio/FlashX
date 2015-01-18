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

namespace fm
{

class bulk_operate;
class set_operate;

enum matrix_layout_t
{
	L_COL,
	L_ROW,
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
public:
	typedef std::shared_ptr<dense_matrix> ptr;

	static ptr create(size_t nrow, size_t ncol, size_t entry_size,
			matrix_layout_t layout, bool in_mem);

	virtual ~dense_matrix() {
	}

	size_t get_entry_size() const {
		return entry_size;
	}

	size_t get_num_rows() const {
		return nrow;
	}

	size_t get_num_cols() const {
		return ncol;
	}

	bool is_in_mem() const {
		return in_mem;
	}

	/**
	 * The shape of a matrix: a tall matrix or a wide matrix.
	 * We care about the shape of a large matrix. We deal with a tall matrix
	 * different from a wide matrix.
	 */
	bool is_wide() const {
		return get_num_cols() > get_num_rows();
	}

	virtual void reset_data() = 0;
	virtual void set_data(const set_operate &op) = 0;

	virtual bool verify_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const;
	virtual dense_matrix::ptr inner_prod(const dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op) const = 0;
};

}

#endif
