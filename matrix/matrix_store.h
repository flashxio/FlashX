#ifndef __MATRIX_STORE_H__
#define __MATRIX_STORE_H__

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

#include <memory>

#include "matrix_header.h"
#include "generic_type.h"

namespace fm
{

class set_operate;

namespace detail
{

class local_matrix_store;
class vec_store;

class matrix_store
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
public:
	typedef std::shared_ptr<matrix_store> ptr;
	typedef std::shared_ptr<const matrix_store> const_ptr;

	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, int num_nodes, bool in_mem);

	matrix_store(size_t nrow, size_t ncol, bool in_mem,
			const scalar_type &_type): type(_type) {
		this->nrow = nrow;
		this->ncol = ncol;
		this->in_mem = in_mem;
		this->entry_size = type.get_size();
	}

	virtual ~matrix_store() {
	}

	void resize(size_t num_rows, size_t num_cols) {
		this->nrow = num_rows;
		this->ncol = num_cols;
	}

	size_t get_num_rows() const {
		return nrow;
	}

	size_t get_num_cols() const {
		return ncol;
	}

	size_t get_entry_size() const {
		return entry_size;
	}

	const scalar_type &get_type() const {
		return type;
	}

	bool is_in_mem() const {
		return in_mem;
	}

	virtual int get_num_nodes() const {
		return -1;
	}

	/*
	 * The shape of a matrix: a tall matrix or a wide matrix.
	 * We care about the shape of a large matrix. We deal with a tall matrix
	 * different from a wide matrix.
	 */
	bool is_wide() const {
		return get_num_cols() > get_num_rows();
	}

	virtual matrix_layout_t store_layout() const = 0;

	virtual void reset_data() = 0;
	virtual void set_data(const set_operate &op) = 0;

	virtual matrix_store::const_ptr transpose() const = 0;
	virtual matrix_store::ptr conv2(matrix_layout_t layout) const {
		// TODO
		assert(0);
		return matrix_store::ptr();
	}

	/*
	 * When matrix data is move to faster memory, data is moved in one chunk
	 * at a time. The chunk size is defined by a specific implementation of
	 * matrix store. Each chunk is assigned with an identifier, which is
	 * defined sequentially.
	 */
	size_t get_num_portions() const;
	virtual std::pair<size_t, size_t> get_portion_size() const = 0;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const = 0;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) = 0;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t id) const;

	virtual matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const {
		// TODO
		assert(0);
		return matrix_store::const_ptr();
	}
	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const {
		// TODO
		assert(0);
		return matrix_store::const_ptr();
	}
	virtual std::shared_ptr<const vec_store> get_col_vec(off_t idx) const {
		assert(0);
		return std::shared_ptr<const vec_store>();
	}
	virtual std::shared_ptr<const vec_store> get_row_vec(off_t idx) const {
		assert(0);
		return std::shared_ptr<const vec_store>();
	}

	virtual bool is_virtual() const {
		return false;
	}
	virtual void materialize_self() const {
	}

	virtual matrix_store::const_ptr append_cols(
			const std::vector<matrix_store::const_ptr> &mats) const = 0;
};

}

}

#endif
