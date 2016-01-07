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
#include <atomic>
#include <unordered_map>

#include "safs_file.h"

#include "matrix_header.h"
#include "generic_type.h"

namespace fm
{

class set_operate;

namespace detail
{

class portion_compute;
class local_matrix_store;
class vec_store;

typedef std::pair<bool, std::shared_ptr<local_matrix_store> > async_res_t;
typedef std::pair<bool, std::shared_ptr<const local_matrix_store> > async_cres_t;

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
protected:
	static std::atomic<size_t> mat_counter;
public:
	typedef std::shared_ptr<matrix_store> ptr;
	typedef std::shared_ptr<const matrix_store> const_ptr;

	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, int num_nodes, bool in_mem,
			safs::safs_file_group::ptr group = NULL);

	matrix_store(size_t nrow, size_t ncol, bool in_mem,
			const scalar_type &_type);

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

	/*
	 * This method gets underlying materialized matrix IDs and the number of
	 * elements in each of these materialized matrices.
	 */
	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const = 0;
	virtual std::string get_name() const = 0;

	virtual matrix_layout_t store_layout() const = 0;

	virtual void reset_data();
	virtual void set_data(const set_operate &op);

	virtual matrix_store::const_ptr transpose() const = 0;

	/*
	 * When matrix data is move to faster memory, data is moved in one chunk
	 * at a time. The chunk size is defined by a specific implementation of
	 * matrix store. Each chunk is assigned with an identifier, which is
	 * defined sequentially.
	 */
	size_t get_num_portions() const;
	virtual std::pair<size_t, size_t> get_portion_size() const = 0;
	/*
	 * These two versions get a portion of data from the matrix asynchronously.
	 * When a local matrix store is returned, it's not guaranteed that the
	 * data in the local matrix store is valid. A status in the returned value
	 * indicates whether the data is valid. If the data is invalid when it's
	 * returned from the two methods, the computation passed to
	 * these two methods are invoked when the portion of data is loaded
	 * in memory. During the time between returning from the methods and
	 * the portion of data becomes available, it's users' responsibility
	 * of keep the local matrix store alive.
	 */
	virtual async_cres_t get_portion_async(size_t start_row, size_t start_col,
			size_t num_rows, size_t num_cols,
			std::shared_ptr<portion_compute> compute) const = 0;
	virtual async_res_t get_portion_async(size_t start_row, size_t start_col,
			size_t num_rows, size_t num_cols,
			std::shared_ptr<portion_compute> compute) = 0;
	/*
	 * These versions fetches the portion of data. It's guaranteed that
	 * the data in the returned local matrix store is valid.
	 */
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const = 0;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) = 0;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t id) const;
	/*
	 * This gets the node Id of the specified portion.
	 */
	virtual int get_portion_node_id(size_t id) const = 0;
	virtual void write_portion_async(
			std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col) = 0;

	virtual matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const {
		matrix_store::const_ptr tm = transpose();
		matrix_store::const_ptr rows = tm->get_rows(idxs);
		if (rows == NULL)
			return matrix_store::const_ptr();
		else
			return rows->transpose();
	}
	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const {
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

	/*
	 * This allow users to enable/disable caching data portions.
	 * It's used by EM matrix and mapply virtual matrix.
	 */
	virtual void set_cache_portion(bool cache_portion) {
	}
};

}

}

#endif
