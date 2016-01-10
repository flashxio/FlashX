#ifndef __VIRTUAL_MATRIX_STORE_H__
#define __VIRTUAL_MATRIX_STORE_H__

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

#include "comm_exception.h"

#include "mem_matrix_store.h"

namespace fm
{

enum materialize_level
{
	// Materialize in the CPU cache.
	// This delivers the best performance when materializing a sequence of
	// matrix operations. However, if a virtual matrix appears in multiple
	// operands, it requires to materialize the virtual matrix multiple times.
	// This is the default level for materialization.
	MATER_CPU,
	// Materialize the entire portion of a virtual matrix.
	// This requires to allocate a much larger memory buffer to keep
	// the materialization result, but it can avoid recomputation if a virtual
	// matrix is used in multiple operands of a function.
	MATER_MEM,
	// Materialize the entire matrix.
	// This is especially expensive if the matrix is stored on disks.
	// We should avoid it as much as possible.
	MATER_FULL,
};

namespace detail
{

/*
 * This class is the base class for the classes whose data isn't stored
 * physically. The data of these classes materializes data when being
 * requested. All these classes are read-only.
 */
class virtual_matrix_store: public matrix_store
{
	materialize_level mater_level;
public:
	typedef std::shared_ptr<const virtual_matrix_store> const_ptr;

	static const_ptr cast(matrix_store::const_ptr mat) {
		assert(mat->is_virtual());
		return std::static_pointer_cast<const virtual_matrix_store>(mat);
	}

	virtual_matrix_store(size_t nrow, size_t ncol, bool in_mem,
			const scalar_type &type): matrix_store(nrow, ncol, in_mem, type) {
		this->mater_level = materialize_level::MATER_CPU;
	}

	virtual void set_materialize_level(materialize_level level,
			detail::matrix_store::ptr materialize_buf) {
		this->mater_level = level;
	}

	materialize_level get_materialize_level() const {
		return mater_level;
	}

	virtual bool is_virtual() const {
		return true;
	}

	/*
	 * When we materialize the matrix, we can specify where the materialized
	 * matrix is stored.
	 */
	virtual matrix_store::const_ptr materialize(bool in_mem,
			int num_nodes) const = 0;

	virtual void reset_data() {
		assert(0);
	}

	virtual void set_data(const set_operate &op) {
		assert(0);
	}

	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id) {
		assert(0);
		return std::shared_ptr<local_matrix_store>();
	}
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) {
		assert(0);
		return std::shared_ptr<local_matrix_store>();
	}
	virtual async_res_t get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, std::shared_ptr<portion_compute> compute) {
		assert(0);
		return async_res_t();
	}
	virtual void write_portion_async(
			std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col) {
		assert(0);
	}

	virtual bool write2file(const std::string &file_name) const {
		assert(0);
		return false;
	}
};

}

}

#endif
