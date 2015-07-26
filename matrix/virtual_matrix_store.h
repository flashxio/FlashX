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

namespace detail
{

/*
 * This class is the base class for the classes whose data isn't stored
 * physically. The data of these classes materializes data when being
 * requested. All these classes are read-only.
 */
class virtual_matrix_store: public matrix_store
{
public:
	typedef std::shared_ptr<const virtual_matrix_store> const_ptr;

	static const_ptr cast(matrix_store::const_ptr mat) {
		assert(mat->is_virtual());
		return std::static_pointer_cast<const virtual_matrix_store>(mat);
	}

	virtual_matrix_store(size_t nrow, size_t ncol, bool in_mem,
			const scalar_type &type): matrix_store(nrow, ncol, in_mem, type) {
	}

	virtual bool is_virtual() const {
		return true;
	}

	virtual matrix_store::const_ptr materialize() const = 0;

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
	virtual std::shared_ptr<local_matrix_store> get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, std::shared_ptr<portion_compute> compute) {
		assert(0);
		return std::shared_ptr<local_matrix_store>();
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
