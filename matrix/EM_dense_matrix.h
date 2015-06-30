#ifndef __EM_DENSE_MATRIX_H__
#define __EM_DENSE_MATRIX_H__

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
#include <boost/format.hpp>

#include "log.h"

#include "bulk_operate.h"
#include "matrix_store.h"
#include "EM_object.h"
#include "mem_worker_thread.h"

namespace fm
{

namespace detail
{

class local_matrix_store;
class mem_matrix_store;

class EM_matrix_store: public matrix_store, public EM_object
{
	matrix_layout_t layout;
	file_holder::ptr holder;
	io_set::ptr ios;

	EM_matrix_store(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type);
	EM_matrix_store(file_holder::ptr holder, io_set::ptr ios,
			size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type): matrix_store(nrow, ncol, false, type) {
		this->layout = layout;
		this->holder = holder;
		this->ios = ios;
	}
public:
	typedef std::shared_ptr<EM_matrix_store> ptr;
	typedef std::shared_ptr<const EM_matrix_store> const_ptr;

	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type) {
		return ptr(new EM_matrix_store(nrow, ncol, layout, type));
	}

	static ptr cast(matrix_store::ptr store) {
		return std::dynamic_pointer_cast<EM_matrix_store>(store);
	}

	static const_ptr cast(matrix_store::const_ptr store) {
		return std::dynamic_pointer_cast<const EM_matrix_store>(store);
	}

	/*
	 * Load the EM matrix to memory.
	 */
	std::shared_ptr<mem_matrix_store> load() const;

	virtual void reset_data();
	virtual void set_data(const set_operate &op);

	virtual matrix_layout_t store_layout() const {
		return layout;
	}

	virtual matrix_store::const_ptr transpose() const {
		matrix_layout_t new_layout;
		if (layout == matrix_layout_t::L_ROW)
			new_layout = matrix_layout_t::L_COL;
		else
			new_layout = matrix_layout_t::L_ROW;
		return matrix_store::const_ptr(new EM_matrix_store(holder, ios,
					get_num_cols(), get_num_rows(), new_layout, get_type()));
	}

	virtual std::vector<safs::io_interface::ptr> create_ios() const;

	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols);

	virtual std::pair<size_t, size_t> get_portion_size() const;
	virtual std::shared_ptr<const local_matrix_store> get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, portion_compute::ptr compute) const;
	virtual std::shared_ptr<local_matrix_store> get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, portion_compute::ptr compute);
	virtual void write_portion_async(
			std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col);

	virtual matrix_store::const_ptr append_cols(
			const std::vector<matrix_store::const_ptr> &mats) const;

	virtual std::shared_ptr<const vec_store> get_col_vec(off_t idx) const;
	virtual std::shared_ptr<const vec_store> get_row_vec(off_t idx) const;
};

}

}

#endif
