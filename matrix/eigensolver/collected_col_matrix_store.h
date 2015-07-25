#ifndef __COLLECTED_COL_MATRIX_STORE_H__
#define __COLLECTED_COL_MATRIX_STORE_H__

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

#include "virtual_matrix_store.h"
#include "local_matrix_store.h"
#include "EM_object.h"

namespace fm
{

namespace eigen
{

class collected_matrix_store: public detail::virtual_matrix_store, public detail::EM_object
{
	// Indicates a row/column in the original matrices.
	struct rc_idx {
		size_t mat_idx;
		size_t inner_idx;

		rc_idx(size_t mat_idx, size_t inner_idx) {
			this->mat_idx = mat_idx;
			this->inner_idx = inner_idx;
		}
	};

	int num_nodes;
	std::vector<rc_idx> rc_idxs;
	std::vector<detail::matrix_store::const_ptr> orig_mats;
	detail::matrix_store::const_ptr merged_mat;

	collected_matrix_store(
			const std::vector<detail::matrix_store::const_ptr> &mats,
			size_t num_cols);
	collected_matrix_store(size_t num_rows, size_t num_cols, bool in_mem,
			const scalar_type &type, size_t num_nodes);
public:
	typedef std::shared_ptr<collected_matrix_store> ptr;

	static ptr create(const std::vector<detail::matrix_store::const_ptr> &mats,
			size_t num_cols);

	virtual std::string get_name() const {
		assert(merged_mat);
		return merged_mat->get_name();
	}

	virtual detail::matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const;
	virtual detail::matrix_store::const_ptr transpose() const;

	virtual detail::local_matrix_store::const_ptr get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows, size_t num_cols,
			detail::portion_compute::ptr compute) const {
		assert(merged_mat);
		return merged_mat->get_portion_async(start_row, start_col,
				num_rows, num_cols, compute);
	}
	virtual detail::local_matrix_store::const_ptr get_portion(size_t start_row,
			size_t start_col, size_t num_rows, size_t num_cols) const {
		assert(merged_mat);
		return merged_mat->get_portion(start_row, start_col, num_rows, num_cols);
	}
	virtual std::pair<size_t, size_t> get_portion_size() const {
		assert(merged_mat);
		return merged_mat->get_portion_size();
	}

	virtual int get_num_nodes() const {
		return num_nodes;
	}

	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		return merged_mat->get_underlying_mats();
	}

	matrix_layout_t store_layout() const {
		assert(merged_mat);
		return merged_mat->store_layout();
	}

	virtual matrix_store::ptr materialize() const {
		return static_cast<const detail::virtual_matrix_store *>(
				merged_mat.get())->materialize();
	}

	virtual std::vector<safs::io_interface::ptr> create_ios() const {
		assert(merged_mat);
		const detail::EM_object *obj = dynamic_cast<const detail::EM_object *>(
				merged_mat.get());
		assert(obj);
		return obj->create_ios();
	}
};

}

}

#endif
