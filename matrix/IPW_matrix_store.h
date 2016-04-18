/*
 * Copyright 2016 Open Connectome Project (http://openconnecto.me)
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

#include "io_interface.h"

#include "virtual_matrix_store.h"
#include "bulk_operate.h"
#include "EM_object.h"

#ifndef __FM_IPW_MATRIX_STORE_H__
#define __FM_IPW_MATRIX_STORE_H__

namespace fm
{
namespace detail
{

class portion_mapply_op;

/*
 * This matrix store is to enable lazy evaluation on the inner product
 * on a wide matrix.
 */
class IPW_matrix_store: public virtual_matrix_store, public EM_object
{
	matrix_store::const_ptr left_mat;
	matrix_store::const_ptr right_mat;
	bulk_operate::const_ptr left_op;
	bulk_operate::const_ptr right_op;
	std::shared_ptr<portion_mapply_op> portion_op;
	matrix_layout_t layout;

	matrix_store::ptr get_combine_res() const;
public:
	IPW_matrix_store(matrix_store::const_ptr left, matrix_store::const_ptr right,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t layout);

	bool has_materialized() const;

	virtual void materialize_self() const;

	virtual matrix_store::const_ptr materialize(bool in_mem,
		int num_nodes) const;

	virtual std::shared_ptr<const vec_store> get_col_vec(off_t idx) const;
	virtual std::shared_ptr<const vec_store> get_row_vec(off_t idx) const;
	virtual matrix_store::const_ptr get_cols(const std::vector<off_t> &idxs) const;
	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const;

	using virtual_matrix_store::get_portion;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t id) const;
	using virtual_matrix_store::get_portion_async;
	virtual async_cres_t get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, std::shared_ptr<portion_compute> compute) const;

	virtual matrix_store::const_ptr transpose() const;

	virtual int get_portion_node_id(size_t id) const {
		// If both matrices are stored in NUMA memory, the portion must be
		// stored on the same NUMA node. Otherwise, we need to return
		// the node Id from the matrix stored in NUMA.
		if (left_mat->get_num_nodes() > 0)
			return left_mat->get_portion_node_id(id);
		else
			return right_mat->get_portion_node_id(id);
	}

	virtual std::pair<size_t, size_t> get_portion_size() const {
		// TODO what do we do if the two input matrices aren't both
		// in memory or both on SSDS.
		assert(left_mat->is_in_mem() == right_mat->is_in_mem());
		return left_mat->get_portion_size();
	}

	virtual int get_num_nodes() const {
		if (left_mat->get_num_nodes() > 0)
			return left_mat->get_num_nodes();
		else
			return right_mat->get_num_nodes();
	}

	virtual matrix_layout_t store_layout() const {
		return layout;
	}

	virtual std::vector<safs::io_interface::ptr> create_ios() const;

	virtual std::string get_name() const;
	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const;
};

}

}

#endif
