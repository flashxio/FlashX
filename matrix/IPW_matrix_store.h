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

#include "sink_matrix.h"
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
class IPW_matrix_store: public sink_store
{
	matrix_store::const_ptr left_mat;
	matrix_store::const_ptr right_mat;
	bulk_operate::const_ptr left_op;
	bulk_operate::const_ptr right_op;
	std::shared_ptr<const portion_mapply_op> portion_op;
	matrix_layout_t layout;
	// We often need to check the underlying matrices of a sink matrix
	// when materializing a virtual matrix.
	// This caches the underlying matrices that are used to compute
	// this sink matrix. We rarely materialize any non-sink matrices,
	// so the underlying matrices usually remain the same.
	std::unordered_map<size_t, size_t> underlying;

	matrix_store::const_ptr get_combine_res() const;
	IPW_matrix_store(matrix_store::const_ptr left, matrix_store::const_ptr right,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t layout);
	// This constructor is used to construct the transpose of the matrix.
	IPW_matrix_store(matrix_store::const_ptr left, matrix_store::const_ptr right,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t layout,
			std::shared_ptr<const portion_mapply_op> portion_op);
public:
	static ptr create(matrix_store::const_ptr left, matrix_store::const_ptr right,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t layout = matrix_layout_t::L_NONE) {
		ptr ret(new IPW_matrix_store(left, right, left_op, right_op, layout));
		sink_store::register_sink_matrices(ret);
		return ret;
	}
	virtual size_t get_data_id() const;

	virtual void inc_dag_ref(size_t id) {
		if (left_mat)
			const_cast<matrix_store &>(*left_mat).inc_dag_ref(get_data_id());
		if (right_mat)
			const_cast<matrix_store &>(*right_mat).inc_dag_ref(get_data_id());
		// We don't need to increase the ref count of a sink matrix
		// because we never get a portion from a sink matrix.
	}
	virtual void reset_dag_ref() {
		if (left_mat)
			const_cast<matrix_store &>(*left_mat).reset_dag_ref();
		if (right_mat)
			const_cast<matrix_store &>(*right_mat).reset_dag_ref();
	}

	virtual bool has_materialized() const;

	virtual void materialize_self() const;

	virtual std::vector<virtual_matrix_store::const_ptr> get_compute_matrices() const;
	virtual matrix_store::const_ptr get_result() const {
		if (has_materialized())
			return get_combine_res();
		else
			return matrix_store::const_ptr();
	}
	virtual matrix_store::const_ptr materialize(bool in_mem,
		int num_nodes) const;

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_layout_t store_layout() const {
		// TODO what is the right layout?
		return layout;
	}

	std::unordered_map<size_t, size_t> get_underlying_mats() const;

	virtual std::string get_name() const;
};

}

}

#endif
