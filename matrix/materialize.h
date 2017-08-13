#ifndef __FM_MATERIALIZE_H__
#define __FM_MATERIALIZE_H__

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

#include "matrix_store.h"

namespace fm
{

class dense_matrix;

namespace detail
{

class local_matrix_store;
class sparse_project_matrix_store;

class portion_mapply_op
{
	size_t out_num_rows;
	size_t out_num_cols;
	const scalar_type &type;
public:
	typedef std::shared_ptr<const portion_mapply_op> const_ptr;

	portion_mapply_op(size_t out_num_rows, size_t out_num_cols,
			const scalar_type &_type): type(_type) {
		this->out_num_rows = out_num_rows;
		this->out_num_cols = out_num_cols;
	}
	virtual ~portion_mapply_op() {
	}

	virtual portion_mapply_op::const_ptr transpose() const = 0;

	/*
	 * There are three versions of performing computation on the portions.
	 * The first version performs computation only on input portions;
	 * the second version performs computation on input portions and
	 * outputs only one matrix;
	 * the third version performs computation on input portions and outputs
	 * multiple matrices.
	 */

	virtual void run(
			const std::vector<std::shared_ptr<const local_matrix_store> > &ins) const;
	virtual void run(
			const std::vector<std::shared_ptr<const local_matrix_store> > &ins,
			local_matrix_store &out) const;
	virtual void run(
			const std::vector<std::shared_ptr<const local_matrix_store> > &ins,
			const std::vector<std::shared_ptr<local_matrix_store> > &outs) const;

	virtual std::string to_string(
			const std::vector<matrix_store::const_ptr> &mats) const = 0;

	// A portion operation may fail. If so, the portion operation should
	// notify the main thread that the operation failed by override this method.
	virtual bool is_success() const {
		return true;
	}

	/*
	 * Give a hint if this operation is aggregation, so we can optimize
	 * the backend accordingly. When this is an aggregation operation,
	 * the second `run' method has to be implemented.
	 */
	virtual bool is_agg() const {
		return false;
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		return true;
	}

	bool is_wide() const {
		return get_out_num_cols() > get_out_num_rows();
	}

	size_t get_out_num_rows() const {
		return out_num_rows;
	}
	size_t get_out_num_cols() const {
		return out_num_cols;
	}
	const scalar_type &get_output_type() const {
		return type;
	}
};

/*
 * These two functions return a virtual matrix that records the computation.
 */

std::shared_ptr<dense_matrix> mapply_portion(
		const std::vector<std::shared_ptr<const dense_matrix> > &mats,
		// A user can specify the layout of the output dense matrix.
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout,
		bool par_access = true);
matrix_store::ptr __mapply_portion_virtual(
		const std::vector<matrix_store::const_ptr> &store,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout,
		bool par_access = true);

/*
 * These three functions return a materialized matrix.
 * The first version determines the storage of the output matrix automatically.
 * The second version allows users to specify the storage for the output matrix.
 * The third version not only allows users to specify the storage for the output
 * matrix, but also allows multiple output matrices.
 */

matrix_store::ptr __mapply_portion(
		const std::vector<matrix_store::const_ptr> &mats,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout,
		bool par_access = true);
matrix_store::ptr __mapply_portion(
		const std::vector<matrix_store::const_ptr> &mats,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout,
		bool out_in_mem, int out_num_nodes, bool par_access = true);
bool __mapply_portion(
		const std::vector<matrix_store::const_ptr> &mats,
		portion_mapply_op::const_ptr op,
		const std::vector<matrix_store::ptr> &out_mats, bool par_access = true);

}

bool materialize(std::vector<std::shared_ptr<dense_matrix> > &mats,
		bool par_access = true, bool mater_self=true);

}

#endif
