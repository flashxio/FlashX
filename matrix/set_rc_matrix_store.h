#ifndef __FM_SET_RC_MATRIX_STORE_H__
#define __FM_SET_RC_MATRIX_STORE_H__

/*
 * Copyright 2017 Open Connectome Project (http://openconnecto.me)
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

#include <vector>

#include "mem_matrix_store.h"
#include "materialize.h"

namespace fm
{

namespace detail
{

/*
 * This class is passed to mapply_matrix_store to replace
 * some rows of a tall matrix.
 */
class set_row_mapply_op: public portion_mapply_op
{
	std::shared_ptr<std::vector<off_t> > idx;
	mem_row_matrix_store::const_ptr rows;
public:
	set_row_mapply_op(std::shared_ptr<std::vector<off_t> > idx,
			mem_row_matrix_store::const_ptr rows,
			size_t num_rows, size_t num_cols,
			const scalar_type &type): portion_mapply_op(num_rows, num_cols,
				type) {
		this->idx = idx;
		this->rows = rows;
	}

	virtual portion_mapply_op::const_ptr transpose() const;

	virtual void run(
			const std::vector<std::shared_ptr<const local_matrix_store> > &ins,
			local_matrix_store &out) const;

	virtual std::string to_string(
			const std::vector<matrix_store::const_ptr> &mats) const {
		return mats[0]->get_name() + "[row_idxs]=" + rows->get_name();
	}
};

/*
 * This class is passed to mapply_matrix_store to replace
 * some cols of a wide matrix.
 */
class set_col_mapply_op: public portion_mapply_op
{
	std::shared_ptr<std::vector<off_t> > idx;
	mem_col_matrix_store::const_ptr cols;
public:
	set_col_mapply_op(std::shared_ptr<std::vector<off_t> > idx,
			mem_col_matrix_store::const_ptr cols,
			size_t num_rows, size_t num_cols,
			const scalar_type &type): portion_mapply_op(num_rows, num_cols,
				type) {
		this->idx = idx;
		this->cols = cols;
	}

	virtual portion_mapply_op::const_ptr transpose() const;

	virtual void run(
			const std::vector<std::shared_ptr<const local_matrix_store> > &ins,
			local_matrix_store &out) const;

	virtual std::string to_string(
			const std::vector<matrix_store::const_ptr> &mats) const {
		return mats[0]->get_name() + "[col_idxs]=" + cols->get_name();
	}
};

}

}

#endif
