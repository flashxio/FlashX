#ifndef __MAPPLY_MATRIX_STORE_H__
#define __MAPPLY_MATRIX_STORE_H__

/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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
#include "dense_matrix.h"
#include "mem_matrix_store.h"

namespace fm
{

namespace detail
{

class portion_mapply_op;

class mapply_matrix_store: public virtual_matrix_store
{
	matrix_layout_t layout;
	std::vector<mem_matrix_store::const_ptr> in_mats;
	portion_mapply_op::const_ptr op;
	// The materialized result matrix.
	mem_matrix_store::const_ptr res;
public:
	typedef std::shared_ptr<const mapply_matrix_store> const_ptr;

	mapply_matrix_store(
			const std::vector<mem_matrix_store::const_ptr> &in_mats,
			portion_mapply_op::const_ptr op, matrix_layout_t layout,
			size_t nrow, size_t ncol);

	virtual void materialize_self() const;

	virtual matrix_store::ptr materialize() const;

	using virtual_matrix_store::get;
	virtual const char *get(size_t row, size_t col) const;

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
	virtual std::pair<size_t, size_t> get_portion_size() const {
		return in_mats.front()->get_portion_size();
	}
	virtual int get_num_nodes() const {
		return in_mats.front()->get_num_nodes();
	}

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_layout_t store_layout() const {
		return layout;
	}
};

}

}

#endif
