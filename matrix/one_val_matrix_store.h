#ifndef __ONE_VAL_MATRIX_STORE_H__
#define __ONE_VAL_MATRIX_STORE_H__

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
#include "NUMA_mapper.h"

namespace fm
{

namespace detail
{

class one_val_matrix_store: public virtual_matrix_store
{
	const size_t mat_id;
	scalar_variable::ptr val;
	matrix_layout_t layout;
	std::vector<simple_raw_array> portion_bufs;
	int num_nodes;
	std::shared_ptr<NUMA_mapper> mapper;
public:
	one_val_matrix_store(scalar_variable::ptr val, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes);

	virtual bool has_materialized() const {
		// We never need to materialize this matrix.
		return false;
	}

	virtual size_t get_data_id() const {
		return INVALID_MAT_ID;
	}
	virtual bool share_data(const matrix_store &store) const;

	virtual std::string get_name() const {
		return (boost::format("one_val_mat(%1%,%2%)") % get_num_rows()
			% get_num_cols()).str();
	}
	virtual matrix_store::const_ptr materialize(bool in_mem,
			int num_nodes) const;
	virtual matrix_store::const_ptr get_cols(const std::vector<off_t> &idxs) const;
	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const;

	virtual std::pair<size_t, size_t> get_portion_size() const;

	using virtual_matrix_store::get_portion;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	using virtual_matrix_store::get_portion_async;
	virtual async_cres_t get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, std::shared_ptr<portion_compute> compute) const {
		return async_cres_t(true,
				get_portion(start_row, start_col, num_rows, num_cols));
	}
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t id) const;
	virtual int get_portion_node_id(size_t id) const;

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_layout_t store_layout() const {
		return layout;
	}
	virtual int get_num_nodes() const {
		return num_nodes;
	}

	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		std::unordered_map<size_t, size_t> ret;
		// TODO right now we only indicate the matrix. We set the number of
		// bytes to 0
		ret.insert(std::pair<size_t, size_t>(mat_id, 0));
		return ret;
	}
};

}

}

#endif
