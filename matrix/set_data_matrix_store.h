#ifndef __SET_DATA_MATRIX_STORE_H__
#define __SET_DATA_MATRIX_STORE_H__

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

#include "matrix_config.h"
#include "bulk_operate.h"
#include "NUMA_mapper.h"
#include "virtual_matrix_store.h"

namespace fm
{

namespace detail
{

/*
 * This matrix assumes that set_operate can deterministically set the data
 * of the matrix and set_operate doesn't need to maintain a large volume of
 * data. This matrix is perfect for the matrices that has sequential numbers
 * or repeated data.
 */
class set_data_matrix_store: public virtual_matrix_store
{
	const size_t mat_id;
	const data_id_t::ptr data_id;
	set_operate::const_ptr row_op;
	set_operate::const_ptr col_op;
	matrix_layout_t layout;
	std::shared_ptr<NUMA_mapper> mapper;
	int num_nodes;

	set_data_matrix_store(set_operate::const_ptr row_op,
			set_operate::const_ptr col_op, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes): virtual_matrix_store(nrow,
				ncol, true, row_op->get_type()), mat_id(mat_counter++),
			data_id(data_id_t::create(mat_id)) {
		this->row_op = row_op;
		this->col_op = col_op;
		this->layout = layout;
		this->num_nodes = num_nodes;
		if (num_nodes > 0)
			this->mapper = std::shared_ptr<NUMA_mapper>(new NUMA_mapper(num_nodes,
						NUMA_range_size_log));
	}
	set_data_matrix_store(set_operate::const_ptr row_op,
			set_operate::const_ptr col_op, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes,
			data_id_t::ptr _data_id): virtual_matrix_store(nrow, ncol, true,
				row_op->get_type()), mat_id(mat_counter++), data_id(_data_id) {
		this->row_op = row_op;
		this->col_op = col_op;
		this->layout = layout;
		this->num_nodes = num_nodes;
		if (num_nodes > 0)
			this->mapper = std::shared_ptr<NUMA_mapper>(new NUMA_mapper(num_nodes,
						NUMA_range_size_log));
	}
public:
	static ptr create(set_operate::const_ptr row_op,
			set_operate::const_ptr col_op, size_t nrow, size_t ncol,
			matrix_layout_t layout, int num_nodes) {
		if (row_op->transpose() == NULL || col_op->transpose() == NULL) {
			BOOST_LOG_TRIVIAL(error) << "set_operate doesn't have transpose";
			return ptr();
		}
		return ptr(new set_data_matrix_store(row_op, col_op, nrow, ncol,
					layout, num_nodes));
	}

	virtual void inc_dag_ref(size_t id) {
		data_id->inc_ref(id);
	}
	virtual void reset_dag_ref() {
		data_id->reset_ref();
	}
	virtual size_t get_dag_ref() const {
		return data_id->get_ref();
	}

	virtual bool has_materialized() const {
		return false;
	}
	virtual size_t get_data_id() const {
		return data_id->get_id();
	}

	virtual std::string get_name() const {
		return (boost::format("set_data_mat(%1%,%2%)") % get_num_rows()
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
