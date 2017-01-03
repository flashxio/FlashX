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

#include "common.h"
#include "one_val_matrix_store.h"
#include "local_matrix_store.h"
#include "bulk_operate.h"
#include "vec_store.h"
#include "matrix_config.h"

namespace fm
{

namespace detail
{

one_val_matrix_store::one_val_matrix_store(scalar_variable::ptr val,
		size_t nrow, size_t ncol, matrix_layout_t layout,
		int num_nodes): virtual_matrix_store(nrow, ncol, true,
			val->get_type()), mat_id(mat_counter++)
{
	if (num_nodes > 0)
		this->mapper = std::shared_ptr<NUMA_mapper>(new NUMA_mapper(num_nodes,
					NUMA_range_size_log));
	this->val = val;
	this->layout = layout;
	this->num_nodes = num_nodes;

	size_t buf_size;
	if (is_wide()) {
		size_t part_num_cols = std::min(mem_matrix_store::CHUNK_SIZE,
				get_num_cols());
		buf_size = get_num_rows() * part_num_cols * get_entry_size();
	}
	else {
		size_t part_num_rows = std::min(mem_matrix_store::CHUNK_SIZE,
				get_num_rows());
		buf_size = get_num_cols() * part_num_rows * get_entry_size();
	}
	set_operate::const_ptr set = get_type().get_set_const(*val);
	if (num_nodes < 0) {
		this->portion_bufs.resize(1);
		this->portion_bufs[0] = simple_raw_array(buf_size, -1);
		set->set(portion_bufs[0].get_raw(), buf_size / get_entry_size(), 0, 0);
	}
	else {
		this->portion_bufs.resize(num_nodes);
		for (int i = 0; i < num_nodes; i++) {
			this->portion_bufs[i] = simple_raw_array(buf_size, i);
			set->set(portion_bufs[i].get_raw(), buf_size / get_entry_size(), 0, 0);
		}
	}
}

bool one_val_matrix_store::share_data(const matrix_store &store) const
{
	if ((store.get_num_rows() == get_num_rows()
				&& store.get_num_cols() == get_num_cols())
			|| (store.get_num_rows() == get_num_cols()
				&& store.get_num_cols() == get_num_rows())) {
		const one_val_matrix_store *mat1
			= dynamic_cast<const one_val_matrix_store *>(&store);
		if (mat1)
			return memcmp(mat1->val->get_raw(), val->get_raw(),
					mat1->val->get_size()) == 0;
		else
			return false;
	}
	else
		return false;
}

matrix_store::const_ptr one_val_matrix_store::materialize(bool in_mem,
			int num_nodes) const
{
	matrix_store::ptr ret = matrix_store::create(get_num_rows(),
			get_num_cols(), store_layout(), get_type(), num_nodes, in_mem);
	if (ret == NULL)
		return matrix_store::const_ptr();

	auto set = get_type().get_set_const(*val);
	ret->set_data(*set);
	return ret;
}

matrix_store::const_ptr one_val_matrix_store::get_cols(
		const std::vector<off_t> &idxs) const
{
	return matrix_store::const_ptr(new one_val_matrix_store(val,
				get_num_rows(), idxs.size(), layout, get_num_nodes()));
}

matrix_store::const_ptr one_val_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	return matrix_store::const_ptr(new one_val_matrix_store(val,
				idxs.size(), get_num_cols(), layout, get_num_nodes()));
}

int one_val_matrix_store::get_portion_node_id(size_t id) const
{
	if (get_num_nodes() < 0)
		return -1;

	if (is_wide()) {
		size_t start_col = mem_matrix_store::CHUNK_SIZE * id;
		return mapper->map2physical(start_col).first;
	}
	else {
		size_t start_row = mem_matrix_store::CHUNK_SIZE * id;
		return mapper->map2physical(start_row).first;
	}
}

local_matrix_store::const_ptr one_val_matrix_store::get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const
{
	int node_id = -1;
	const char *buf;
	if (get_num_nodes() > 0) {
		size_t range_size = mapper->get_range_size();
		if (is_wide()) {
			node_id = mapper->map2physical(start_col).first;
			assert(ROUND(start_col, range_size)
					== ROUND(start_col + num_cols - 1, range_size));
		}
		else {
			node_id = mapper->map2physical(start_row).first;
			assert(ROUND(start_row, range_size)
					== ROUND(start_row + num_rows - 1, range_size));
		}
		buf = portion_bufs[node_id].get_raw();
	}
	else
		buf = portion_bufs[0].get_raw();

	local_matrix_store::ptr ret;
	if (layout == matrix_layout_t::L_ROW) {
		ret = local_matrix_store::ptr(new local_cref_contig_row_matrix_store(
					buf, start_row, start_col, num_rows, num_cols, get_type(),
					node_id));
	}
	else if (layout == matrix_layout_t::L_COL) {
		ret = local_matrix_store::ptr(new local_cref_contig_col_matrix_store(
					buf, start_row, start_col, num_rows, num_cols, get_type(),
					node_id));
	}
	return ret;
}

local_matrix_store::const_ptr one_val_matrix_store::get_portion(size_t id) const
{
	size_t start_row;
	size_t start_col;
	size_t num_rows;
	size_t num_cols;
	if (is_wide()) {
		start_row = 0;
		start_col = mem_matrix_store::CHUNK_SIZE * id;
		num_rows = get_num_rows();
		num_cols = std::min(mem_matrix_store::CHUNK_SIZE,
				get_num_cols() - start_col);
	}
	else {
		start_row = mem_matrix_store::CHUNK_SIZE * id;
		start_col = 0;
		num_rows = std::min(mem_matrix_store::CHUNK_SIZE,
				get_num_rows() - start_row);
		num_cols = get_num_cols();
	}
	return get_portion(start_row, start_col, num_rows, num_cols);
}

matrix_store::const_ptr one_val_matrix_store::transpose() const
{
	matrix_layout_t new_layout;
	if (layout == matrix_layout_t::L_ROW)
		new_layout = matrix_layout_t::L_COL;
	else if (layout == matrix_layout_t::L_COL)
		new_layout = matrix_layout_t::L_ROW;
	else
		assert(0);
	return matrix_store::const_ptr(new one_val_matrix_store(val,
				get_num_cols(), get_num_rows(), new_layout, get_num_nodes()));
}

std::pair<size_t, size_t> one_val_matrix_store::get_portion_size() const
{
	if (is_wide())
		return std::pair<size_t, size_t>(get_num_rows(),
				mem_matrix_store::CHUNK_SIZE);
	else
		return std::pair<size_t, size_t>(mem_matrix_store::CHUNK_SIZE,
				get_num_cols());
}

}

}
