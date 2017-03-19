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

#include "set_data_matrix_store.h"
#include "mem_matrix_store.h"
#include "local_matrix_store.h"

namespace fm
{

namespace detail
{

matrix_store::const_ptr set_data_matrix_store::materialize(bool in_mem,
		int num_nodes) const
{
	matrix_store::ptr ret = matrix_store::create(get_num_rows(),
			get_num_cols(), store_layout(), get_type(), num_nodes, in_mem);
	if (ret == NULL)
		return matrix_store::const_ptr();

	if (store_layout() == matrix_layout_t::L_ROW)
		ret->set_data(*row_op);
	else
		ret->set_data(*col_op);
	return ret;
}

matrix_store::const_ptr set_data_matrix_store::get_cols(
		const std::vector<off_t> &idxs) const
{
	// We let the default solution to deal with this.
	// TODO we might need to handle this case better.
	if (!is_wide())
		return matrix_store::const_ptr();

	mem_matrix_store::ptr store = mem_matrix_store::create(get_num_rows(),
			idxs.size(), matrix_layout_t::L_COL, get_type(), -1);
	for (size_t i = 0; i < store->get_num_cols(); i++)
		col_op->set(store->get_col(i), store->get_num_rows(), 0, idxs[i]);
	return store;
}

matrix_store::const_ptr set_data_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	// We let the default solution to deal with this.
	// TODO we might need to handle this case better.
	if (is_wide())
		return matrix_store::const_ptr();

	mem_matrix_store::ptr store = mem_matrix_store::create(idxs.size(),
			get_num_cols(), matrix_layout_t::L_ROW, get_type(), -1);
	for (size_t i = 0; i < store->get_num_rows(); i++)
		row_op->set(store->get_row(i), store->get_num_cols(), idxs[i], 0);
	return store;
}

std::pair<size_t, size_t> set_data_matrix_store::get_portion_size() const
{
	if (is_wide())
		return std::pair<size_t, size_t>(get_num_rows(),
				mem_matrix_store::CHUNK_SIZE);
	else
		return std::pair<size_t, size_t>(mem_matrix_store::CHUNK_SIZE,
				get_num_cols());
}

local_matrix_store::const_ptr set_data_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	int node_id;
	if (get_num_nodes() < 0)
		node_id = -1;
	else if (is_wide())
		node_id = mapper->map2physical(start_col).first;
	else
		node_id = mapper->map2physical(start_row).first;

	if (store_layout() == matrix_layout_t::L_ROW) {
		local_row_matrix_store::ptr buf(new local_buf_row_matrix_store(
					start_row, start_col, num_rows, num_cols, get_type(),
					node_id));
		for (size_t i = 0; i < buf->get_num_rows(); i++)
			row_op->set(buf->get_row(i), buf->get_num_cols(), start_row + i,
					start_col);
		return buf;
	}
	else {
		local_col_matrix_store::ptr buf(new local_buf_col_matrix_store(
					start_row, start_col, num_rows, num_cols, get_type(),
					node_id));
		for (size_t i = 0; i < buf->get_num_cols(); i++)
			col_op->set(buf->get_col(i), buf->get_num_rows(), start_row,
					start_col + i);
		return buf;
	}
}

local_matrix_store::const_ptr set_data_matrix_store::get_portion(size_t id) const
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

int set_data_matrix_store::get_portion_node_id(size_t id) const
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

matrix_store::const_ptr set_data_matrix_store::transpose() const
{
	matrix_layout_t new_layout = layout
		== matrix_layout_t::L_ROW ? matrix_layout_t::L_COL : matrix_layout_t::L_ROW;
	return matrix_store::const_ptr(new set_data_matrix_store(col_op->transpose(),
				row_op->transpose(), get_num_cols(), get_num_rows(), new_layout,
				get_num_nodes(), data_id));
}

}

}
