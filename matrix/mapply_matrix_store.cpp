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

#include "mapply_matrix_store.h"
#include "local_matrix_store.h"

namespace fm
{

namespace detail
{

namespace
{

class lmapply_col_matrix_store: public lvirtual_col_matrix_store
{
	std::vector<local_matrix_store::const_ptr> ins;
	const portion_mapply_op &op;
	raw_data_array data;

	void materialize() const {
		local_buf_col_matrix_store res(get_global_start_row(),
				get_global_start_col(), get_num_rows(), get_num_cols(),
				get_type(), -1);
		op.run(ins, res);
		lmapply_col_matrix_store *mutable_this
			= const_cast<lmapply_col_matrix_store *>(this);
		mutable_this->data = res.get_data();
	}
public:
	lmapply_col_matrix_store(
			const std::vector<local_matrix_store::const_ptr> &ins,
			const portion_mapply_op &_op,
			off_t global_start_row, off_t global_start_col,
			size_t nrow, size_t ncol, const scalar_type &type,
			int node_id): lvirtual_col_matrix_store(global_start_row,
				global_start_col, nrow, ncol, type, node_id), op(_op) {
		this->ins = ins;
	}

	virtual const char *get_raw_arr() const {
		if (data.is_empty())
			materialize();
		return data.get_raw();
	}

	virtual matrix_store::const_ptr transpose() const {
		// TODO
		assert(0);
		return matrix_store::const_ptr();
	}

	virtual const char *get_col(size_t col) const {
		if (data.is_empty())
			materialize();
		return data.get_raw() + get_num_rows() * col * get_type().get_size();
	}
};

class lmapply_row_matrix_store: public lvirtual_row_matrix_store
{
	std::vector<local_matrix_store::const_ptr> ins;
	const portion_mapply_op &op;
	raw_data_array data;

	void materialize() const {
		local_buf_row_matrix_store res(get_global_start_row(),
				get_global_start_col(), get_num_rows(), get_num_cols(),
				get_type(), -1);
		op.run(ins, res);
		lmapply_row_matrix_store *mutable_this
			= const_cast<lmapply_row_matrix_store *>(this);
		mutable_this->data = res.get_data();
	}
public:
	lmapply_row_matrix_store(
			const std::vector<local_matrix_store::const_ptr> &ins,
			const portion_mapply_op &_op,
			off_t global_start_row, off_t global_start_col,
			size_t nrow, size_t ncol, const scalar_type &type,
			int node_id): lvirtual_row_matrix_store(global_start_row,
				global_start_col, nrow, ncol, type, node_id), op(_op) {
		this->ins = ins;
	}

	virtual const char *get_raw_arr() const {
		if (data.is_empty())
			materialize();
		return data.get_raw();
	}

	virtual matrix_store::const_ptr transpose() const {
		// TODO
		assert(0);
		return matrix_store::const_ptr();
	}

	virtual const char *get_row(size_t row) const {
		if (data.is_empty())
			materialize();
		return data.get_raw() + get_num_cols() * row * get_type().get_size();
	}

	virtual const char *get_rows(size_t row_start, size_t row_end) const {
		// TODO
		assert(0);
		return NULL;
	}
};

}

mapply_matrix_store::mapply_matrix_store(
		const std::vector<mem_dense_matrix::const_ptr> &in_mats,
		portion_mapply_op::const_ptr op, matrix_layout_t layout,
		size_t nrow, size_t ncol): virtual_matrix_store(nrow, ncol,
			op->get_output_type())
{
	this->layout = layout;
	this->in_mats = in_mats;
	this->op = op;
}

matrix_store::ptr mapply_matrix_store::materialize() const
{
	return _mapply_portion(in_mats, op, layout);
}

matrix_store::const_ptr mapply_matrix_store::get_cols(
		const std::vector<off_t> &idxs) const
{
	// TODO
	assert(0);
	return matrix_store::const_ptr();
}

matrix_store::const_ptr mapply_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	// TODO
	assert(0);
	return matrix_store::const_ptr();
}

local_matrix_store::const_ptr mapply_matrix_store::get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const
{
	std::vector<local_matrix_store::const_ptr> parts(in_mats.size());
	if (is_wide()) {
		assert(start_row == 0);
		assert(num_rows == get_num_rows());
		for (size_t i = 0; i < in_mats.size(); i++)
			parts[i] = static_cast<const mem_matrix_store &>(
					in_mats[i]->get_data()).get_portion(start_row, start_col,
					in_mats[i]->get_num_rows(), num_cols);
	}
	else {
		assert(start_col == 0);
		assert(num_cols == get_num_cols());
		for (size_t i = 0; i < in_mats.size(); i++)
			parts[i] = static_cast<const mem_matrix_store &>(
					in_mats[i]->get_data()).get_portion(start_row, start_col,
					num_rows, in_mats[i]->get_num_cols());
	}

	if (store_layout() == matrix_layout_t::L_ROW)
		return local_matrix_store::const_ptr(new lmapply_row_matrix_store(
					parts, *op, start_row, start_col, num_rows, num_cols,
					get_type(), parts.front()->get_node_id()));
	else
		return local_matrix_store::const_ptr(new lmapply_col_matrix_store(
					parts, *op, start_row, start_col, num_rows, num_cols,
					get_type(), parts.front()->get_node_id()));
}

local_matrix_store::const_ptr mapply_matrix_store::get_portion(
			size_t id) const
{
	size_t start_row;
	size_t start_col;
	size_t num_rows;
	size_t num_cols;
	std::pair<size_t, size_t> chunk_size = get_portion_size();
	if (is_wide()) {
		start_row = 0;
		start_col = chunk_size.second * id;
		num_rows = get_num_rows();
		num_cols = std::min(chunk_size.second, get_num_cols() - start_col);
	}
	else {
		start_row = chunk_size.first * id;
		start_col = 0;
		num_rows = std::min(chunk_size.first, get_num_rows() - start_row);
		num_cols = get_num_cols();
	}
	return get_portion(start_row, start_col, num_rows, num_cols);
}

matrix_store::const_ptr mapply_matrix_store::transpose() const
{
	// TODO
	assert(0);
	return matrix_store::const_ptr();
}

}

}
