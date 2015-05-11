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

#include "one_val_matrix_store.h"
#include "local_matrix_store.h"
#include "bulk_operate.h"

namespace fm
{

namespace detail
{

one_val_matrix_store::one_val_matrix_store(scalar_variable::ptr val,
		size_t nrow, size_t ncol, matrix_layout_t layout): virtual_matrix_store(
			nrow, ncol, val->get_type())
{
	this->val = val;
	this->layout = layout;

	size_t buf_size;
	if (is_wide()) {
		size_t part_num_cols = std::min(CHUNK_SIZE, get_num_cols());
		buf_size = get_num_rows() * part_num_cols * get_entry_size();
	}
	else {
		size_t part_num_rows = std::min(CHUNK_SIZE, get_num_rows());
		buf_size = get_num_cols() * part_num_rows * get_entry_size();
	}
	this->portion_buf = std::unique_ptr<char[]>(new char[buf_size]);
	const set_operate &set = get_type().get_set_const(*val);
	set.set(portion_buf.get(), buf_size / get_entry_size(), 0, 0);
}

matrix_store::const_ptr one_val_matrix_store::get_cols(
		const std::vector<off_t> &idxs) const
{
	return matrix_store::const_ptr(new one_val_matrix_store(val,
				get_num_rows(), idxs.size(), layout));
}

matrix_store::const_ptr one_val_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	return matrix_store::const_ptr(new one_val_matrix_store(val,
				idxs.size(), get_num_cols(), layout));
}

local_matrix_store::const_ptr one_val_matrix_store::get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const
{
	local_matrix_store::ptr ret;
	if (layout == matrix_layout_t::L_ROW) {
		ret = local_matrix_store::ptr(new local_cref_contig_row_matrix_store(
					start_row, start_col, portion_buf.get(),
					num_rows, num_cols, get_type(), -1));
	}
	else if (layout == matrix_layout_t::L_COL) {
		ret = local_matrix_store::ptr(new local_cref_contig_col_matrix_store(
					start_row, start_col, portion_buf.get(),
					num_rows, num_cols, get_type(), -1));
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
		start_col = CHUNK_SIZE * id;
		num_rows = get_num_rows();
		num_cols = std::min(CHUNK_SIZE, get_num_cols() - start_col);
	}
	else {
		start_row = CHUNK_SIZE * id;
		start_col = 0;
		num_rows = std::min(CHUNK_SIZE, get_num_rows() - start_row);
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
				get_num_cols(), get_num_rows(), new_layout));
}

}

}
