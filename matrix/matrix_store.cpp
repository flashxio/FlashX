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

#include "matrix_store.h"
#include "local_matrix_store.h"
#include "mem_matrix_store.h"
#include "EM_dense_matrix.h"

namespace fm
{

namespace detail
{

std::atomic<size_t> matrix_store::mat_counter;

matrix_store::ptr matrix_store::create(size_t nrow, size_t ncol,
		matrix_layout_t layout, const scalar_type &type, int num_nodes,
		bool in_mem, safs::safs_file_group::ptr group)
{
	if (in_mem)
		return mem_matrix_store::create(nrow, ncol, layout, type, num_nodes);
	else
		return EM_matrix_store::create(nrow, ncol, layout, type, group);
}

matrix_store::matrix_store(size_t nrow, size_t ncol, bool in_mem,
		const scalar_type &_type): type(_type)
{
	this->nrow = nrow;
	this->ncol = ncol;
	this->in_mem = in_mem;
	this->entry_size = type.get_size();
}

size_t matrix_store::get_num_portions() const
{
	std::pair<size_t, size_t> chunk_size = get_portion_size();
	if (is_wide())
		return ceil(((double) get_num_cols()) / chunk_size.second);
	else
		return ceil(((double) get_num_rows()) / chunk_size.first);
}

local_matrix_store::ptr matrix_store::get_portion(size_t id)
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

local_matrix_store::const_ptr matrix_store::get_portion(size_t id) const
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

}

}
