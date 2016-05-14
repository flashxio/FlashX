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

#include "project_matrix_store.h"

namespace fm
{

namespace detail
{

sparse_project_matrix_store::sparse_project_matrix_store(size_t nrow,
		size_t ncol, matrix_layout_t layout, const scalar_type &type): mem_matrix_store(
			nrow, ncol, type)
{
	this->layout = layout;
}

namespace
{

template<class T>
class rand_init: public type_set_operate<T>
{
public:
	virtual void set(T *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			if (random() % 2 == 0)
				arr[i] = 1;
			else
				arr[i] = -1;
		}
	}
};

}

sparse_project_matrix_store::ptr sparse_project_matrix_store::create_sparse_rand(
		size_t nrow, size_t ncol, matrix_layout_t layout,
		const scalar_type &type, double density)
{
	if (type != get_scalar_type<double>() && type != get_scalar_type<int>()) {
		BOOST_LOG_TRIVIAL(error)
			<< "sparse project matrix only supports int and double";
		return sparse_project_matrix_store::ptr();
	}

	sparse_project_matrix_store::ptr ret(new sparse_project_matrix_store(nrow,
				ncol, layout, type));
	// I intentially increase the density to generate the number of nnz.
	// It's possible I'll pick the same slot twice in the for loop below.
	size_t required = density * nrow * ncol;
	ret->nnz_idxs.resize((density * 1.1) * nrow * ncol);
	for (size_t i = 0; i < ret->nnz_idxs.size(); i++) {
		ret->nnz_idxs[i].row_idx = random() % nrow;
		ret->nnz_idxs[i].col_idx = random() % ncol;
	}
	if (layout == matrix_layout_t::L_ROW && ret->is_wide())
		std::sort(ret->nnz_idxs.begin(), ret->nnz_idxs.end(),
				wide_row_first_comp());
	else if (ret->is_wide())
		std::sort(ret->nnz_idxs.begin(), ret->nnz_idxs.end(),
				wide_col_first_comp());
	else if (layout == matrix_layout_t::L_ROW)
		std::sort(ret->nnz_idxs.begin(), ret->nnz_idxs.end(),
				tall_row_first_comp());
	else
		std::sort(ret->nnz_idxs.begin(), ret->nnz_idxs.end(),
				tall_col_first_comp());

	// There might be duplicates in the vector. We should remove them.
	auto end = std::unique(ret->nnz_idxs.begin(), ret->nnz_idxs.end());
	if (end - ret->nnz_idxs.begin() > (off_t) required)
		ret->nnz_idxs.resize(required);

	mem_col_matrix_store::ptr vals = mem_col_matrix_store::create(
			ret->nnz_idxs.size(), 1, type);
	if (type == get_scalar_type<double>())
		vals->set_data(rand_init<double>());
	else if (type == get_scalar_type<int>())
		vals->set_data(rand_init<int>());
	ret->vals = vals;
	return ret;
}

const char *sparse_project_matrix_store::get(size_t row, size_t col) const
{
	throw unsupported_exception(
			"don't support getting an entry in a sparse matrix");
}

matrix_store::const_ptr sparse_project_matrix_store::transpose() const
{
	matrix_layout_t new_layout;
	if (store_layout() == matrix_layout_t::L_COL)
		new_layout = matrix_layout_t::L_ROW;
	else
		new_layout = matrix_layout_t::L_COL;

	sparse_project_matrix_store *ret = new sparse_project_matrix_store(
			get_num_cols(), get_num_rows(), new_layout, get_type());
	ret->nnz_idxs.resize(nnz_idxs.size());
	for (size_t i = 0; i < nnz_idxs.size(); i++) {
		ret->nnz_idxs[i].row_idx = nnz_idxs[i].col_idx;
		ret->nnz_idxs[i].col_idx = nnz_idxs[i].row_idx;
	}
	ret->vals = vals;
	return matrix_store::const_ptr(ret);
}

local_matrix_store::const_ptr sparse_project_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows, size_t num_cols) const
{
	// TODO we current only support get the entire portion.
	if (is_wide() && (start_row > 0 || num_rows != get_num_rows())) {
		BOOST_LOG_TRIVIAL(error) << "we only support getting the entire portion";
		return local_matrix_store::const_ptr();
	}
	else if (!is_wide() && (start_col > 0 || num_cols != get_num_cols())) {
		BOOST_LOG_TRIVIAL(error) << "we only support getting the entire portion";
		return local_matrix_store::const_ptr();
	}

	// The requested portion is in the same physical portion.
	if (is_wide() && start_col / mem_matrix_store::CHUNK_SIZE
			!= (start_col + num_cols - 1) / mem_matrix_store::CHUNK_SIZE) {
		BOOST_LOG_TRIVIAL(error)
			<< "the requested portion should be in the same physical portion";
		return local_matrix_store::const_ptr();
	}
	else if (!is_wide() && start_row / mem_matrix_store::CHUNK_SIZE
			!= (start_row + num_rows - 1) / mem_matrix_store::CHUNK_SIZE) {
		BOOST_LOG_TRIVIAL(error)
			<< "the requested portion should be in the same physical portion";
		return local_matrix_store::const_ptr();
	}

	if (start_row + num_rows > get_num_rows()
			|| start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "out of boundary";
		return local_matrix_store::const_ptr();
	}

	// Find the starting and ending points in the vector of offsets.
	// Search for [start, end]
	sparse_project_matrix_store::nnz_idx start, end;
	start.row_idx = start_row;
	start.col_idx = start_col;
	end.row_idx = start_row + num_rows - 1;
	end.col_idx = start_col + num_cols - 1;
	std::vector<nnz_idx>::const_iterator start_it, end_it;
	if (layout == matrix_layout_t::L_ROW && is_wide()) {
		start_it = std::lower_bound(nnz_idxs.begin(), nnz_idxs.end(), start,
				wide_row_first_comp());
		end_it = std::lower_bound(nnz_idxs.begin(), nnz_idxs.end(), end,
				wide_row_first_comp());
	}
	else if (is_wide()) {
		start_it = std::lower_bound(nnz_idxs.begin(), nnz_idxs.end(), start,
				wide_col_first_comp());
		end_it = std::lower_bound(nnz_idxs.begin(), nnz_idxs.end(), end,
				wide_col_first_comp());
	}
	else if (layout == matrix_layout_t::L_ROW) {
		start_it = std::lower_bound(nnz_idxs.begin(), nnz_idxs.end(), start,
				tall_row_first_comp());
		end_it = std::lower_bound(nnz_idxs.begin(), nnz_idxs.end(), end,
				tall_row_first_comp());
	}
	else {
		start_it = std::lower_bound(nnz_idxs.begin(), nnz_idxs.end(), start,
				tall_col_first_comp());
		end_it = std::lower_bound(nnz_idxs.begin(), nnz_idxs.end(), end,
				tall_col_first_comp());
	}

	// There isn't nnz in the portion.
	if (start_it == nnz_idxs.end() || start_it == end_it)
		return local_matrix_store::const_ptr();

	// We search for [start, end].
	// If the last element has value `end', we need to move the end iterator
	// to include the last element.
	if (end_it != nnz_idxs.end()
			&& end_it->row_idx == end.row_idx && end_it->col_idx == end.col_idx)
		end_it++;
	std::vector<sparse_project_matrix_store::nnz_idx> local_idxs(start_it,
			end_it);
	for (size_t i = 0; i < local_idxs.size(); i++) {
		assert(local_idxs[i].row_idx >= start_row
				&& local_idxs[i].row_idx < start_row + num_rows);
		local_idxs[i].row_idx -= start_row;
		assert(local_idxs[i].col_idx >= start_col
				&& local_idxs[i].col_idx < start_col + num_cols);
		local_idxs[i].col_idx -= start_col;
	}
	local_matrix_store::const_ptr local_vals = vals->get_portion(
			start_it - nnz_idxs.begin(), 0, end_it - start_it, 1);
	assert(local_vals);
	return local_matrix_store::const_ptr(new local_sparse_matrix_store(
				start_row, start_col, num_rows, num_cols, local_idxs,
				local_vals, store_layout()));
}

void sparse_project_matrix_store::write_portion_async(
		local_matrix_store::const_ptr portion, off_t start_row, off_t start_col)
{
	throw unsupported_exception("don't support write_portion_async");
}

bool local_sparse_matrix_store::resize(off_t local_start_row,
		off_t local_start_col, size_t local_num_rows, size_t local_num_cols)
{
	throw unsupported_exception("don't support resizing a local sparse matrix");
}

local_matrix_store::ptr local_sparse_matrix_store::conv2(
		matrix_layout_t layout) const
{
	throw unsupported_exception(
			"don't support conv layout of a local sparse matrix");
}

size_t local_sparse_matrix_store::get_all_rows(std::vector<const char *> &rows) const
{
	throw unsupported_exception(
			"don't support getting all rows of a local sparse matrix");
}

size_t local_sparse_matrix_store::get_all_cols(std::vector<const char *> &cols) const
{
	throw unsupported_exception(
			"don't support getting all cols of a local sparse matrix");
}

local_matrix_store::const_ptr local_sparse_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	throw unsupported_exception(
			"don't support getting a portion from a local sparse matrix");
}

const char *local_sparse_matrix_store::get(size_t row, size_t col) const
{
	return NULL;
}

matrix_store::const_ptr local_sparse_matrix_store::transpose() const
{
	return matrix_store::const_ptr();
}

local_sparse_matrix_store::local_sparse_matrix_store(off_t global_start_row,
		off_t global_start_col, size_t nrow, size_t ncol,
		const std::vector<sparse_project_matrix_store::nnz_idx> &local_idxs,
		local_matrix_store::const_ptr vals, matrix_layout_t layout): local_matrix_store(
			global_start_row, global_start_col, nrow, ncol, vals->get_type(),
			vals->get_node_id())
{
	this->local_idxs = local_idxs;
	this->vals = vals;
	this->layout = layout;
}

const char *local_sparse_matrix_store::get_col_nnz(off_t col_idx,
		std::vector<off_t> &row_idxs) const
{
	if (store_layout() == matrix_layout_t::L_ROW) {
		row_idxs.clear();
		return NULL;
	}

	sparse_project_matrix_store::nnz_idx start;
	start.row_idx = 0;
	start.col_idx = col_idx;
	auto it = std::lower_bound(local_idxs.begin(), local_idxs.end(), start,
			col_comp());
	if (it == local_idxs.end()) {
		row_idxs.clear();
		return NULL;
	}
	// There doesn't exist the col.
	if (it->col_idx != col_idx) {
		row_idxs.clear();
		return NULL;
	}

	off_t val_off = it - local_idxs.begin();

	row_idxs.clear();
	for (; it != local_idxs.end() && it->col_idx == col_idx; it++) {
		row_idxs.push_back(it->row_idx);
	}
	return vals->get(val_off, 0);
}

const char *local_sparse_matrix_store::get_row_nnz(off_t row_idx,
		std::vector<off_t> &col_idxs) const
{
	if (store_layout() == matrix_layout_t::L_COL) {
		col_idxs.clear();
		return NULL;
	}

	sparse_project_matrix_store::nnz_idx start;
	start.row_idx = row_idx;
	start.col_idx = 0;
	auto it = std::lower_bound(local_idxs.begin(), local_idxs.end(), start,
			row_comp());
	if (it == local_idxs.end()) {
		col_idxs.clear();
		return NULL;
	}
	// There doesn't exist the row.
	if (it->row_idx != row_idx) {
		col_idxs.clear();
		return NULL;
	}
	off_t val_off = it - local_idxs.begin();

	col_idxs.clear();
	for (; it != local_idxs.end() && it->row_idx == row_idx; it++) {
		col_idxs.push_back(it->col_idx);
	}
	return vals->get(val_off, 0);
}

}

}
