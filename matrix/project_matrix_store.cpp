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
#include "dense_matrix.h"

namespace fm
{

namespace detail
{

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
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

}

sparse_project_matrix_store::sparse_project_matrix_store(size_t nrow,
		size_t ncol, matrix_layout_t layout,
		const scalar_type &type): virtual_matrix_store(nrow, ncol, true,
			type), mat_id(mat_counter++)
{
	this->layout = layout;
}

sparse_project_matrix_store::sparse_project_matrix_store(size_t nrow,
		size_t ncol, matrix_layout_t layout, const scalar_type &type,
		double density): virtual_matrix_store(nrow, ncol, true,
			type), mat_id(mat_counter++)
{
	this->layout = layout;

	// I intentially increase the density to generate the number of nnz.
	// It's possible I'll pick the same slot twice in the for loop below.
	size_t required = density * nrow * ncol;
	nz_idxs.resize((density * 1.1) * nrow * ncol);
	for (size_t i = 0; i < nz_idxs.size(); i++) {
		nz_idxs[i].row_idx = random() % nrow;
		nz_idxs[i].col_idx = random() % ncol;
	}
	if (layout == matrix_layout_t::L_ROW && is_wide())
		std::sort(nz_idxs.begin(), nz_idxs.end(), wide_row_first_comp());
	else if (is_wide())
		std::sort(nz_idxs.begin(), nz_idxs.end(), wide_col_first_comp());
	else if (layout == matrix_layout_t::L_ROW)
		std::sort(nz_idxs.begin(), nz_idxs.end(), tall_row_first_comp());
	else
		std::sort(nz_idxs.begin(), nz_idxs.end(), tall_col_first_comp());

	// There might be duplicates in the vector. We should remove them.
	auto end = std::unique(nz_idxs.begin(), nz_idxs.end());
	if (end - nz_idxs.begin() > (off_t) required)
		nz_idxs.resize(required);

	mem_col_matrix_store::ptr vals = mem_col_matrix_store::create(
			nz_idxs.size(), 1, type);
	this->vals = vals;

	size_t num_portions = get_num_portions();
	for (size_t i = 0; i < num_portions; i++) {
		// Get the starting and ending coordinates of a portion.
		sparse_project_matrix_store::nz_idx start;
		if (is_wide()) {
			start.row_idx = 0;
			start.col_idx = mem_matrix_store::CHUNK_SIZE * i;
		}
		else {
			start.row_idx = mem_matrix_store::CHUNK_SIZE * i;
			start.col_idx = 0;
		}

		// Get the starting and ending locations in the nz vector.
		std::vector<nz_idx>::const_iterator start_it;
		if (layout == matrix_layout_t::L_ROW && is_wide())
			start_it = std::lower_bound(nz_idxs.begin(), nz_idxs.end(), start,
					wide_row_first_comp());
		else if (is_wide())
			start_it = std::lower_bound(nz_idxs.begin(), nz_idxs.end(), start,
					wide_col_first_comp());
		else if (layout == matrix_layout_t::L_ROW)
			start_it = std::lower_bound(nz_idxs.begin(), nz_idxs.end(), start,
					tall_row_first_comp());
		else
			start_it = std::lower_bound(nz_idxs.begin(), nz_idxs.end(), start,
					tall_col_first_comp());

		// There aren't any data in the portion and the following portions.
		if (start_it == nz_idxs.end())
			break;

		portion_offs.push_back(start_it - nz_idxs.begin());
	}
	// We use the length of the nz vector to indicates the end of the last
	// portion.
	while (portion_offs.size() < num_portions + 1)
		portion_offs.push_back(nz_idxs.size());
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
				ncol, layout, type, density));
	if (type == get_scalar_type<double>())
		ret->vals->set_data(rand_init<double>());
	else if (type == get_scalar_type<int>())
		ret->vals->set_data(rand_init<int>());
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
	ret->portion_offs = this->portion_offs;
	ret->nz_idxs.resize(nz_idxs.size());
	for (size_t i = 0; i < nz_idxs.size(); i++) {
		ret->nz_idxs[i].row_idx = nz_idxs[i].col_idx;
		ret->nz_idxs[i].col_idx = nz_idxs[i].row_idx;
	}
	ret->vals = vals;
	return matrix_store::const_ptr(ret);
}

bool sparse_project_matrix_store::is_entire_portion(size_t start_row,
		size_t start_col, size_t num_rows, size_t num_cols) const
{
	if (is_wide())
		return start_row == 0 && num_rows == get_num_rows()
			&& start_col % mem_matrix_store::CHUNK_SIZE == 0
			&& (num_cols == mem_matrix_store::CHUNK_SIZE
					|| start_col + num_cols == get_num_cols());
	else
		return start_col == 0 && num_cols == get_num_cols()
			&& start_row % mem_matrix_store::CHUNK_SIZE == 0
			&& (num_rows == mem_matrix_store::CHUNK_SIZE
					|| start_row + num_rows == get_num_rows());
}

local_matrix_store::const_ptr sparse_project_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows, size_t num_cols) const
{
	// TODO we current only support get the entire portion.
	if (!is_entire_portion(start_row, start_col, num_rows, num_cols)) {
		BOOST_LOG_TRIVIAL(error) << "we only support getting the entire portion";
		return local_matrix_store::const_ptr();
	}

	if (start_row + num_rows > get_num_rows()
			|| start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "out of boundary";
		return local_matrix_store::const_ptr();
	}

	size_t portion_id;
	if (is_wide())
		portion_id = start_col / mem_matrix_store::CHUNK_SIZE;
	else
		portion_id = start_row / mem_matrix_store::CHUNK_SIZE;

	std::vector<nz_idx>::const_iterator start_it
		= nz_idxs.begin() + portion_offs[portion_id];
	std::vector<nz_idx>::const_iterator end_it
		= nz_idxs.begin() + portion_offs[portion_id + 1];

	// There isn't nnz in the portion.
	if (start_it == nz_idxs.end() || start_it == end_it) {
		std::vector<sparse_project_matrix_store::nz_idx> empty_idxs;
		local_matrix_store::ptr empty_vals(new local_buf_col_matrix_store(
					0, 0, 0, 0, vals->get_type(), -1));
		if (store_layout() == matrix_layout_t::L_COL)
			return local_matrix_store::const_ptr(new lsparse_col_matrix_store(
						start_row, start_col, num_rows, num_cols, empty_idxs,
						empty_vals));
		else
			return local_matrix_store::const_ptr(new lsparse_row_matrix_store(
						start_row, start_col, num_rows, num_cols, empty_idxs,
						empty_vals));
	}

	// Right now, the vector contains the global row and col indices.
	// But we want local indices. We'll convert them to local indices.
	std::vector<sparse_project_matrix_store::nz_idx> local_idxs(start_it,
			end_it);
	for (size_t i = 0; i < local_idxs.size(); i++) {
		assert(local_idxs[i].row_idx >= (off_t) start_row
				&& local_idxs[i].row_idx < (off_t) (start_row + num_rows));
		local_idxs[i].row_idx -= start_row;
		assert(local_idxs[i].col_idx >= (off_t) start_col
				&& local_idxs[i].col_idx < (off_t) (start_col + num_cols));
		local_idxs[i].col_idx -= start_col;
	}
	local_matrix_store::const_ptr local_vals = vals->get_portion(
			start_it - nz_idxs.begin(), 0, end_it - start_it, 1);
	assert(local_vals);

	if (store_layout() == matrix_layout_t::L_COL)
		return local_matrix_store::const_ptr(new lsparse_col_matrix_store(
					start_row, start_col, num_rows, num_cols, local_idxs,
					local_vals));
	else
		return local_matrix_store::const_ptr(new lsparse_row_matrix_store(
					start_row, start_col, num_rows, num_cols, local_idxs,
					local_vals));
}

matrix_store::const_ptr sparse_project_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	return matrix_store::const_ptr();
}

namespace
{

class conv_dense_op: public portion_mapply_op
{
public:
	conv_dense_op(size_t num_rows, size_t num_cols,
			const scalar_type &type): portion_mapply_op(num_rows,
				num_cols, type) {
	}

	virtual void run(const std::vector<local_matrix_store::const_ptr> &ins,
			local_matrix_store &out) const {
		assert(ins.size() == 1);
		// A sparse matrix might have no nnz at all.
		if (ins[0] == NULL)
			return;
		out.copy_from(*ins[0]);
	}

	virtual portion_mapply_op::const_ptr transpose() const {
		assert(0);
		return portion_mapply_op::const_ptr();
	}
	virtual std::string to_string(
			const std::vector<matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 1);
		return std::string("conv_dense(") + mats[0]->get_name() + ")";
	}
};

struct empty_deleter {
	void operator()(const matrix_store *addr) {
	}
};

}

matrix_store::const_ptr sparse_project_matrix_store::materialize(bool in_mem,
		int num_nodes) const
{
	std::vector<matrix_store::const_ptr> stores(1);
	stores[0] = matrix_store::const_ptr(this, empty_deleter());
	portion_mapply_op::const_ptr op(new conv_dense_op(get_num_rows(), get_num_cols(),
				get_type()));
	return __mapply_portion(stores, op, store_layout(), in_mem, num_nodes, true);
}

std::pair<size_t, size_t> sparse_project_matrix_store::get_portion_size() const
{
	if (is_wide())
		return std::pair<size_t, size_t>(get_num_rows(),
				mem_matrix_store::CHUNK_SIZE);
	else
		return std::pair<size_t, size_t>(mem_matrix_store::CHUNK_SIZE,
				get_num_cols());
}

matrix_store::const_ptr lsparse_col_matrix_store::transpose() const
{
	matrix_info info = get_global_transpose_info();
	std::vector<sparse_project_matrix_store::nz_idx> new_local_idxs(local_idxs.size());
	for (size_t i = 0; i < local_idxs.size(); i++) {
		new_local_idxs[i].row_idx = local_idxs[i].col_idx;
		new_local_idxs[i].col_idx = local_idxs[i].row_idx;
	}
	lsparse_row_matrix_store::ptr ret(new lsparse_row_matrix_store(
				info.start_row, info.start_col, info.num_rows, info.num_cols,
				new_local_idxs, vals));
	// If the matrix is smaller than its original size, we should resize it.
	if (get_num_rows() != info.num_cols || get_num_cols() != info.num_rows) {
		info = get_local_transpose_info();
		ret->resize(info.start_row, info.start_col, info.num_rows,
				info.num_cols);
	}
	return ret;
}

local_matrix_store::const_ptr lsparse_col_matrix_store::get_portion(
		size_t local_start_row, size_t local_start_col, size_t num_rows,
		size_t num_cols) const
{
	assert(0);
	return local_matrix_store::const_ptr();
}

const char *lsparse_col_matrix_store::get_col_nnz(off_t col_idx,
		std::vector<off_t> &row_idxs) const
{
	if (store_layout() == matrix_layout_t::L_ROW) {
		row_idxs.clear();
		return NULL;
	}

	sparse_project_matrix_store::nz_idx start;
	start.row_idx = get_local_start_row();
	start.col_idx = col_idx + get_local_start_col();
	auto it = std::lower_bound(local_idxs.begin(), local_idxs.end(), start,
			col_first_comp());
	if (it == local_idxs.end()) {
		row_idxs.clear();
		return NULL;
	}
	// There doesn't exist the col.
	if (it->col_idx != start.col_idx) {
		row_idxs.clear();
		return NULL;
	}

	// Find the right location in a col.
	if (it->row_idx < get_local_start_row())
		it++;
	if (it->row_idx >= (off_t) (get_local_start_row() + get_num_rows())) {
		row_idxs.clear();
		return NULL;
	}

	off_t val_off = it - local_idxs.begin();

	row_idxs.clear();
	for (; it != local_idxs.end() && it->col_idx == start.col_idx
			&& it->row_idx < (off_t) (get_local_start_row() + get_num_rows());
			it++) {
		row_idxs.push_back(it->row_idx - get_local_start_row());
	}
	return vals->get(val_off, 0);
}

void lsparse_col_matrix_store::materialize_self() const
{
	if (buf)
		return;

	local_col_matrix_store::ptr col_store(new local_buf_col_matrix_store(
				get_global_start_row(), get_global_start_col(),
				get_num_rows(), get_num_cols(), get_type(), -1));
	col_store->reset_data();
	const_cast<lsparse_col_matrix_store *>(this)->buf = col_store;

	std::vector<off_t> idxs;
	std::vector<char *> dst_ptrs;
	for (size_t i = 0; i < this->get_num_cols(); i++) {
		const char *in_col = this->get_col_nnz(i, idxs);
		// TODO this might be expensive. I should optimize it.
		char *out_col = col_store->get_col(i);
		dst_ptrs.resize(idxs.size());
		for (size_t j = 0; j < idxs.size(); j++)
			dst_ptrs[j] = out_col + idxs[j] * this->get_entry_size();
		this->get_type().get_sg().scatter(in_col, dst_ptrs);
	}
}

matrix_store::const_ptr lsparse_row_matrix_store::transpose() const
{
	matrix_info info = get_global_transpose_info();
	std::vector<sparse_project_matrix_store::nz_idx> new_local_idxs(local_idxs.size());
	for (size_t i = 0; i < local_idxs.size(); i++) {
		new_local_idxs[i].row_idx = local_idxs[i].col_idx;
		new_local_idxs[i].col_idx = local_idxs[i].row_idx;
	}
	lsparse_col_matrix_store::ptr ret(new lsparse_col_matrix_store(
				info.start_row, info.start_col, info.num_rows, info.num_cols,
				new_local_idxs, vals));
	// If the matrix is smaller than its original size, we should resize it.
	if (get_num_rows() != info.num_cols || get_num_cols() != info.num_rows) {
		info = get_local_transpose_info();
		ret->resize(info.start_row, info.start_col, info.num_rows,
				info.num_cols);
	}
	return ret;
}

local_matrix_store::const_ptr lsparse_row_matrix_store::get_portion(
		size_t local_start_row, size_t local_start_col, size_t num_rows,
		size_t num_cols) const
{
	assert(0);
	return local_matrix_store::const_ptr();
}

void lsparse_row_matrix_store::materialize_self() const
{
	if (buf)
		return;

	local_row_matrix_store::ptr row_store(new local_buf_row_matrix_store(
				get_global_start_row(), get_global_start_col(),
				get_num_rows(), get_num_cols(), get_type(), -1));
	row_store->reset_data();
	const_cast<lsparse_row_matrix_store *>(this)->buf = row_store;

	std::vector<off_t> idxs;
	std::vector<char *> dst_ptrs;
	for (size_t i = 0; i < this->get_num_rows(); i++) {
		const char *in_row = this->get_row_nnz(i, idxs);
		char *out_row = row_store->get_row(i);
		dst_ptrs.resize(idxs.size());
		for (size_t j = 0; j < idxs.size(); j++)
			dst_ptrs[j] = out_row + idxs[j] * this->get_entry_size();
		this->get_type().get_sg().scatter(in_row, dst_ptrs);
	}
}

const char *lsparse_row_matrix_store::get_row_nnz(off_t row_idx,
		std::vector<off_t> &col_idxs) const
{
	if (store_layout() == matrix_layout_t::L_COL) {
		col_idxs.clear();
		return NULL;
	}

	sparse_project_matrix_store::nz_idx start;
	start.row_idx = row_idx + get_local_start_row();
	start.col_idx = get_local_start_col();
	auto it = std::lower_bound(local_idxs.begin(), local_idxs.end(), start,
			row_first_comp());
	if (it == local_idxs.end()) {
		col_idxs.clear();
		return NULL;
	}
	// There doesn't exist the row.
	if (it->row_idx != start.row_idx) {
		col_idxs.clear();
		return NULL;
	}

	// Find the right location in a row.
	if (it->col_idx < get_local_start_col())
		it++;
	if (it->col_idx >= (off_t) (get_local_start_col() + get_num_cols())) {
		col_idxs.clear();
		return NULL;
	}

	off_t val_off = it - local_idxs.begin();

	col_idxs.clear();
	for (; it != local_idxs.end() && it->row_idx == start.row_idx
			&& it->col_idx < (off_t) (get_local_start_col() + get_num_cols());
			it++) {
		col_idxs.push_back(it->col_idx - get_local_start_col());
	}
	return vals->get(val_off, 0);
}

// Test if the current location is inside the specified range.
static inline bool inside_range(const sparse_project_matrix_store::nz_idx &start,
		const sparse_project_matrix_store::nz_idx &end,
		const sparse_project_matrix_store::nz_idx &curr)
{
	return curr.row_idx >= start.row_idx && curr.row_idx <= end.row_idx
		&& curr.col_idx >= start.col_idx && curr.col_idx < end.col_idx;
}

const char *lsparse_row_matrix_store::get_rows_nnz(off_t start_row,
		off_t end_row,
		std::vector<sparse_project_matrix_store::nz_idx> &idxs) const
{
	if (store_layout() == matrix_layout_t::L_COL) {
		idxs.clear();
		return NULL;
	}
	// If the portion is resized and rows aren't stored contiguously.
	if (get_orig_num_cols() != get_num_cols()) {
		idxs.clear();
		return NULL;
	}
	if (local_idxs.empty()) {
		idxs.clear();
		return NULL;
	}

	sparse_project_matrix_store::nz_idx start, end;
	start.row_idx = start_row + get_local_start_row();
	start.col_idx = get_local_start_col();
	end.row_idx = end_row - 1 + get_local_start_row();
	end.col_idx = get_local_start_col() + get_num_cols();
	auto start_it = std::lower_bound(local_idxs.begin(), local_idxs.end(),
			start, row_first_comp());
	// All non-empty rows are smaller than the starting row.
	if (start_it == local_idxs.end()) {
		idxs.clear();
		return NULL;
	}
	// The first row is out of the range.
	if (start_it->row_idx > end.row_idx) {
		idxs.clear();
		return NULL;
	}

	auto end_it = std::lower_bound(local_idxs.begin(), local_idxs.end(),
			end, row_first_comp());
	if (end_it != local_idxs.end()) {
		// If the last entry is inside the range.
		if (inside_range(start, end, *end_it))
			end_it++;
	}

	idxs.clear();
	for (auto it = start_it; it != end_it; it++) {
		// Rows have been resized. They can't be stored contiguously.
		if (it->col_idx >= (off_t) get_num_cols()) {
			idxs.clear();
			return NULL;
		}
		// The stored indexes are also relative.
		// The local store might have been resized.
		// The returned indexes should be relative to the current resized
		// local store.
		idxs.push_back(sparse_project_matrix_store::nz_idx(
					it->row_idx - get_local_start_row(),
					it->col_idx - get_local_start_col()));
		assert(inside_range(start, end, *it));
	}

	return vals->get(start_it - local_idxs.begin(), 0);
}

struct nz_loc
{
	sparse_project_matrix_store::nz_idx idx;
	off_t loc;
};

struct comp_row_first
{
	bool operator()(const nz_loc &e1, const nz_loc &e2) const {
		if (e1.idx.row_idx < e2.idx.row_idx)
			return true;
		else if (e1.idx.row_idx > e2.idx.row_idx)
			return false;
		else
			return e1.idx.col_idx < e2.idx.col_idx;
	}
};

struct comp_col_first
{
	bool operator()(const nz_loc &e1, const nz_loc &e2) const {
		if (e1.idx.col_idx < e2.idx.col_idx)
			return true;
		else if (e1.idx.col_idx > e2.idx.col_idx)
			return false;
		else
			return e1.idx.row_idx < e2.idx.row_idx;
	}
};

static void reshuffle(
		const std::vector<sparse_project_matrix_store::nz_idx> &local_idxs,
		local_matrix_store::const_ptr vals, matrix_layout_t layout,
		std::vector<sparse_project_matrix_store::nz_idx> &new_local_idxs,
		local_matrix_store::ptr new_vals)
{
	// sort the non-zero entries based on their row idxs.
	std::vector<nz_loc> nz_locs(local_idxs.size());
	for (size_t i = 0; i < nz_locs.size(); i++) {
		nz_locs[i].idx = local_idxs[i];
		nz_locs[i].loc = i;
	}
	if (layout == matrix_layout_t::L_ROW)
		std::sort(nz_locs.begin(), nz_locs.end(), comp_row_first());
	else
		std::sort(nz_locs.begin(), nz_locs.end(), comp_col_first());

	// Reshuffle data according to the row-major order.
	std::vector<char *> ptrs(local_idxs.size());
	for (size_t i = 0; i < nz_locs.size(); i++) {
		new_local_idxs[i] = nz_locs[i].idx;
		ptrs[i] = new_vals->get_raw_arr()
			+ nz_locs[i].loc * vals->get_type().get_size();
	}
	vals->get_type().get_sg().scatter(vals->get_raw_arr(), ptrs);
}

local_matrix_store::ptr lsparse_col_matrix_store::conv2(
		matrix_layout_t layout) const
{
	// If the layout is the same as the original one,
	// copy the original portion.
	if (layout == matrix_layout_t::L_COL)
		return local_matrix_store::ptr(new lsparse_col_matrix_store(*this));

	// We can't handle the resized portion.
	assert(get_orig_num_rows() == get_num_rows());
	assert(get_orig_num_cols() == get_num_cols());

	local_col_matrix_store::ptr new_vals(new local_buf_col_matrix_store(
				0, 0, vals->get_num_rows(), 1, vals->get_type(), -1));
	std::vector<sparse_project_matrix_store::nz_idx> new_local_idxs(
			local_idxs.size());
	reshuffle(local_idxs, vals, layout, new_local_idxs, new_vals);
	return local_matrix_store::ptr(new lsparse_row_matrix_store(
				get_global_start_col(), get_global_start_row(),
				get_num_cols(), get_num_rows(), new_local_idxs, new_vals));
}

local_matrix_store::ptr lsparse_row_matrix_store::conv2(
		matrix_layout_t layout) const
{
	// If the layout is the same as the original one,
	// copy the original portion.
	if (layout == matrix_layout_t::L_ROW)
		return local_matrix_store::ptr(new lsparse_row_matrix_store(*this));

	// We can't handle the resized portion.
	assert(get_orig_num_rows() == get_num_rows());
	assert(get_orig_num_cols() == get_num_cols());

	local_col_matrix_store::ptr new_vals(new local_buf_col_matrix_store(
				0, 0, vals->get_num_rows(), 1, vals->get_type(), -1));
	std::vector<sparse_project_matrix_store::nz_idx> new_local_idxs(
			local_idxs.size());
	reshuffle(local_idxs, vals, layout, new_local_idxs, new_vals);
	return local_matrix_store::ptr(new lsparse_row_matrix_store(
				get_global_start_col(), get_global_start_row(),
				get_num_cols(), get_num_rows(), new_local_idxs, new_vals));
}

}

}
