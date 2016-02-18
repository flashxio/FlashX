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

#include <boost/format.hpp>

#include <immintrin.h>

#include <cblas.h>

#include "local_matrix_store.h"
#include "bulk_operate.h"
#include "dense_matrix.h"
#include "local_vec_store.h"
#include "factor.h"

namespace fm
{

namespace detail
{

static const size_t LONG_DIM_LEN = 1024;

#ifdef __AVX__
static void memcpy256(char *dest, const char *src, size_t num256)
{
	__m256i *dst_arr = reinterpret_cast<__m256i *>(dest);
	const __m256i *src_arr = reinterpret_cast<const __m256i *>(src);
	for (size_t i = 0; i < num256; i++)
		_mm256_stream_si256(dst_arr + i, *(src_arr + i));
}
#endif

bool local_matrix_store::large_copy_from(const local_matrix_store &store)
{
#ifdef __AVX__
	char *dst_arr = get_raw_arr();
	const char *src_arr = get_raw_arr();
	size_t num_bytes = get_num_rows() * get_num_cols() * get_entry_size();
	if (dst_arr && src_arr && num_bytes % sizeof(__m256i) == 0) {
		memcpy256(dst_arr, src_arr, num_bytes / sizeof(__m256i));
		return true;
	}
	else
#endif
		return false;
}

local_matrix_store::ptr local_matrix_store::conv2(matrix_layout_t layout) const
{
	local_matrix_store::ptr ret;
	if (layout == matrix_layout_t::L_ROW)
		ret = local_matrix_store::ptr(new local_buf_row_matrix_store(
					get_global_start_row(), get_global_start_col(),
					get_num_rows(), get_num_cols(), get_type(), get_node_id()));
	else
		ret = local_matrix_store::ptr(new local_buf_col_matrix_store(
					get_global_start_row(), get_global_start_col(),
					get_num_rows(), get_num_cols(), get_type(), get_node_id()));
	ret->copy_from(*this);
	return ret;
}

bool local_matrix_store::resize(off_t local_start_row,
		off_t local_start_col, size_t local_num_rows, size_t local_num_cols)
{
	assert(local_start_row + local_num_rows <= get_orig_num_rows()
			|| local_start_col + local_num_cols <= get_orig_num_cols());
	this->local_start_row = local_start_row;
	this->local_start_col = local_start_col;
	matrix_store::resize(local_num_rows, local_num_cols);
	return false;
}

void local_matrix_store::resize_transpose(local_matrix_store &store) const
{
	if (get_local_start_row() > 0 || get_local_start_col() > 0
			|| get_num_rows() < get_orig_num_rows()
			|| get_num_cols() < get_orig_num_cols()) {
		struct matrix_info t_info = get_local_transpose_info();
		store.resize(t_info.start_row, t_info.start_col, t_info.num_rows,
				t_info.num_cols);
	}
}

matrix_store::const_ptr local_buf_col_matrix_store::transpose() const
{
	return create_transpose<local_buf_row_matrix_store,
		   local_buf_col_matrix_store>(*this);
}

matrix_store::ptr local_buf_col_matrix_store::transpose()
{
	return create_transpose<local_buf_row_matrix_store,
		   local_buf_col_matrix_store>(*this);
}

matrix_store::const_ptr local_buf_row_matrix_store::transpose() const
{
	return create_transpose<local_buf_col_matrix_store,
		   local_buf_row_matrix_store>(*this);
}

matrix_store::ptr local_buf_row_matrix_store::transpose()
{
	return create_transpose<local_buf_col_matrix_store,
		   local_buf_row_matrix_store>(*this);
}

matrix_store::const_ptr local_ref_contig_col_matrix_store::transpose() const
{
	return create_transpose<local_cref_contig_row_matrix_store,
		   local_ref_contig_col_matrix_store>(*this);
}

matrix_store::ptr local_ref_contig_col_matrix_store::transpose()
{
	return create_transpose<local_ref_contig_row_matrix_store,
		   local_ref_contig_col_matrix_store>(*this);
}

matrix_store::const_ptr local_ref_contig_row_matrix_store::transpose() const
{
	return create_transpose<local_cref_contig_col_matrix_store,
		   local_ref_contig_row_matrix_store>(*this);
}

matrix_store::ptr local_ref_contig_row_matrix_store::transpose()
{
	return create_transpose<local_ref_contig_col_matrix_store,
		   local_ref_contig_row_matrix_store>(*this);
}

matrix_store::const_ptr local_ref_col_matrix_store::transpose() const
{
	return create_transpose<local_cref_row_matrix_store,
		   local_ref_col_matrix_store>(*this);
}

matrix_store::ptr local_ref_col_matrix_store::transpose()
{
	return create_transpose<local_ref_row_matrix_store,
		   local_ref_col_matrix_store>(*this);
}

matrix_store::const_ptr local_ref_row_matrix_store::transpose() const
{
	return create_transpose<local_cref_col_matrix_store,
		   local_ref_row_matrix_store>(*this);
}

matrix_store::ptr local_ref_row_matrix_store::transpose()
{
	return create_transpose<local_ref_col_matrix_store,
		   local_ref_row_matrix_store>(*this);
}

matrix_store::const_ptr local_cref_contig_col_matrix_store::transpose() const
{
	return create_transpose<local_cref_contig_row_matrix_store,
		   local_cref_contig_col_matrix_store>(*this);
}

matrix_store::ptr local_cref_contig_col_matrix_store::transpose()
{
	return create_transpose<local_cref_contig_row_matrix_store,
		   local_cref_contig_col_matrix_store>(*this);
}

matrix_store::const_ptr local_cref_contig_row_matrix_store::transpose() const
{
	return create_transpose<local_cref_contig_col_matrix_store,
		   local_cref_contig_row_matrix_store>(*this);
}

matrix_store::ptr local_cref_contig_row_matrix_store::transpose()
{
	return create_transpose<local_cref_contig_col_matrix_store,
		   local_cref_contig_row_matrix_store>(*this);
}

matrix_store::const_ptr local_cref_col_matrix_store::transpose() const
{
	return create_transpose<local_cref_row_matrix_store,
		   local_cref_col_matrix_store>(*this);
}

matrix_store::ptr local_cref_col_matrix_store::transpose()
{
	return create_transpose<local_cref_row_matrix_store,
		   local_cref_col_matrix_store>(*this);
}

matrix_store::const_ptr local_cref_row_matrix_store::transpose() const
{
	return create_transpose<local_cref_col_matrix_store,
		   local_cref_row_matrix_store>(*this);
}

matrix_store::ptr local_cref_row_matrix_store::transpose()
{
	return create_transpose<local_cref_col_matrix_store,
		   local_cref_row_matrix_store>(*this);
}

void local_row_matrix_store::reset_data()
{
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	// If the store has data stored contiguously.
	char *raw_arr = get_raw_arr();
	if (raw_arr)
		memset(raw_arr, 0, nrow * ncol * get_entry_size());
	else {
		for (size_t i = 0; i < nrow; i++)
			memset(get_row(i), 0, ncol * get_entry_size());
	}
}

void local_row_matrix_store::set_data(const set_operate &op)
{
	assert(op.get_type() == get_type());
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	char *raw_arr = get_raw_arr();
	if (raw_arr) {
		// If data is stored in contiguous memory, we can calculate
		// the address of each row directly.
		for (size_t i = 0; i < nrow; i++) {
			char *row = raw_arr + ncol * i * get_entry_size();
			op.set(row, ncol, get_global_start_row() + i,
					get_global_start_col());
		}
	}
	else {
		for (size_t i = 0; i < nrow; i++)
			op.set(get_row(i), ncol, get_global_start_row() + i,
					get_global_start_col());
	}
}

bool local_row_matrix_store::copy_from(const local_matrix_store &store)
{
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	assert(nrow == store.get_num_rows());
	assert(ncol == store.get_num_cols());
	assert(store.get_type() == this->get_type());

	char *this_arr = get_raw_arr();
	const char *other_arr = store.get_raw_arr();
	// If this is a one-row or one col matrix, and data is stored contiguously,
	// we only need to copy data.
	if ((nrow == 1 || ncol == 1) && (this_arr && other_arr))
		memcpy(this_arr, other_arr, nrow * ncol * get_entry_size());
	else if (store.store_layout() == matrix_layout_t::L_ROW) {
		// If the store has data stored contiguously.
		if (this_arr && other_arr)
			memcpy(this_arr, other_arr, nrow * ncol * get_entry_size());
		else {
			const local_row_matrix_store &row_store
				= static_cast<const local_row_matrix_store &>(store);
			for (size_t i = 0; i < nrow; i++)
				memcpy(get_row(i), row_store.get_row(i),
						ncol * get_entry_size());
		}
	}
	// We should handle wide matrix and tall matrix differently
	// to minimize virtual function calls.
	// The idea is to reduce cache misses.
	else if (is_wide() && store.get_raw_arr() == NULL) {
		const local_col_matrix_store &col_store
			= static_cast<const local_col_matrix_store &>(store);
		std::vector<char *> dest_col(get_num_rows());
		size_t entry_size = store.get_type().get_size();
		const scatter_gather &sg = store.get_type().get_sg();
		for (size_t j = 0; j < store.get_num_rows(); j++)
			dest_col[j] = this->get_row(j);
		// This copies from column by column.
		for (size_t i = 0; i < store.get_num_cols(); i++) {
			const char *src_col = col_store.get_col(i);
			sg.scatter(src_col, dest_col);
			for (size_t j = 0; j < store.get_num_rows(); j++)
				dest_col[j] += entry_size;
		}
	}
	else if (is_wide()) {
		std::vector<char *> dst_rows(get_num_rows());
		for (size_t i = 0; i < get_num_rows(); i++)
			dst_rows[i] = get_row(i);
		get_type().get_conv().conv(store.get_raw_arr(),
				get_num_rows() * get_num_cols(), dst_rows);
	}
	else if (get_raw_arr() == NULL) {
		const local_col_matrix_store &col_store
			= static_cast<const local_col_matrix_store &>(store);
		std::vector<const char *> src_row(get_num_cols());
		size_t entry_size = store.get_type().get_size();
		const scatter_gather &sg = store.get_type().get_sg();
		for (size_t j = 0; j < store.get_num_cols(); j++)
			src_row[j] = col_store.get_col(j);
		// This copies from row by row.
		for (size_t i = 0; i < store.get_num_rows(); i++) {
			char *dst_row = get_row(i);
			sg.gather(src_row, dst_row);
			for (size_t j = 0; j < src_row.size(); j++)
				src_row[j] += entry_size;
		}
	}
	else {
		std::vector<const char *> src_cols(get_num_cols());
		const local_col_matrix_store &col_store
			= static_cast<const local_col_matrix_store &>(store);
		for (size_t i = 0; i < get_num_cols(); i++)
			src_cols[i] = col_store.get_col(i);
		get_type().get_conv().conv(src_cols, get_num_rows(), get_raw_arr());
	}
	return true;
}

void local_col_matrix_store::reset_data()
{
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	// If the store has data stored contiguously.
	char *raw_arr = get_raw_arr();
	if (raw_arr)
		memset(raw_arr, 0, nrow * ncol * get_entry_size());
	else {
		for (size_t i = 0; i < ncol; i++)
			memset(get_col(i), 0, nrow * get_entry_size());
	}
}

void local_col_matrix_store::set_data(const set_operate &op)
{
	assert(op.get_type() == get_type());
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	char *raw_arr = get_raw_arr();
	if (raw_arr) {
		// If data is stored in contiguous memory, we can calculate
		// the address of each row directly.
		for (size_t i = 0; i < ncol; i++) {
			char *col = raw_arr + nrow * i * get_entry_size();
			op.set(col, nrow, get_global_start_row(),
					get_global_start_col() + i);
		}
	}
	else {
		for (size_t i = 0; i < ncol; i++)
			op.set(get_col(i), nrow, get_global_start_row(),
					get_global_start_col() + i);
	}
}

bool local_col_matrix_store::copy_from(const local_matrix_store &store)
{
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	assert(nrow == store.get_num_rows());
	assert(ncol == store.get_num_cols());
	assert(store.get_type() == this->get_type());

	char *this_arr = get_raw_arr();
	const char *other_arr = store.get_raw_arr();
	// If this is a one-row or one col matrix, and data is stored contiguously,
	// we only need to copy data.
	if ((nrow == 1 || ncol == 1) && (this_arr && other_arr))
		memcpy(this_arr, other_arr, nrow * ncol * get_entry_size());
	else if (store.store_layout() == matrix_layout_t::L_COL) {
		// If the store has data stored contiguously.
		if (this_arr && other_arr)
			memcpy(this_arr, other_arr, nrow * ncol * get_entry_size());
		else {
			const local_col_matrix_store &col_store
				= static_cast<const local_col_matrix_store &>(store);
			for (size_t i = 0; i < ncol; i++)
				memcpy(get_col(i), col_store.get_col(i), nrow * get_entry_size());
		}
	}
	// We should handle wide matrix and tall matrix differently
	// to minimize virtual function calls.
	else if (is_wide() && get_raw_arr() == NULL) {
		const local_row_matrix_store &row_store
			= static_cast<const local_row_matrix_store &>(store);
		std::vector<const char *> src_col(get_num_rows());
		size_t entry_size = store.get_type().get_size();
		const scatter_gather &sg = store.get_type().get_sg();
		for (size_t j = 0; j < store.get_num_rows(); j++)
			src_col[j] = row_store.get_row(j);
		// This copies from column by column.
		for (size_t i = 0; i < store.get_num_cols(); i++) {
			char *dest_col = this->get_col(i);
			sg.gather(src_col, dest_col);
			for (size_t j = 0; j < store.get_num_rows(); j++)
				src_col[j] += entry_size;
		}
	}
	else if (is_wide()) {
		const local_row_matrix_store &row_store
			= static_cast<const local_row_matrix_store &>(store);
		std::vector<const char *> src_rows(get_num_rows());
		for (size_t i = 0; i < get_num_rows(); i++)
			src_rows[i] = row_store.get_row(i);
		get_type().get_conv().conv(src_rows, get_num_cols(), get_raw_arr());
	}
	else if (store.get_raw_arr() == NULL) {
		const local_row_matrix_store &row_store
			= static_cast<const local_row_matrix_store &>(store);
		std::vector<char *> dst_row(get_num_cols());
		size_t entry_size = store.get_type().get_size();
		const scatter_gather &sg = store.get_type().get_sg();
		for (size_t j = 0; j < store.get_num_cols(); j++)
			dst_row[j] = get_col(j);

		// This copies from row by row.
		for (size_t i = 0; i < store.get_num_rows(); i++) {
			const char *src_row = row_store.get_row(i);
			sg.scatter(src_row, dst_row);
			for (size_t j = 0; j < dst_row.size(); j++)
				dst_row[j] += entry_size;
		}
	}
	else {
		std::vector<char *> dst_cols(get_num_cols());
		for (size_t i = 0; i < get_num_cols(); i++)
			dst_cols[i] = get_col(i);
		get_type().get_conv().conv(store.get_raw_arr(),
				get_num_rows() * get_num_cols(), dst_cols);
	}
	return true;
}

/*
 * The two functions below work for matrices with specific data layout.
 * It works for different matrix shapes, but it may have better performance
 * for matrices with certain shapes. We require the method in dense_matrix
 * to determine the right data layout for a matrix with a specific data layout.
 */

/*
 * In this case, the left matrix is row-major. The right matrix is stored
 * column-wise.
 */
static void inner_prod_row(const local_row_matrix_store &m1,
		const local_col_matrix_store &m2, const bulk_operate &left_op,
		const bulk_operate &right_op, local_matrix_store &res,
		std::vector<char> &tmp_buf)
{
	size_t ncol = m1.get_num_cols();
	size_t nrow = m1.get_num_rows();
	tmp_buf.resize(ncol * left_op.output_entry_size());
	for (size_t i = 0; i < nrow; i++) {
		for (size_t j = 0; j < m2.get_num_cols(); j++) {
			left_op.runAA(ncol, m1.get_row(i), m2.get_col(j), tmp_buf.data());
			right_op.runAgg(ncol, tmp_buf.data(), NULL, res.get(i, j));
		}
	}
}

/*
 * In this case, the left matrix is in column major and we don't assume
 * the layout of the right matrix.
 */
static void inner_prod_col(const local_col_matrix_store &m1,
		const local_matrix_store &m2, const bulk_operate &left_op,
		const bulk_operate &right_op, local_col_matrix_store &res,
		std::vector<char> &tmp_buf)
{
	size_t ncol = m1.get_num_cols();
	size_t nrow = m1.get_num_rows();
	tmp_buf.resize(nrow * left_op.output_entry_size());

	// For the first col from the left matrix and the first row from
	// the right matrix.
	const char *m1_col = m1.get_col(0);
	for (size_t j = 0; j < m2.get_num_cols(); j++)
		left_op.runAE(nrow, m1_col, m2.get(0, j), res.get_col(j));

	// For the rest of the matrix.
	for (size_t i = 1; i < ncol; i++) {
		const char *m1_col = m1.get_col(i);
		for (size_t j = 0; j < m2.get_num_cols(); j++) {
			left_op.runAE(nrow, m1_col, m2.get(i, j), tmp_buf.data());
			char *store_col = res.get_col(j);
			right_op.runAA(nrow, tmp_buf.data(), store_col, store_col);
		}
	}
}

static void _inner_prod(const local_matrix_store &m1,
		const local_matrix_store &m2, const bulk_operate &left_op,
		const bulk_operate &right_op, local_matrix_store &res,
		std::vector<char> &tmp_buf)
{
	if (m1.store_layout() == matrix_layout_t::L_ROW) {
		assert(m2.store_layout() == matrix_layout_t::L_COL);
		inner_prod_row(static_cast<const local_row_matrix_store &>(m1),
				static_cast<const local_col_matrix_store &>(m2), left_op,
				right_op, res, tmp_buf);
	}
	else {
		assert(res.store_layout() == matrix_layout_t::L_COL);
		inner_prod_col(static_cast<const local_col_matrix_store &>(m1),
				m2, left_op, right_op,
				static_cast<local_col_matrix_store &>(res), tmp_buf);
	}
}

void inner_prod(const local_matrix_store &left, const local_matrix_store &right,
		const bulk_operate &left_op, const bulk_operate &right_op,
		local_matrix_store &res)
{
	std::vector<char> tmp_buf;
	// If the matrix is small.
	if ((left.is_wide() && left.get_num_cols() <= LONG_DIM_LEN)
			|| (!left.is_wide() && left.get_num_rows() <= LONG_DIM_LEN))
		_inner_prod(left, right, left_op, right_op, res, tmp_buf);
	// resize the wide matrix.
	else if (left.is_wide()) {
		size_t orig_num_cols = left.get_num_cols();
		local_matrix_store::exposed_area orig_left = left.get_exposed_area();
		local_matrix_store::exposed_area orig_right = right.get_exposed_area();
		local_matrix_store &mutable_left = const_cast<local_matrix_store &>(
				left);
		local_matrix_store &mutable_right = const_cast<local_matrix_store &>(
				right);
		local_matrix_store::ptr tmp_res;
		if (res.store_layout() == matrix_layout_t::L_ROW)
			tmp_res = local_matrix_store::ptr(new local_buf_row_matrix_store(
						0, 0, res.get_num_rows(), res.get_num_cols(),
						res.get_type(), -1));
		else
			tmp_res = local_matrix_store::ptr(new local_buf_col_matrix_store(
						0, 0, res.get_num_rows(), res.get_num_cols(),
						res.get_type(), -1));
		for (size_t col_idx = 0; col_idx < orig_num_cols;
				col_idx += LONG_DIM_LEN) {
			size_t llen = std::min(orig_num_cols - col_idx, LONG_DIM_LEN);
			mutable_left.resize(orig_left.local_start_row,
					orig_left.local_start_col + col_idx, left.get_num_rows(),
					llen);
			mutable_right.resize(orig_right.local_start_row + col_idx,
					orig_right.local_start_col, llen, right.get_num_cols());
			// We accumulate the product on the result matrix directly.
			_inner_prod(left, right, left_op, right_op, *tmp_res, tmp_buf);
			mapply2(res, *tmp_res, right_op, res);
		}
		mutable_left.restore_size(orig_left);
		mutable_right.restore_size(orig_right);
	}
	// resize the tall matrix
	else {
		size_t orig_num_rows = left.get_num_rows();
		local_matrix_store::exposed_area orig_left = left.get_exposed_area();
		local_matrix_store::exposed_area orig_res = res.get_exposed_area();
		local_matrix_store &mutable_left = const_cast<local_matrix_store &>(
				left);
		for (size_t row_idx = 0; row_idx < orig_num_rows;
				row_idx += LONG_DIM_LEN) {
			size_t llen = std::min(orig_num_rows - row_idx, LONG_DIM_LEN);
			mutable_left.resize(orig_left.local_start_row + row_idx,
					orig_left.local_start_col, llen, left.get_num_cols());
			res.resize(orig_res.local_start_row + row_idx,
					orig_res.local_start_col, llen, res.get_num_cols());
			_inner_prod(left, right, left_op, right_op, res, tmp_buf);
		}
		mutable_left.restore_size(orig_left);
		res.restore_size(orig_res);
	}
}

// This case is to aggregate all elements to a single value.
static void agg_both(const local_matrix_store &store, const agg_operate &op,
		local_vec_store &res)
{
	size_t ncol = store.get_num_cols();
	size_t nrow = store.get_num_rows();
	assert(res.get_length() == 1);
	// If the store has data stored contiguously.
	if (store.get_raw_arr())
		op.runAgg(ncol * nrow, store.get_raw_arr(), NULL, res.get_raw_arr());
	// For row-major matrix and the agg op allows to combine partial agg res.
	else if (store.store_layout() == matrix_layout_t::L_ROW && op.has_combine()) {
		const local_row_matrix_store &row_store
			= static_cast<const local_row_matrix_store &>(store);
		local_buf_vec_store part_res(0, store.get_num_rows(),
				op.get_output_type(), -1);
		for (size_t i = 0; i < nrow; i++)
			op.runAgg(ncol, row_store.get_row(i), NULL, part_res.get(i));
		op.runCombine(part_res.get_length(), part_res.get_raw_arr(), NULL,
				res.get_raw_arr());
	}
	// For col-major matrix and the agg op allows to combine partial agg res.
	else if (op.has_combine()) {
		assert(store.store_layout() == matrix_layout_t::L_COL);
		const local_col_matrix_store &col_store
			= static_cast<const local_col_matrix_store &>(store);
		local_buf_vec_store part_res(0, store.get_num_cols(),
				op.get_output_type(), -1);
		for (size_t i = 0; i < ncol; i++)
			op.runAgg(nrow, col_store.get_col(i), NULL, part_res.get(i));
		op.runCombine(part_res.get_length(), part_res.get_raw_arr(), NULL,
				res.get_raw_arr());
	}
	// For row-major matrix and the agg op doesn't allow to combine partial
	// agg res.
	else if (store.store_layout() == matrix_layout_t::L_ROW) {
		local_buf_row_matrix_store buf(0, 0, store.get_num_rows(),
				store.get_num_cols(), store.get_type(), -1);
		buf.copy_from(store);
		op.runAgg(ncol * nrow, buf.get_raw_arr(), NULL, res.get_raw_arr());
	}
	// For col-major matrix and the agg op doesn't allow to combine partial
	// agg res.
	else {
		local_buf_col_matrix_store buf(0, 0, store.get_num_rows(),
				store.get_num_cols(), store.get_type(), -1);
		buf.copy_from(store);
		op.runAgg(ncol * nrow, buf.get_raw_arr(), NULL, res.get_raw_arr());
	}
}

static void agg_rows(const local_matrix_store &store, const agg_operate &op,
		local_vec_store &res)
{
	size_t ncol = store.get_num_cols();
	size_t nrow = store.get_num_rows();
	assert(ncol > 1);
	// Aggregate on rows, but the matrix is stored in col-major.
	// Instead of running aggregation on each row directly, we compute partial
	// aggregation on columns. This only works if the aggregation operation
	// allows partial aggregation and the agg and combine operations are
	// the same.
	if (store.store_layout() == matrix_layout_t::L_COL && op.is_same()) {
		assert(res.get_length() == store.get_num_rows());
		const local_col_matrix_store &col_store
			= static_cast<const local_col_matrix_store &>(store);
		if (nrow <= LONG_DIM_LEN) {
			char *res_arr = res.get_raw_arr();
			op.get_agg().runAA(res.get_length(), col_store.get_col(0),
					col_store.get_col(1), res_arr);
			for (size_t i = 2; i < ncol; i++)
				op.get_agg().runAA(res.get_length(), res_arr,
						col_store.get_col(i), res_arr);
		}
		// In this case, we need to break the matrix into smaller partitions
		// so that the agg result is always in L1/L2 cache.
		else {
			char *res_arr = res.get_raw_arr();
			size_t res_entry_size = res.get_type().get_size();
			size_t orig_num_rows = col_store.get_num_rows();
			local_matrix_store::exposed_area orig = col_store.get_exposed_area();
			local_col_matrix_store &mutable_store
				= const_cast<local_col_matrix_store &>(col_store);
			for (size_t row_idx = 0; row_idx < orig_num_rows;
					row_idx += LONG_DIM_LEN) {
				size_t llen = std::min(orig_num_rows - row_idx, LONG_DIM_LEN);
				mutable_store.resize(row_idx, 0, llen, col_store.get_num_cols());

				char *lres_arr = res_arr + row_idx * res_entry_size;
				op.get_agg().runAA(llen, col_store.get_col(0),
						col_store.get_col(1), lres_arr);
				for (size_t i = 2; i < col_store.get_num_cols(); i++)
					op.get_agg().runAA(llen, lres_arr, col_store.get_col(i),
							lres_arr);
			}
			mutable_store.restore_size(orig);
		}
	}
	// This is the default solution to compute aggregation.
	else {
		local_matrix_store::const_ptr buf_mat;
		assert(res.get_length() == store.get_num_rows());
		const local_row_matrix_store *row_store;
		if (store.store_layout() == matrix_layout_t::L_COL) {
			buf_mat = store.conv2(matrix_layout_t::L_ROW);
			assert(buf_mat);
			row_store = static_cast<const local_row_matrix_store *>(
					buf_mat.get());
		}
		else
			row_store = static_cast<const local_row_matrix_store *>(&store);
		// Aggregate on rows, and the matrix is stored in row-major.
		for (size_t i = 0; i < nrow; i++)
			op.runAgg(ncol, row_store->get_row(i), NULL, res.get(i));
	}
}

// This is implemented very similarly to above.
static void agg_cols(const local_matrix_store &store, const agg_operate &op,
		local_vec_store &res)
{
	size_t ncol = store.get_num_cols();
	size_t nrow = store.get_num_rows();
	assert(nrow > 1);
	if (store.store_layout() == matrix_layout_t::L_ROW && op.is_same()) {
		assert(res.get_length() == store.get_num_cols());
		const local_row_matrix_store &row_store
			= static_cast<const local_row_matrix_store &>(store);
		if (ncol <= LONG_DIM_LEN) {
			char *res_arr = res.get_raw_arr();
			op.get_agg().runAA(res.get_length(), row_store.get_row(0),
					row_store.get_row(1), res_arr);
			for (size_t i = 2; i < nrow; i++)
				op.get_agg().runAA(res.get_length(), res_arr,
						row_store.get_row(i), res_arr);
		}
		else {
			char *res_arr = res.get_raw_arr();
			size_t res_entry_size = res.get_type().get_size();
			size_t orig_num_cols = row_store.get_num_cols();
			local_matrix_store::exposed_area orig = row_store.get_exposed_area();
			local_row_matrix_store &mutable_store
				= const_cast<local_row_matrix_store &>(row_store);
			for (size_t col_idx = 0; col_idx < orig_num_cols;
					col_idx += LONG_DIM_LEN) {
				size_t llen = std::min(orig_num_cols - col_idx, LONG_DIM_LEN);
				mutable_store.resize(0, col_idx, row_store.get_num_rows(), llen);

				char *lres_arr = res_arr + col_idx * res_entry_size;
				op.get_agg().runAA(llen, row_store.get_row(0),
						row_store.get_row(1), lres_arr);
				for (size_t i = 2; i < row_store.get_num_rows(); i++)
					op.get_agg().runAA(llen, lres_arr, row_store.get_row(i),
							lres_arr);
			}
			mutable_store.restore_size(orig);
		}
	}
	else {
		local_matrix_store::const_ptr buf_mat;
		assert(res.get_length() == store.get_num_cols());
		const local_col_matrix_store *col_store;
		if (store.store_layout() == matrix_layout_t::L_ROW) {
			buf_mat = store.conv2(matrix_layout_t::L_COL);
			assert(buf_mat);
			col_store = static_cast<const local_col_matrix_store *>(
					buf_mat.get());
		}
		else
			col_store = static_cast<const local_col_matrix_store *>(&store);
		for (size_t i = 0; i < ncol; i++)
			op.runAgg(nrow, col_store->get_col(i), NULL, res.get(i));
	}
}

void aggregate(const local_matrix_store &store, const agg_operate &op,
		int margin, local_vec_store &res)
{
	if (margin == matrix_margin::BOTH)
		agg_both(store, op, res);
	else if (margin == matrix_margin::MAR_ROW)
		agg_rows(store, op, res);
	else if (margin == matrix_margin::MAR_COL)
		agg_cols(store, op, res);
	else {
		// This shouldn't happen.
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"aggregate on an unknown margin %1%") % margin;
		assert(0);
	}
}

static void _mapply2(const local_matrix_store &m1, const local_matrix_store &m2,
			const bulk_operate &op, local_matrix_store &res)
{
	assert(m1.store_layout() == m2.store_layout()
			&& m1.store_layout() == res.store_layout());
	size_t ncol = m1.get_num_cols();
	size_t nrow = m1.get_num_rows();
	// If the store has data stored contiguously.
	if (m1.get_raw_arr() && m2.get_raw_arr() && res.get_raw_arr())
		op.runAA(ncol * nrow, m1.get_raw_arr(), m2.get_raw_arr(),
				res.get_raw_arr());
	else if (m1.store_layout() == matrix_layout_t::L_ROW) {
		const local_row_matrix_store &row_m1 = (const local_row_matrix_store &) m1;
		const local_row_matrix_store &row_m2 = (const local_row_matrix_store &) m2;
		local_row_matrix_store &row_res = (local_row_matrix_store &) res;
		for (size_t i = 0; i < nrow; i++)
			op.runAA(ncol, row_m1.get_row(i), row_m2.get_row(i),
					row_res.get_row(i));
	}
	else {
		assert(m1.store_layout() == matrix_layout_t::L_COL);
		const local_col_matrix_store &col_m1 = (const local_col_matrix_store &) m1;
		const local_col_matrix_store &col_m2 = (const local_col_matrix_store &) m2;
		local_col_matrix_store &col_res = (local_col_matrix_store &) res;
		for (size_t i = 0; i < ncol; i++)
			op.runAA(nrow, col_m1.get_col(i), col_m2.get_col(i),
					col_res.get_col(i));
	}
}

void mapply2(const local_matrix_store &m1, const local_matrix_store &m2,
			const bulk_operate &op, local_matrix_store &res)
{
	bool is_virt = m1.is_virtual() || m2.is_virtual();
	// resize the wide matrix.
	if (is_virt && m1.is_wide() && m1.get_num_cols() > LONG_DIM_LEN) {
		size_t orig_num_cols = m1.get_num_cols();
		local_matrix_store::exposed_area orig_m1 = m1.get_exposed_area();
		local_matrix_store::exposed_area orig_m2 = m2.get_exposed_area();
		local_matrix_store::exposed_area orig_res = res.get_exposed_area();
		local_matrix_store &mutable_m1 = const_cast<local_matrix_store &>(m1);
		local_matrix_store &mutable_m2 = const_cast<local_matrix_store &>(m2);
		for (size_t col_idx = 0; col_idx < orig_num_cols;
				col_idx += LONG_DIM_LEN) {
			size_t llen = std::min(orig_num_cols - col_idx, LONG_DIM_LEN);
			mutable_m1.resize(orig_m1.local_start_row,
					orig_m1.local_start_col + col_idx, m1.get_num_rows(), llen);
			mutable_m2.resize(orig_m2.local_start_row,
					orig_m2.local_start_col + col_idx, m2.get_num_rows(), llen);
			res.resize(orig_res.local_start_row,
					orig_res.local_start_col + col_idx, res.get_num_rows(), llen);
			_mapply2(m1, m2, op, res);
		}
		mutable_m1.restore_size(orig_m1);
		mutable_m2.restore_size(orig_m2);
		res.restore_size(orig_res);
	}
	// resize the tall matrix
	else if (is_virt && m1.get_num_rows() > LONG_DIM_LEN) {
		size_t orig_num_rows = m1.get_num_rows();
		local_matrix_store::exposed_area orig_m1 = m1.get_exposed_area();
		local_matrix_store::exposed_area orig_m2 = m2.get_exposed_area();
		local_matrix_store::exposed_area orig_res = res.get_exposed_area();
		local_matrix_store &mutable_m1 = const_cast<local_matrix_store &>(m1);
		local_matrix_store &mutable_m2 = const_cast<local_matrix_store &>(m2);
		for (size_t row_idx = 0; row_idx < orig_num_rows;
				row_idx += LONG_DIM_LEN) {
			size_t llen = std::min(orig_num_rows - row_idx, LONG_DIM_LEN);
			mutable_m1.resize(orig_m1.local_start_row + row_idx,
					orig_m1.local_start_col, llen, m1.get_num_cols());
			mutable_m2.resize(orig_m2.local_start_row + row_idx,
					orig_m2.local_start_col, llen, m2.get_num_cols());
			res.resize(orig_res.local_start_row + row_idx,
					orig_res.local_start_col, llen, res.get_num_cols());
			_mapply2(m1, m2, op, res);
		}
		mutable_m1.restore_size(orig_m1);
		mutable_m2.restore_size(orig_m2);
		res.restore_size(orig_res);
	}
	else
		// If the local matrix isn't virtual, we don't need to resize it
		// to increase CPU cache hits.
		_mapply2(m1, m2, op, res);
}

static void _sapply(const local_matrix_store &store, const bulk_uoperate &op,
		local_matrix_store &res)
{
	assert(res.store_layout() == store.store_layout());
	size_t ncol = store.get_num_cols();
	size_t nrow = store.get_num_rows();
	// If the store has data stored contiguously.
	if (store.get_raw_arr() && res.get_raw_arr())
		op.runA(ncol * nrow, store.get_raw_arr(), res.get_raw_arr());
	else if (store.store_layout() == matrix_layout_t::L_ROW) {
		const local_row_matrix_store &row_store
			= (const local_row_matrix_store &) store;
		local_row_matrix_store &row_res = (local_row_matrix_store &) res;
		for (size_t i = 0; i < nrow; i++)
			op.runA(ncol, row_store.get_row(i), row_res.get_row(i));
	}
	else {
		assert(store.store_layout() == matrix_layout_t::L_COL);
		const local_col_matrix_store &col_store
			= (const local_col_matrix_store &) store;
		local_col_matrix_store &col_res = (local_col_matrix_store &) res;
		for (size_t i = 0; i < ncol; i++)
			op.runA(nrow, col_store.get_col(i), col_res.get_col(i));
	}
}

void sapply(const local_matrix_store &store, const bulk_uoperate &op,
		local_matrix_store &res)
{
	// resize the wide matrix.
	if (store.is_virtual() && store.is_wide()
			&& store.get_num_cols() > LONG_DIM_LEN) {
		size_t orig_num_cols = store.get_num_cols();
		local_matrix_store::exposed_area orig_in = store.get_exposed_area();
		local_matrix_store::exposed_area orig_res = res.get_exposed_area();
		local_matrix_store &mutable_store = const_cast<local_matrix_store &>(
				store);
		for (size_t col_idx = 0; col_idx < orig_num_cols;
				col_idx += LONG_DIM_LEN) {
			size_t llen = std::min(orig_num_cols - col_idx, LONG_DIM_LEN);
			mutable_store.resize(orig_in.local_start_row,
					orig_in.local_start_col + col_idx, store.get_num_rows(), llen);
			res.resize(orig_res.local_start_row,
					orig_res.local_start_col + col_idx, res.get_num_rows(), llen);
			_sapply(store, op, res);
		}
		mutable_store.restore_size(orig_in);
		res.restore_size(orig_res);
	}
	// resize the tall matrix
	else if (store.is_virtual() && store.get_num_rows() > LONG_DIM_LEN) {
		size_t orig_num_rows = store.get_num_rows();
		local_matrix_store::exposed_area orig_in = store.get_exposed_area();
		local_matrix_store::exposed_area orig_res = res.get_exposed_area();
		local_matrix_store &mutable_store = const_cast<local_matrix_store &>(
				store);
		for (size_t row_idx = 0; row_idx < orig_num_rows;
				row_idx += LONG_DIM_LEN) {
			size_t llen = std::min(orig_num_rows - row_idx, LONG_DIM_LEN);
			mutable_store.resize(orig_in.local_start_row + row_idx,
					orig_in.local_start_col, llen, store.get_num_cols());
			res.resize(orig_res.local_start_row + row_idx,
					orig_res.local_start_col, llen, res.get_num_cols());
			_sapply(store, op, res);
		}
		mutable_store.restore_size(orig_in);
		res.restore_size(orig_res);
	}
	else
		// If the local matrix isn't virtual, we don't need to resize it
		// to increase CPU cache hits.
		_sapply(store, op, res);
}

void apply(int margin, const arr_apply_operate &op,
		const local_matrix_store &in_mat, local_matrix_store &out_mat)
{
	assert(margin == matrix_margin::MAR_ROW || margin == matrix_margin::MAR_COL);
	// In these two cases, we need to convert the matrix store layout
	// before we can apply the function to the matrix.
	local_matrix_store::const_ptr buf_mat;
	if (in_mat.store_layout() == matrix_layout_t::L_COL
			&& margin == matrix_margin::MAR_ROW) {
		buf_mat = in_mat.conv2(matrix_layout_t::L_ROW);
		assert(buf_mat);
	}
	else if (in_mat.store_layout() == matrix_layout_t::L_ROW
			&& margin == matrix_margin::MAR_COL) {
		buf_mat = in_mat.conv2(matrix_layout_t::L_COL);
		assert(buf_mat);
	}

	const local_matrix_store *this_mat;
	if (buf_mat)
		this_mat = buf_mat.get();
	else
		this_mat = &in_mat;

	if (margin == matrix_margin::MAR_ROW) {
		assert(this_mat->store_layout() == matrix_layout_t::L_ROW);
		assert(out_mat.store_layout() == matrix_layout_t::L_ROW);
		const local_row_matrix_store &row_mat
			= static_cast<const local_row_matrix_store &>(*this_mat);
		local_row_matrix_store &ret_row_mat
			= static_cast<local_row_matrix_store &>(out_mat);
		for (size_t i = 0; i < row_mat.get_num_rows(); i++) {
			local_cref_vec_store in_vec(row_mat.get_row(i), 0,
					row_mat.get_num_cols(), row_mat.get_type(), -1);
			local_ref_vec_store out_vec(ret_row_mat.get_row(i), 0,
					ret_row_mat.get_num_cols(), ret_row_mat.get_type(), -1);
			op.run(in_vec, out_vec);
		}
	}
	else {
		assert(this_mat->store_layout() == matrix_layout_t::L_COL);
		assert(out_mat.store_layout() == matrix_layout_t::L_COL);
		const local_col_matrix_store &col_mat
			= static_cast<const local_col_matrix_store &>(*this_mat);
		local_col_matrix_store &ret_col_mat
			= static_cast<local_col_matrix_store &>(out_mat);
		for (size_t i = 0; i < col_mat.get_num_cols(); i++) {
			local_cref_vec_store in_vec(col_mat.get_col(i), 0,
					col_mat.get_num_rows(), col_mat.get_type(), -1);
			local_ref_vec_store out_vec(ret_col_mat.get_col(i), 0,
					ret_col_mat.get_num_rows(), ret_col_mat.get_type(), -1);
			op.run(in_vec, out_vec);
		}
	}
}

namespace
{

/*
 * There are rounding error problems when multiplying float-points.
 * We deal with float-point multiplication specially to reduce rounding error.
 */
class double_multiply_operate: public bulk_operate
{
public:
	virtual void runAA(size_t num_eles, const void *left_arr,
			const void *right_arr, void *output_arr) const {
		const double *a = static_cast<const double *>(left_arr);
		const double *b = static_cast<const double *>(right_arr);
		double *c = static_cast<double *>(output_arr);
		for (size_t i = 0; i < num_eles; i++)
			c[i] = ((long double) a[i]) * b[i];
	}
	virtual void runAE(size_t num_eles, const void *left_arr,
			const void *right, void *output_arr) const {
		double a = *static_cast<const double *>(right);
		const double *x = static_cast<const double *>(left_arr);
		double *c = static_cast<double *>(output_arr);
		cblas_dcopy(num_eles, x, 1, c, 1);
		cblas_dscal(num_eles, a, c, 1);
	}
	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr, void *output_arr) const {
		double a = *static_cast<const double *>(left);
		const double *x = static_cast<const double *>(right_arr);
		double *c = static_cast<double *>(output_arr);
		cblas_dcopy(num_eles, x, 1, c, 1);
		cblas_dscal(num_eles, a, c, 1);
	}
	virtual void runAgg(size_t num_eles, const void *left_arr,
			const void *orig, void *output) const {
		assert(0);
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<double>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<double>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<double>();
	}
	virtual std::string get_name() const {
		return "dmultiply";
	}
};

double_multiply_operate dm_op;

}

void mapply_rows(const local_matrix_store &store, const local_vec_store &vals,
		const bulk_operate &_op, local_matrix_store &res)
{
	assert(res.store_layout() == store.store_layout());
	assert(store.get_num_rows() == res.get_num_rows());
	assert(store.get_num_cols() == res.get_num_cols());
	assert(store.get_num_cols() == vals.get_length());
	size_t ncol = store.get_num_cols();
	size_t nrow = store.get_num_rows();
	const bulk_operate *op = &_op;
	if (op == &get_scalar_type<double>().get_basic_ops().get_multiply())
		op = &dm_op;
	if (store.store_layout() == matrix_layout_t::L_ROW) {
		const local_row_matrix_store &row_store
			= static_cast<const local_row_matrix_store &>(store);
		local_row_matrix_store &row_res
			= static_cast<local_row_matrix_store &>(res);
		for (size_t i = 0; i < nrow; i++)
			op->runAA(ncol, row_store.get_row(i), vals.get_raw_arr(),
					row_res.get_row(i));
	}
	else {
		assert(store.store_layout() == matrix_layout_t::L_COL);
		const local_col_matrix_store &col_store
			= static_cast<const local_col_matrix_store &>(store);
		local_col_matrix_store &col_res
			= static_cast<local_col_matrix_store &>(res);
		for (size_t i = 0; i < ncol; i++)
			op->runAE(nrow, col_store.get_col(i), vals.get(i),
					col_res.get_col(i));
	}
}

void mapply_cols(const local_matrix_store &store, const local_vec_store &vals,
		const bulk_operate &_op, local_matrix_store &res)
{
	assert(res.store_layout() == store.store_layout());
	assert(store.get_num_rows() == res.get_num_rows());
	assert(store.get_num_cols() == res.get_num_cols());
	assert(store.get_num_rows() == vals.get_length());
	size_t ncol = store.get_num_cols();
	size_t nrow = store.get_num_rows();
	const bulk_operate *op = &_op;
	if (op == &get_scalar_type<double>().get_basic_ops().get_multiply())
		op = &dm_op;
	if (store.store_layout() == matrix_layout_t::L_ROW) {
		const local_row_matrix_store &row_store
			= static_cast<const local_row_matrix_store &>(store);
		local_row_matrix_store &row_res
			= static_cast<local_row_matrix_store &>(res);
		for (size_t i = 0; i < nrow; i++)
			op->runAE(ncol, row_store.get_row(i), vals.get(i),
					row_res.get_row(i));
	}
	else {
		assert(store.store_layout() == matrix_layout_t::L_COL);
		const local_col_matrix_store &col_store
			= static_cast<const local_col_matrix_store &>(store);
		local_col_matrix_store &col_res
			= static_cast<local_col_matrix_store &>(res);
		for (size_t i = 0; i < ncol; i++)
			op->runAA(nrow, col_store.get_col(i), vals.get_raw_arr(),
					col_res.get_col(i));
	}
}

static bool _groupby_row(const detail::local_matrix_store &labels,
		const detail::local_row_matrix_store &mat, const agg_operate &op,
		detail::local_row_matrix_store &results, std::vector<bool> &agg_flags)
{
	for (size_t i = 0; i < mat.get_num_rows(); i++) {
		factor_value_t label_id = labels.get<factor_value_t>(i, 0);
		if ((size_t) label_id >= agg_flags.size())
			return false;
		// If we never get partially aggregated result for a label, we should
		// copy the data to the corresponding row.
		if (!agg_flags[label_id])
			memcpy(results.get_row(label_id), mat.get_row(i),
					mat.get_num_cols() * mat.get_entry_size());
		else
			op.get_agg().runAA(mat.get_num_cols(), mat.get_row(i),
					results.get_row(label_id), results.get_row(label_id));
		auto bit = agg_flags[label_id];
		bit = true;
	}
	return true;
}

bool groupby_row(const detail::local_matrix_store &labels,
		const detail::local_row_matrix_store &mat, const agg_operate &op,
		detail::local_row_matrix_store &results, std::vector<bool> &agg_flags)
{
	assert(!mat.is_wide());
	// resize the tall matrix
	if (mat.is_virtual() && mat.get_num_rows() > LONG_DIM_LEN) {
		size_t orig_num_rows = mat.get_num_rows();
		local_matrix_store::exposed_area orig_in = mat.get_exposed_area();
		local_matrix_store::exposed_area orig_labels = labels.get_exposed_area();
		local_row_matrix_store &mutable_mat
			= const_cast<local_row_matrix_store &>(mat);
		local_matrix_store &mutable_labels = const_cast<local_matrix_store &>(
				labels);
		for (size_t row_idx = 0; row_idx < orig_num_rows;
				row_idx += LONG_DIM_LEN) {
			size_t llen = std::min(orig_num_rows - row_idx, LONG_DIM_LEN);
			mutable_mat.resize(orig_in.local_start_row + row_idx,
					orig_in.local_start_col, llen, mat.get_num_cols());
			mutable_labels.resize(orig_labels.local_start_row + row_idx,
					orig_labels.local_start_col, llen, labels.get_num_cols());
			bool ret = _groupby_row(labels, mat, op, results, agg_flags);
			if (!ret) {
				mutable_mat.restore_size(orig_in);
				mutable_labels.restore_size(orig_labels);
				return false;
			}
		}
		mutable_mat.restore_size(orig_in);
		mutable_labels.restore_size(orig_labels);
		return true;
	}
	else
		return _groupby_row(labels, mat, op, results, agg_flags);
}

local_matrix_store::const_ptr local_row_matrix_store::get_portion(
		size_t local_start_row, size_t local_start_col, size_t num_rows,
		size_t num_cols) const
{
	if (local_start_row + num_rows > get_num_rows()
			|| local_start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) <<
			"get row portion from local matrix: out of boundary";
		return local_matrix_store::const_ptr();
	}

	local_row_matrix_store::const_ptr ret;
	size_t global_start_row = get_global_start_row() + local_start_row;
	size_t global_start_col = get_global_start_col() + local_start_col;
	if (is_wide()) {
		std::vector<const char *> rows(num_rows);
		for (size_t i = 0; i < num_rows; i++)
			rows[i] = get_row(local_start_row + i)
				+ local_start_col * get_entry_size();
		ret = local_row_matrix_store::const_ptr(new local_cref_row_matrix_store(
					orig_data_ref, rows, global_start_row, global_start_col,
					num_rows, num_cols, get_type(), get_node_id()));
	}
	else {
		assert(local_start_col == 0 && num_cols == get_num_cols());
		off_t start_off = get_num_cols() * local_start_row * get_entry_size();
		ret = local_row_matrix_store::const_ptr(new local_cref_contig_row_matrix_store(
					orig_data_ref, get_raw_arr() + start_off,
					global_start_row, global_start_col,
					num_rows, num_cols, get_type(), get_node_id()));
	}
	// If the current local matrix has the original data, the returned local matrix
	// should also hold the data.
	if (this->hold_orig_data())
		assert(ret->hold_orig_data());
	return ret;
}

local_matrix_store::const_ptr local_col_matrix_store::get_portion(
		size_t local_start_row, size_t local_start_col, size_t num_rows,
		size_t num_cols) const
{
	if (local_start_row + num_rows > get_num_rows()
			|| local_start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error)
			<< "get col portion from local matrix: out of boundary";
		return local_matrix_store::const_ptr();
	}

	local_col_matrix_store::const_ptr ret;
	size_t global_start_row = get_global_start_row() + local_start_row;
	size_t global_start_col = get_global_start_col() + local_start_col;
	if (is_wide()) {
		assert(local_start_row == 0 && num_rows == get_num_rows());
		off_t start_off = get_num_rows() * local_start_col * get_entry_size();
		ret = local_col_matrix_store::const_ptr(new local_cref_contig_col_matrix_store(
					orig_data_ref, get_raw_arr() + start_off,
					global_start_row, global_start_col,
					num_rows, num_cols, get_type(), get_node_id()));
	}
	else {
		std::vector<const char *> cols(num_cols);
		for (size_t i = 0; i < num_cols; i++)
			cols[i] = get_col(local_start_col + i)
				+ local_start_row * get_entry_size();
		ret = local_col_matrix_store::const_ptr(new local_cref_col_matrix_store(
					orig_data_ref, cols, global_start_row, global_start_col,
					num_rows, num_cols, get_type(), get_node_id()));
	}
	// If the current local matrix has the original data, the returned local matrix
	// should also hold the data.
	if (this->hold_orig_data())
		assert(ret->hold_orig_data());
	return ret;
}

local_matrix_store::ptr local_row_matrix_store::get_portion(
		size_t local_start_row, size_t local_start_col, size_t num_rows,
		size_t num_cols)
{
	if (read_only())
		return local_matrix_store::ptr();

	if (local_start_row + num_rows > get_num_rows()
			|| local_start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error)
			<< "get const row portion from local matrix: out of boundary";
		return local_matrix_store::ptr();
	}

	local_row_matrix_store::ptr ret;
	size_t global_start_row = get_global_start_row() + local_start_row;
	size_t global_start_col = get_global_start_col() + local_start_col;
	if (is_wide()) {
		std::vector<char *> rows(num_rows);
		for (size_t i = 0; i < num_rows; i++)
			rows[i] = get_row(local_start_row + i)
				+ local_start_col * get_entry_size();
		ret = local_row_matrix_store::ptr(new local_ref_row_matrix_store(
					orig_data_ref, rows, global_start_row, global_start_col,
					num_rows, num_cols, get_type(), get_node_id()));
	}
	else {
		assert(local_start_col == 0 && num_cols == get_num_cols());
		off_t start_off = get_num_cols() * local_start_row * get_entry_size();
		ret = local_row_matrix_store::ptr(new local_ref_contig_row_matrix_store(
					orig_data_ref, get_raw_arr() + start_off,
					global_start_row, global_start_col,
					num_rows, num_cols, get_type(), get_node_id()));
	}
	// If the current local matrix has the original data, the returned local matrix
	// should also hold the data.
	if (this->hold_orig_data())
		assert(ret->hold_orig_data());
	return ret;
}

local_matrix_store::ptr local_col_matrix_store::get_portion(
		size_t local_start_row, size_t local_start_col, size_t num_rows,
		size_t num_cols)
{
	if (read_only())
		return local_matrix_store::ptr();

	if (local_start_row + num_rows > get_num_rows()
			|| local_start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error)
			<< "get const col portion from local matrix: out of boundary";
		return local_matrix_store::ptr();
	}

	local_col_matrix_store::ptr ret;
	size_t global_start_row = get_global_start_row() + local_start_row;
	size_t global_start_col = get_global_start_col() + local_start_col;
	if (is_wide()) {
		assert(local_start_row == 0 && num_rows == get_num_rows());
		off_t start_off = get_num_rows() * local_start_col * get_entry_size();
		ret = local_col_matrix_store::ptr(new local_ref_contig_col_matrix_store(
					orig_data_ref, get_raw_arr() + start_off,
					global_start_row, global_start_col,
					num_rows, num_cols, get_type(), get_node_id()));
	}
	else {
		std::vector<char *> cols(num_cols);
		for (size_t i = 0; i < num_cols; i++)
			cols[i] = get_col(local_start_col + i)
				+ local_start_row * get_entry_size();
		ret = local_col_matrix_store::ptr(new local_ref_col_matrix_store(
					orig_data_ref, cols, global_start_row, global_start_col,
					num_rows, num_cols, get_type(), get_node_id()));
	}
	// If the current local matrix has the original data, the returned local matrix
	// should also hold the data.
	if (this->hold_orig_data())
		assert(ret->hold_orig_data());
	return ret;
}

}

}
