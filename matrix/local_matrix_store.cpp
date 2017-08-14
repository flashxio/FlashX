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
#include "materialize.h"
#include "factor.h"

namespace fm
{

namespace detail
{

static const size_t L1_SIZE = 1024 * 32;

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
	return true;
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

size_t local_matrix_store::get_all_rows(std::vector<const char *> &rows) const
{
	if (store_layout() == matrix_layout_t::L_COL)
		return 0;
	else if (get_raw_arr()) {
		rows.resize(get_num_rows());
		const char *arr = get_raw_arr();
		size_t row_len = get_num_cols() * get_entry_size();
		for (size_t i = 0; i < get_num_rows(); i++)
			rows[i] = arr + row_len * i;
	}
	else {
		const local_row_matrix_store *row_store
			= static_cast<const local_row_matrix_store *>(this);
		rows.resize(get_num_rows());
		for (size_t i = 0; i < get_num_rows(); i++)
			rows[i] = row_store->get_row(i);
	}
	return rows.size();
}

size_t local_matrix_store::get_all_cols(std::vector<const char *> &cols) const
{
	if (store_layout() == matrix_layout_t::L_ROW)
		return 0;
	else if (get_raw_arr()) {
		cols.resize(get_num_cols());
		const char *arr = get_raw_arr();
		size_t col_len = get_num_rows() * get_entry_size();
		for (size_t i = 0; i < get_num_cols(); i++)
			cols[i] = arr + col_len * i;
	}
	else {
		const local_col_matrix_store *col_store
			= static_cast<const local_col_matrix_store *>(this);
		cols.resize(get_num_cols());
		for (size_t i = 0; i < get_num_cols(); i++)
			cols[i] = col_store->get_col(i);
	}
	return cols.size();
}

size_t local_matrix_store::get_all_rows(std::vector<char *> &rows)
{
	if (store_layout() == matrix_layout_t::L_COL)
		return 0;
	else if (get_raw_arr()) {
		rows.resize(get_num_rows());
		char *arr = get_raw_arr();
		size_t row_len = get_num_cols() * get_entry_size();
		for (size_t i = 0; i < get_num_rows(); i++)
			rows[i] = arr + row_len * i;
	}
	else {
		local_row_matrix_store *row_store
			= static_cast<local_row_matrix_store *>(this);
		rows.resize(get_num_rows());
		for (size_t i = 0; i < get_num_rows(); i++)
			rows[i] = row_store->get_row(i);
	}
	return rows.size();
}

size_t local_matrix_store::get_all_cols(std::vector<char *> &cols)
{
	if (store_layout() == matrix_layout_t::L_ROW)
		return 0;
	else if (get_raw_arr()) {
		cols.resize(get_num_cols());
		char *arr = get_raw_arr();
		size_t col_len = get_num_rows() * get_entry_size();
		for (size_t i = 0; i < get_num_cols(); i++)
			cols[i] = arr + col_len * i;
	}
	else {
		local_col_matrix_store *col_store
			= static_cast<local_col_matrix_store *>(this);
		cols.resize(get_num_cols());
		for (size_t i = 0; i < get_num_cols(); i++)
			cols[i] = col_store->get_col(i);
	}
	return cols.size();
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
		get_type().get_conv().conv2(store.get_raw_arr(),
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
		const size_t LONG_DIM_LEN = get_long_dim_len(*this);
		char *dst_arr = get_raw_arr();
		for (size_t row_idx = 0; row_idx < get_num_rows();
				row_idx += LONG_DIM_LEN) {
			size_t llen = std::min(LONG_DIM_LEN, get_num_rows() - row_idx);
			get_type().get_conv().conv1(src_cols, llen, dst_arr);
			for (size_t i = 0; i < get_num_cols(); i++)
				src_cols[i] += llen * get_entry_size();
			dst_arr += llen * src_cols.size() * get_entry_size();
		}
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
		const size_t LONG_DIM_LEN = get_long_dim_len(*this);
		char *dst_arr = get_raw_arr();
		for (size_t col_idx = 0; col_idx < get_num_cols();
				col_idx += LONG_DIM_LEN) {
			size_t llen = std::min(LONG_DIM_LEN, get_num_cols() - col_idx);
			get_type().get_conv().conv1(src_rows, llen, dst_arr);
			for (size_t i = 0; i < get_num_rows(); i++)
				src_rows[i] += llen * get_entry_size();
			dst_arr += src_rows.size() * llen * get_entry_size();
		}
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
		get_type().get_conv().conv2(store.get_raw_arr(),
				get_num_rows() * get_num_cols(), dst_cols);
	}
	return true;
}

/*
 * We optimize inner product on matrices with different shapes and different
 * data layout differently. Specifically, we optimize inner product for
 * row-major tall matrix, row-major wide matrix, col-major tall matrix and
 * col-major wide matrix. We prefer row-major wide matrix and col-major
 * tall matrix because they are easier to be optimized.
 */

/*
 * In this case, the left matrix is row-major tall matrix. The right matrix
 * is stored column-wise.
 */
static void inner_prod_tall_row(const local_row_matrix_store &m1,
		const local_col_matrix_store &m2, const bulk_operate &left_op,
		const bulk_operate &right_op, local_matrix_store &res,
		std::vector<char> &tmp_buf)
{
	size_t ncol = m1.get_num_cols();
	size_t nrow = m1.get_num_rows();
	tmp_buf.resize(ncol * left_op.output_entry_size());

	// Get all rows in m1.
	std::vector<const char *> m1_rows;
	size_t m1_num_rows = m1.get_all_rows(m1_rows);
	assert(m1_num_rows > 0);

	// Get all cols in m2.
	std::vector<const char *> m2_cols;
	size_t m2_num_cols = m2.get_all_cols(m2_cols);
	assert(m2_num_cols > 0);

	if (res.store_layout() == matrix_layout_t::L_ROW) {
		std::vector<char *> res_rows;
		size_t res_num_rows = res.get_all_rows(res_rows);
		assert(res_num_rows > 0);

		for (size_t i = 0; i < nrow; i++) {
			for (size_t j = 0; j < m2.get_num_cols(); j++) {
				left_op.runAA(ncol, m1_rows[i], m2_cols[j], tmp_buf.data());
				right_op.runAgg(ncol, tmp_buf.data(),
						res_rows[i] + j * res.get_entry_size());
			}
		}
	}
	else {
		std::vector<char *> res_cols;
		size_t res_num_cols = res.get_all_rows(res_cols);
		assert(res_num_cols > 0);

		for (size_t i = 0; i < nrow; i++) {
			for (size_t j = 0; j < m2.get_num_cols(); j++) {
				left_op.runAA(ncol, m1_rows[i], m2_cols[j], tmp_buf.data());
				right_op.runAgg(ncol, tmp_buf.data(),
						res_cols[j] + i * res.get_entry_size());
			}
		}
	}
}

/*
 * In this case, the left matrix is row-major wide matrix. The right matrix
 * is stored column-wise.
 */
static void inner_prod_wide_row(const local_row_matrix_store &m1,
		const local_col_matrix_store &m2, const bulk_operate &left_op,
		const bulk_operate &right_op, local_matrix_store &res,
		std::vector<char> &tmp_buf)
{
	size_t ncol = m1.get_num_cols();
	size_t nrow = m1.get_num_rows();
	tmp_buf.resize(ncol * left_op.output_entry_size());

	for (size_t i = 0; i < nrow; i++) {
		const char *m1_row = m1.get_row(i);
		for (size_t j = 0; j < m2.get_num_cols(); j++) {
			left_op.runAA(ncol, m1_row, m2.get_col(j), tmp_buf.data());
			right_op.runAgg(ncol, tmp_buf.data(), res.get(i, j));
		}
	}
}

static void inner_prod_col_part(const local_col_matrix_store &m1,
		size_t m1_start_col, size_t m1_num_cols, const local_matrix_store &m2,
		size_t m2_start_col, size_t m2_num_cols, const bulk_operate &left_op,
		const bulk_operate &right_op, local_col_matrix_store &res,
		std::vector<char> &tmp_buf)
{
	size_t nrow = m1.get_num_rows();
	tmp_buf.resize(nrow * left_op.output_entry_size());
	size_t m1_end_col = m1_start_col + m1_num_cols;
	size_t m2_end_col = m2_start_col + m2_num_cols;

	// For the first col from the left matrix and the first row from
	// the right matrix.
	if (m1_start_col == 0) {
		const char *m1_col = m1.get_col(0);
		for (size_t j = m2_start_col; j < m2_start_col + m2_num_cols; j++)
			left_op.runAE(nrow, m1_col, m2.get(0, j), res.get_col(j));
		m1_start_col = 1;
	}

	// For the rest of the matrix.
	for (size_t i = m1_start_col; i < m1_end_col; i++) {
		const char *m1_col = m1.get_col(i);
		for (size_t j = m2_start_col; j < m2_end_col; j++) {
			left_op.runAE(nrow, m1_col, m2.get(i, j), tmp_buf.data());
			char *store_col = res.get_col(j);
			right_op.runAA(nrow, tmp_buf.data(), store_col, store_col);
		}
	}
}

/*
 * In this case, the left matrix is tall and in column major and we don't
 * assume the layout of the right matrix.
 */
static void inner_prod_col(const local_col_matrix_store &m1,
		const local_matrix_store &m2, const bulk_operate &left_op,
		const bulk_operate &right_op, local_col_matrix_store &res,
		std::vector<char> &tmp_buf)
{
	size_t num_cols_part = L1_SIZE / m1.get_num_rows() / m1.get_entry_size() / 2;
	if (num_cols_part <= 1)
		inner_prod_col_part(m1, 0, m1.get_num_cols(), m2, 0, m2.get_num_cols(),
				left_op, right_op, res, tmp_buf);
	else {
		// We further partition the matrix vertically to keep the columns
		// in the L1 CPU cache as much as possible. However, it doesn't seem
		// to have performance improvement.
		// TODO We need more performance profiling to improve its performance.
		for (size_t i = 0; i < m1.get_num_cols(); i += num_cols_part) {
			size_t m1_num_cols = std::min(num_cols_part, m1.get_num_cols() - i);
			for (size_t j = 0; j < m2.get_num_cols(); j += num_cols_part) {
				size_t m2_num_cols = std::min(num_cols_part,
						m2.get_num_cols() - j);
				inner_prod_col_part(m1, i, m1_num_cols, m2, j, m2_num_cols,
						left_op, right_op, res, tmp_buf);
			}
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
		const local_row_matrix_store &row_m1
			= static_cast<const local_row_matrix_store &>(m1);
		const local_col_matrix_store &col_m2
			= static_cast<const local_col_matrix_store &>(m2);
		if (m1.is_wide())
			inner_prod_wide_row(row_m1, col_m2, left_op, right_op, res, tmp_buf);
		else
			inner_prod_tall_row(row_m1, col_m2, left_op, right_op, res, tmp_buf);
	}
	else {
		assert(res.store_layout() == matrix_layout_t::L_COL);
		const local_col_matrix_store &col_m1
			= static_cast<const local_col_matrix_store &>(m1);
		local_col_matrix_store &col_res
			= static_cast<local_col_matrix_store &>(res);
		inner_prod_col(col_m1, m2, left_op, right_op, col_res, tmp_buf);
	}
	m1.complete();
	m2.complete();
}

void inner_prod_wide(const local_matrix_store &left, const local_matrix_store &right,
		const bulk_operate &left_op, const bulk_operate &right_op,
		local_matrix_store &res)
{
	size_t part_len = std::min(get_part_dim_len(left, part_dim_t::PART_DIM2),
			get_part_dim_len(right, part_dim_t::PART_DIM1));
	std::vector<char> tmp_buf;
	// If the matrix is small.
	// Or one of the input matrices has been resized.
	if (left.get_num_cols() <= part_len || !left.is_whole() || !right.is_whole())
		_inner_prod(left, right, left_op, right_op, res, tmp_buf);
	// resize the wide matrix.
	else {
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
		bool success = true;
		for (size_t col_idx = 0; col_idx < orig_num_cols;
				col_idx += part_len) {
			size_t llen = std::min(orig_num_cols - col_idx, part_len);
			success = mutable_left.resize(orig_left.local_start_row,
					orig_left.local_start_col + col_idx, left.get_num_rows(),
					llen);
			if (success)
				success = mutable_right.resize(
						orig_right.local_start_row + col_idx,
						orig_right.local_start_col, llen, right.get_num_cols());
			if (!success) {
				// We assume if resize fails, it fails in the first iteration.
				assert(col_idx == 0);
				break;
			}

			if (col_idx > 0) {
				// We accumulate the product on the result matrix directly.
				_inner_prod(left, right, left_op, right_op, *tmp_res, tmp_buf);
				mapply2(res, *tmp_res, right_op, part_dim_t::PART_NONE, res);
			}
			else
				_inner_prod(left, right, left_op, right_op, res, tmp_buf);
		}
		mutable_left.restore_size(orig_left);
		mutable_right.restore_size(orig_right);
		if (!success)
			_inner_prod(left, right, left_op, right_op, res, tmp_buf);
	}
}

// In this case, the left matrix is a tall matrix and the right matrix
// is a small matrix.
void inner_prod_tall(const local_matrix_store &left, const local_matrix_store &right,
		const bulk_operate &left_op, const bulk_operate &right_op,
		local_matrix_store &res)
{
	size_t part_len = get_part_dim_len(left, part_dim_t::PART_DIM1);
	std::vector<char> tmp_buf;
	// If the matrix is small.
	// Or the left matrices has been resized.
	if (left.get_num_rows() <= part_len || !left.is_whole())
		_inner_prod(left, right, left_op, right_op, res, tmp_buf);
	else {
		size_t orig_num_rows = left.get_num_rows();
		local_matrix_store::exposed_area orig_left = left.get_exposed_area();
		local_matrix_store::exposed_area orig_res = res.get_exposed_area();
		local_matrix_store &mutable_left = const_cast<local_matrix_store &>(
				left);
		bool success = true;
		for (size_t row_idx = 0; row_idx < orig_num_rows; row_idx += part_len) {
			size_t llen = std::min(orig_num_rows - row_idx, part_len);
			success = mutable_left.resize(orig_left.local_start_row + row_idx,
					orig_left.local_start_col, llen, left.get_num_cols());
			if (success)
				success = res.resize(orig_res.local_start_row + row_idx,
						orig_res.local_start_col, llen, res.get_num_cols());
			if (!success) {
				// We assume if resize fails, it fails in the first iteration.
				assert(row_idx == 0);
				break;
			}
			_inner_prod(left, right, left_op, right_op, res, tmp_buf);
		}
		mutable_left.restore_size(orig_left);
		res.restore_size(orig_res);
		if (!success)
			_inner_prod(left, right, left_op, right_op, res, tmp_buf);
	}
}

// This case is to aggregate all elements to a single value.
static void agg_both(const local_matrix_store &store, const agg_operate &op,
		local_matrix_store &res)
{
	size_t ncol = store.get_num_cols();
	size_t nrow = store.get_num_rows();
	assert(res.get_num_rows() == 1 && res.get_num_cols() == 1);
	// If the store has data stored contiguously.
	if (store.get_raw_arr())
		op.runAgg(ncol * nrow, store.get_raw_arr(), res.get_raw_arr());
	// For row-major matrix and the agg op allows to combine partial agg res.
	else if (store.store_layout() == matrix_layout_t::L_ROW && op.has_combine()) {
		const local_row_matrix_store &row_store
			= static_cast<const local_row_matrix_store &>(store);
		local_buf_vec_store part_res(0, store.get_num_rows(),
				op.get_output_type(), -1);
		for (size_t i = 0; i < nrow; i++)
			op.runAgg(ncol, row_store.get_row(i), part_res.get(i));
		op.runCombine(part_res.get_length(), part_res.get_raw_arr(),
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
			op.runAgg(nrow, col_store.get_col(i), part_res.get(i));
		op.runCombine(part_res.get_length(), part_res.get_raw_arr(),
				res.get_raw_arr());
	}
	// For row-major matrix and the agg op doesn't allow to combine partial
	// agg res.
	else if (store.store_layout() == matrix_layout_t::L_ROW) {
		local_buf_row_matrix_store buf(0, 0, store.get_num_rows(),
				store.get_num_cols(), store.get_type(), -1);
		buf.copy_from(store);
		op.runAgg(ncol * nrow, buf.get_raw_arr(), res.get_raw_arr());
	}
	// For col-major matrix and the agg op doesn't allow to combine partial
	// agg res.
	else {
		local_buf_col_matrix_store buf(0, 0, store.get_num_rows(),
				store.get_num_cols(), store.get_type(), -1);
		buf.copy_from(store);
		op.runAgg(ncol * nrow, buf.get_raw_arr(), res.get_raw_arr());
	}
}

static void agg_rows(const local_matrix_store &store, const agg_operate &op,
		local_matrix_store &res)
{
	size_t ncol = store.get_num_cols();
	size_t nrow = store.get_num_rows();
	if (ncol == 1) {
		// When we run agg on a single element, we assume we get the same value
		// as the input.
		res.copy_from(store);
		return;
	}

	// We always assume this is a single-column matrix.
	assert(res.get_num_cols() == 1);
	// Aggregate on rows, but the matrix is stored in col-major.
	// Instead of running aggregation on each row directly, we compute partial
	// aggregation on columns. This only works if the aggregation operation
	// allows partial aggregation and the agg and combine operations are
	// the same.
	if (store.store_layout() == matrix_layout_t::L_COL && op.is_same()) {
		size_t res_len = res.get_num_rows();
		assert(res_len == store.get_num_rows());
		assert(res.store_layout() == matrix_layout_t::L_COL);
		local_col_matrix_store &col_res
			= static_cast<local_col_matrix_store &>(res);
		char *res_arr = col_res.get_col(0);
		assert(res_arr);
		const local_col_matrix_store &col_store
			= static_cast<const local_col_matrix_store &>(store);
		op.get_agg().runAA(res_len, col_store.get_col(0),
				col_store.get_col(1), res_arr);
		for (size_t i = 2; i < ncol; i++)
			op.get_agg().runAA(res_len, res_arr, col_store.get_col(i),
					res_arr);
	}
	// This is the default solution to compute aggregation.
	else {
		local_matrix_store::const_ptr buf_mat;
		assert(res.get_num_rows() == store.get_num_rows());
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
			op.runAgg(ncol, row_store->get_row(i), res.get(i, 0));
	}
}

// This is implemented very similarly to above.
static void agg_cols(const local_matrix_store &store, const agg_operate &op,
		local_matrix_store &res)
{
	size_t ncol = store.get_num_cols();
	size_t nrow = store.get_num_rows();
	if (nrow == 1) {
		// When we run agg on a single element, we assume we get the same value
		// as the input.
		res.copy_from(store);
		return;
	}

	// We always assume this is a single-column matrix.
	assert(res.get_num_cols() == 1);
	if (store.store_layout() == matrix_layout_t::L_ROW && op.is_same()) {
		size_t res_len = res.get_num_rows();
		assert(res_len == store.get_num_cols());
		assert(res.store_layout() == matrix_layout_t::L_COL);
		local_col_matrix_store &col_res
			= static_cast<local_col_matrix_store &>(res);
		char *res_arr = col_res.get_col(0);
		assert(res_arr);
		const local_row_matrix_store &row_store
			= static_cast<const local_row_matrix_store &>(store);
		op.get_agg().runAA(res_len, row_store.get_row(0),
				row_store.get_row(1), res_arr);
		for (size_t i = 2; i < nrow; i++)
			op.get_agg().runAA(res_len, res_arr, row_store.get_row(i),
					res_arr);
	}
	else {
		local_matrix_store::const_ptr buf_mat;
		assert(res.get_num_rows() == store.get_num_cols());
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
			op.runAgg(nrow, col_store->get_col(i), res.get(i, 0));
	}
}

static void _aggregate(const local_matrix_store &store, const agg_operate &op,
		int margin, local_matrix_store &res)
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
	store.complete();
}

void aggregate(const local_matrix_store &store, const agg_operate &op,
		int margin, part_dim_t dim, local_matrix_store &res)
{
	if (dim == part_dim_t::PART_NONE) {
		_aggregate(store, op, margin, res);
		return;
	}

	const size_t part_len = get_part_dim_len(store, dim);
	bool agg_long = margin == matrix_margin::BOTH
		|| (margin == matrix_margin::MAR_ROW && dim == part_dim_t::PART_DIM2)
		|| (margin == matrix_margin::MAR_COL && dim == part_dim_t::PART_DIM1);
	if ((dim == part_dim_t::PART_DIM1 && part_len >= store.get_num_rows())
			|| (dim == part_dim_t::PART_DIM2 && part_len >= store.get_num_cols())
			|| (!op.is_same() && agg_long)
			// If the input matrix has been resized.
			|| !store.is_whole())
		_aggregate(store, op, margin, res);
	else {
		local_matrix_store::exposed_area orig_in = store.get_exposed_area();
		local_matrix_store::exposed_area orig_res;
		bool resize_res = false;
		size_t lnum_rows;
		// We only need to resize the res matrix in the following two conditions.
		if ((margin == matrix_margin::MAR_COL && dim == part_dim_t::PART_DIM2)
				|| (margin == matrix_margin::MAR_ROW
					&& dim == part_dim_t::PART_DIM1)) {
			resize_res = true;
			orig_res = res.get_exposed_area();
			lnum_rows = 0;
		}
		else
			lnum_rows = res.get_num_rows();
		local_buf_col_matrix_store lbuf(0, 0, lnum_rows, 1, res.get_type(), -1);
		assert(res.store_layout() == matrix_layout_t::L_COL);
		local_col_matrix_store &col_res
			= static_cast<local_col_matrix_store &>(res);
		local_matrix_store &mutable_store
			= const_cast<local_matrix_store &>(store);
		// The res matrix is a single-column matrix.
		assert(res.get_num_cols() == 1);
		bool success = true;
		if (dim == part_dim_t::PART_DIM2) {
			size_t orig_num_cols = store.get_num_cols();
			for (size_t col_idx = 0; col_idx < orig_num_cols;
					col_idx += part_len) {
				size_t llen = std::min(orig_num_cols - col_idx, part_len);
				success = mutable_store.resize(orig_in.local_start_row,
						orig_in.local_start_col + col_idx,
						store.get_num_rows(), llen);
				if (!success) {
					assert(col_idx == 0);
					break;
				}

				if (resize_res) {
					res.resize(col_idx, 0, llen, 1);
					_aggregate(store, op, margin, res);
				}
				else if (col_idx == 0)
					_aggregate(store, op, margin, res);
				else {
					_aggregate(store, op, margin, lbuf);
					op.get_agg().runAA(res.get_num_rows(), lbuf.get_raw_arr(),
							col_res.get_col(0), col_res.get_col(0));
				}
			}
		}
		else {
			size_t orig_num_rows = store.get_num_rows();
			for (size_t row_idx = 0; row_idx < orig_num_rows;
					row_idx += part_len) {
				size_t llen = std::min(orig_num_rows - row_idx, part_len);
				success = mutable_store.resize(orig_in.local_start_row + row_idx,
						orig_in.local_start_col, llen, store.get_num_cols());
				if (!success) {
					assert(row_idx == 0);
					break;
				}

				if (resize_res) {
					res.resize(row_idx, 0, llen, 1);
					_aggregate(store, op, margin, res);
				}
				else if (row_idx == 0)
					_aggregate(store, op, margin, res);
				else {
					_aggregate(store, op, margin, lbuf);
					op.get_agg().runAA(res.get_num_rows(), lbuf.get_raw_arr(),
							col_res.get_col(0), col_res.get_col(0));
				}
			}
		}
		mutable_store.restore_size(orig_in);
		if (resize_res)
			res.restore_size(orig_res);
		if (!success)
			_aggregate(store, op, margin, res);
	}
}

void copy_last_col(const local_matrix_store &store, local_vec_store &vec)
{
	assert(store.get_num_rows() == vec.get_length());
	if (store.store_layout() == matrix_layout_t::L_COL) {
		const local_col_matrix_store &col_store
			= static_cast<const local_col_matrix_store &>(store);
		memcpy(vec.get_raw_arr(), col_store.get_col(store.get_num_cols() - 1),
				store.get_num_rows() * store.get_type().get_size());
	}
	else {
		size_t entry_size = store.get_type().get_size();
		for (size_t i = 0; i < store.get_num_rows(); i++)
			memcpy(vec.get(i), store.get(i, store.get_num_cols() - 1), entry_size);
	}
}

void copy_last_row(const local_matrix_store &store, local_vec_store &vec)
{
	assert(store.get_num_cols() == vec.get_length());
	if (store.store_layout() == matrix_layout_t::L_ROW) {
		const local_row_matrix_store &row_store
			= static_cast<const local_row_matrix_store &>(store);
		memcpy(vec.get_raw_arr(), row_store.get_row(store.get_num_rows() - 1),
				store.get_num_cols() * store.get_type().get_size());
	}
	else {
		size_t entry_size = store.get_type().get_size();
		for (size_t i = 0; i < store.get_num_cols(); i++)
			memcpy(vec.get(i), store.get(store.get_num_rows() - 1, i), entry_size);
	}
}

static void _cum(const local_matrix_store &store, const local_vec_store *prev_res,
		const agg_operate &op, int margin, local_matrix_store &res)
{
	assert(store.store_layout() == res.store_layout());
	assert(store.get_num_rows() == res.get_num_rows());
	assert(store.get_num_cols() == res.get_num_cols());
	if (prev_res)
		assert((prev_res->get_length() == store.get_num_rows()
					&& margin == matrix_margin::MAR_ROW)
				|| (prev_res->get_length() == store.get_num_cols()
					&& margin == matrix_margin::MAR_COL));
	// If accumulate on rows and the matrix is stored in row major.
	if (margin == matrix_margin::MAR_ROW
			&& store.store_layout() == matrix_layout_t::L_ROW) {
		const local_row_matrix_store &row_store
			= static_cast<const local_row_matrix_store &>(store);
		local_row_matrix_store &row_res
			= static_cast<local_row_matrix_store &>(res);
		for (size_t i = 0; i < store.get_num_rows(); i++) {
			const char *prev = NULL;
			if (prev_res)
				prev = prev_res->get(i);
			op.get_agg().runCum(row_store.get_num_cols(),
					row_store.get_row(i), prev, row_res.get_row(i));
		}
	}
	// If accumulate on rows and the matrix is stored in col major.
	else if (margin == matrix_margin::MAR_ROW
			&& store.store_layout() == matrix_layout_t::L_COL) {
		const local_col_matrix_store &col_store
			= static_cast<const local_col_matrix_store &>(store);
		local_col_matrix_store &col_res
			= static_cast<local_col_matrix_store &>(res);
		if (prev_res)
			op.get_agg().runAA(store.get_num_rows(), prev_res->get_raw_arr(),
					col_store.get_col(0), col_res.get_col(0));
		else
			memcpy(col_res.get_col(0), col_store.get_col(0),
					store.get_num_rows() * store.get_type().get_size());
		for (size_t i = 1; i < store.get_num_cols(); i++)
			op.get_agg().runAA(store.get_num_rows(), col_store.get_col(i),
					col_res.get_col(i - 1), col_res.get_col(i));
	}
	// If accumulate on cols and the matrix is stored in row major.
	else if (margin == matrix_margin::MAR_COL
			&& store.store_layout() == matrix_layout_t::L_ROW) {
		const local_row_matrix_store &row_store
			= static_cast<const local_row_matrix_store &>(store);
		local_row_matrix_store &row_res
			= static_cast<local_row_matrix_store &>(res);
		if (prev_res)
			op.get_agg().runAA(store.get_num_cols(), prev_res->get_raw_arr(),
					row_store.get_row(0), row_res.get_row(0));
		else
			memcpy(row_res.get_row(0), row_store.get_row(0),
					store.get_num_cols() * store.get_type().get_size());
		for (size_t i = 1; i < store.get_num_rows(); i++)
			op.get_agg().runAA(store.get_num_cols(), row_store.get_row(i),
					row_res.get_row(i - 1), row_res.get_row(i));
	}
	// If accumulate on cols and the matrix is stored in col major.
	else {
		const local_col_matrix_store &col_store
			= static_cast<const local_col_matrix_store &>(store);
		local_col_matrix_store &col_res
			= static_cast<local_col_matrix_store &>(res);
		for (size_t i = 0; i < store.get_num_cols(); i++) {
			const char *prev = NULL;
			if (prev_res)
				prev = prev_res->get(i);
			op.get_agg().runCum(col_store.get_num_rows(),
					col_store.get_col(i), prev, col_res.get_col(i));
		}
	}
}

void cum(const local_matrix_store &store, const local_vec_store *prev_res,
		const agg_operate &op, int margin, part_dim_t dim,
		local_matrix_store &res)
{
	const size_t part_len = get_part_dim_len(store, dim);
	local_vec_store::ptr prev_res_buf;
	// resize the wide matrix.
	if (store.is_virtual() && dim == part_dim_t::PART_DIM2
			&& store.get_num_cols() > part_len
			// both input matrices haven't been resized.
			&& store.is_whole()) {
		if (prev_res) {
			prev_res_buf = local_vec_store::ptr(new local_buf_vec_store(0,
						prev_res->get_length(), prev_res->get_type(), -1));
			memcpy(prev_res_buf->get_raw_arr(), prev_res->get_raw_arr(),
					prev_res->get_length() * prev_res->get_type().get_size());
		}
		size_t orig_num_cols = store.get_num_cols();
		local_matrix_store::exposed_area orig_store = store.get_exposed_area();
		local_matrix_store::exposed_area orig_res = res.get_exposed_area();
		local_matrix_store &mutable_store = const_cast<local_matrix_store &>(store);
		bool success = true;
		for (size_t col_idx = 0; col_idx < orig_num_cols;
				col_idx += part_len) {
			size_t llen = std::min(orig_num_cols - col_idx, part_len);
			success = mutable_store.resize(orig_store.local_start_row,
					orig_store.local_start_col + col_idx, store.get_num_rows(), llen);
			if (success)
				success = res.resize(orig_res.local_start_row,
						orig_res.local_start_col + col_idx, res.get_num_rows(),
						llen);
			if (!success) {
				assert(col_idx == 0);
				break;
			}
			_cum(store, prev_res_buf.get(), op, margin, res);
			// If we perform cumulative computation on rows, we need to
			// get the last col from the previous computation.
			if (margin == matrix_margin::MAR_ROW) {
				if (prev_res_buf == NULL)
					prev_res_buf = local_vec_store::ptr(new local_buf_vec_store(
								0, store.get_num_rows(), store.get_type(), -1));
				copy_last_col(res, *prev_res_buf);
			}
		}
		mutable_store.restore_size(orig_store);
		res.restore_size(orig_res);
		if (!success)
			_cum(store, prev_res, op, margin, res);
	}
	// resize the tall matrix
	else if (store.is_virtual() && dim == part_dim_t::PART_DIM1
			&& store.get_num_rows() > part_len
			// both input matrices haven't been resized.
			&& store.is_whole()) {
		if (prev_res) {
			prev_res_buf = local_vec_store::ptr(new local_buf_vec_store(0,
						prev_res->get_length(), prev_res->get_type(), -1));
			memcpy(prev_res_buf->get_raw_arr(), prev_res->get_raw_arr(),
					prev_res->get_length() * prev_res->get_type().get_size());
		}
		size_t orig_num_rows = store.get_num_rows();
		local_matrix_store::exposed_area orig_store = store.get_exposed_area();
		local_matrix_store::exposed_area orig_res = res.get_exposed_area();
		local_matrix_store &mutable_store = const_cast<local_matrix_store &>(store);
		bool success = true;
		for (size_t row_idx = 0; row_idx < orig_num_rows; row_idx += part_len) {
			size_t llen = std::min(orig_num_rows - row_idx, part_len);
			success = mutable_store.resize(orig_store.local_start_row + row_idx,
					orig_store.local_start_col, llen, store.get_num_cols());
			if (success)
				success = res.resize(orig_res.local_start_row + row_idx,
						orig_res.local_start_col, llen, res.get_num_cols());
			if (!success) {
				assert(row_idx == 0);
				break;
			}
			_cum(store, prev_res_buf.get(), op, margin, res);
			// If we perform cumulative computation on rows, we need to
			// get the last col from the previous computation.
			if (margin == matrix_margin::MAR_COL) {
				if (prev_res_buf == NULL)
					prev_res_buf = local_vec_store::ptr(new local_buf_vec_store(
								0, store.get_num_cols(), store.get_type(), -1));
				copy_last_row(res, *prev_res_buf);
			}
		}
		mutable_store.restore_size(orig_store);
		res.restore_size(orig_res);
		if (!success)
			_cum(store, prev_res, op, margin, res);
	}
	else {
		// If the local matrix isn't virtual, we don't need to resize it
		// to increase CPU cache hits.
		_cum(store, prev_res, op, margin, res);
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
	m1.complete();
	m2.complete();
}

void mapply2(const local_matrix_store &m1, const local_matrix_store &m2,
			const bulk_operate &op, part_dim_t dim, local_matrix_store &res)
{
	bool is_virt = m1.is_virtual() || m2.is_virtual();
	const size_t part_len = get_part_dim_len(m1, dim);
	// resize the wide matrix.
	if (is_virt && dim == part_dim_t::PART_DIM2 && m1.get_num_cols() > part_len
			// both input matrices haven't been resized.
			&& m1.is_whole() && m2.is_whole()) {
		size_t orig_num_cols = m1.get_num_cols();
		local_matrix_store::exposed_area orig_m1 = m1.get_exposed_area();
		local_matrix_store::exposed_area orig_m2 = m2.get_exposed_area();
		local_matrix_store::exposed_area orig_res = res.get_exposed_area();
		local_matrix_store &mutable_m1 = const_cast<local_matrix_store &>(m1);
		local_matrix_store &mutable_m2 = const_cast<local_matrix_store &>(m2);
		bool success = true;
		for (size_t col_idx = 0; col_idx < orig_num_cols;
				col_idx += part_len) {
			size_t llen = std::min(orig_num_cols - col_idx, part_len);
			success = mutable_m1.resize(orig_m1.local_start_row,
					orig_m1.local_start_col + col_idx, m1.get_num_rows(), llen);
			if (success)
				success = mutable_m2.resize(orig_m2.local_start_row,
						orig_m2.local_start_col + col_idx, m2.get_num_rows(),
						llen);
			if (success)
				success = res.resize(orig_res.local_start_row,
						orig_res.local_start_col + col_idx, res.get_num_rows(),
						llen);
			if (!success) {
				assert(col_idx == 0);
				break;
			}
			_mapply2(m1, m2, op, res);
		}
		mutable_m1.restore_size(orig_m1);
		mutable_m2.restore_size(orig_m2);
		res.restore_size(orig_res);
		if (!success)
			_mapply2(m1, m2, op, res);
	}
	// resize the tall matrix
	else if (is_virt && dim == part_dim_t::PART_DIM1
			&& m1.get_num_rows() > part_len
			// both input matrices haven't been resized.
			&& m1.is_whole() && m2.is_whole()) {
		size_t orig_num_rows = m1.get_num_rows();
		local_matrix_store::exposed_area orig_m1 = m1.get_exposed_area();
		local_matrix_store::exposed_area orig_m2 = m2.get_exposed_area();
		local_matrix_store::exposed_area orig_res = res.get_exposed_area();
		local_matrix_store &mutable_m1 = const_cast<local_matrix_store &>(m1);
		local_matrix_store &mutable_m2 = const_cast<local_matrix_store &>(m2);
		bool success = true;
		for (size_t row_idx = 0; row_idx < orig_num_rows; row_idx += part_len) {
			size_t llen = std::min(orig_num_rows - row_idx, part_len);
			success = mutable_m1.resize(orig_m1.local_start_row + row_idx,
					orig_m1.local_start_col, llen, m1.get_num_cols());
			if (success)
				success = mutable_m2.resize(orig_m2.local_start_row + row_idx,
						orig_m2.local_start_col, llen, m2.get_num_cols());
			if (success)
				success = res.resize(orig_res.local_start_row + row_idx,
						orig_res.local_start_col, llen, res.get_num_cols());
			if (!success) {
				assert(row_idx == 0);
				break;
			}
			_mapply2(m1, m2, op, res);
		}
		mutable_m1.restore_size(orig_m1);
		mutable_m2.restore_size(orig_m2);
		res.restore_size(orig_res);
		if (!success)
			_mapply2(m1, m2, op, res);
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
	store.complete();
}

void sapply(const local_matrix_store &store, const bulk_uoperate &op,
		part_dim_t dim, local_matrix_store &res)
{
	const size_t part_len = get_part_dim_len(store, dim);
	// resize the wide matrix.
	if (store.is_virtual() && dim == part_dim_t::PART_DIM2
			&& store.get_num_cols() > part_len
			// If the input matrix hasn't been resized.
			&& store.is_whole()) {
		size_t orig_num_cols = store.get_num_cols();
		local_matrix_store::exposed_area orig_in = store.get_exposed_area();
		local_matrix_store::exposed_area orig_res = res.get_exposed_area();
		local_matrix_store &mutable_store = const_cast<local_matrix_store &>(
				store);
		bool success = true;
		for (size_t col_idx = 0; col_idx < orig_num_cols; col_idx += part_len) {
			size_t llen = std::min(orig_num_cols - col_idx, part_len);
			success = mutable_store.resize(orig_in.local_start_row,
					orig_in.local_start_col + col_idx, store.get_num_rows(), llen);
			if (success)
				success = res.resize(orig_res.local_start_row,
						orig_res.local_start_col + col_idx, res.get_num_rows(), llen);
			if (!success) {
				assert(col_idx == 0);
				break;
			}
			_sapply(store, op, res);
		}
		mutable_store.restore_size(orig_in);
		res.restore_size(orig_res);
		if (!success)
			_sapply(store, op, res);
	}
	// resize the tall matrix
	else if (store.is_virtual() && dim == part_dim_t::PART_DIM1
			&& store.get_num_rows() > part_len
			// If the input matrix hasn't been resized.
			&& store.is_whole()) {
		size_t orig_num_rows = store.get_num_rows();
		local_matrix_store::exposed_area orig_in = store.get_exposed_area();
		local_matrix_store::exposed_area orig_res = res.get_exposed_area();
		local_matrix_store &mutable_store = const_cast<local_matrix_store &>(
				store);
		bool success = true;
		for (size_t row_idx = 0; row_idx < orig_num_rows; row_idx += part_len) {
			size_t llen = std::min(orig_num_rows - row_idx, part_len);
			success = mutable_store.resize(orig_in.local_start_row + row_idx,
					orig_in.local_start_col, llen, store.get_num_cols());
			if (success)
				success = res.resize(orig_res.local_start_row + row_idx,
						orig_res.local_start_col, llen, res.get_num_cols());
			if (!success) {
				assert(row_idx == 0);
				break;
			}
			_sapply(store, op, res);
		}
		mutable_store.restore_size(orig_in);
		res.restore_size(orig_res);
		if (!success)
			_sapply(store, op, res);
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
	in_mat.complete();
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
			void *output) const {
		assert(0);
	}
	virtual void runCum(size_t num_eles, const void *left_arr,
			const void *prev, void *output) const {
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
	store.complete();
	// TODO we should split the matrix
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
	store.complete();
	// TODO we should split the matrix
}

static bool _groupby(const detail::local_matrix_store &labels,
		const detail::local_matrix_store &mat, const agg_operate &op,
		matrix_margin margin, detail::local_matrix_store &results,
		std::vector<bool> &agg_flags)
{
	if (margin == matrix_margin::MAR_ROW) {
		assert(mat.store_layout() == matrix_layout_t::L_ROW);
		assert(results.store_layout() == matrix_layout_t::L_ROW);
		assert(labels.get_num_rows() == mat.get_num_rows());
		const detail::local_row_matrix_store &row_mat
			= static_cast<const detail::local_row_matrix_store &>(mat);
		detail::local_row_matrix_store &row_results
			= static_cast<detail::local_row_matrix_store &>(results);
		for (size_t i = 0; i < row_mat.get_num_rows(); i++) {
			factor_value_t label_id = labels.get<factor_value_t>(i, 0);
			if ((size_t) label_id >= agg_flags.size()) {
				BOOST_LOG_TRIVIAL(error) << boost::format(
						"Factor value %1% is larger than max levels %2%")
					% (size_t) label_id % agg_flags.size();
				return false;
			}
			// If we never get partially aggregated result for a label, we should
			// copy the data to the corresponding row.
			if (!agg_flags[label_id])
				memcpy(row_results.get_row(label_id), row_mat.get_row(i),
						row_mat.get_num_cols() * row_mat.get_entry_size());
			else
				op.get_agg().runAA(row_mat.get_num_cols(), row_mat.get_row(i),
						row_results.get_row(label_id),
						row_results.get_row(label_id));
			auto bit = agg_flags[label_id];
			bit = true;
		}
	}
	else {
		assert(mat.store_layout() == matrix_layout_t::L_COL);
		assert(results.store_layout() == matrix_layout_t::L_COL);
		assert(labels.get_num_cols() == mat.get_num_cols());
		const detail::local_col_matrix_store &col_mat
			= static_cast<const detail::local_col_matrix_store &>(mat);
		detail::local_col_matrix_store &col_results
			= static_cast<detail::local_col_matrix_store &>(results);
		for (size_t i = 0; i < col_mat.get_num_cols(); i++) {
			factor_value_t label_id = labels.get<factor_value_t>(0, i);
			if ((size_t) label_id >= agg_flags.size()) {
				BOOST_LOG_TRIVIAL(error) << boost::format(
						"Factor value %1% is larger than max levels %2%")
					% (size_t) label_id % agg_flags.size();
				return false;
			}
			// If we never get partially aggregated result for a label, we should
			// copy the data to the corresponding col.
			if (!agg_flags[label_id])
				memcpy(col_results.get_col(label_id), col_mat.get_col(i),
						col_mat.get_num_rows() * col_mat.get_entry_size());
			else
				op.get_agg().runAA(col_mat.get_num_rows(), col_mat.get_col(i),
						col_results.get_col(label_id),
						col_results.get_col(label_id));
			auto bit = agg_flags[label_id];
			bit = true;
		}
	}
	labels.complete();
	mat.complete();
	return true;
}

bool groupby(const detail::local_matrix_store &labels,
		const detail::local_matrix_store &mat, const agg_operate &op,
		matrix_margin margin, part_dim_t dim,
		detail::local_matrix_store &results, std::vector<bool> &agg_flags)
{
	const size_t part_len = get_part_dim_len(mat, dim);
	// Group by rows on a wide matrix.
	if (dim == part_dim_t::PART_DIM2 && margin == matrix_margin::MAR_ROW
			&& mat.get_num_cols() > part_len
			// the input matrix hasn't been resized.
			&& mat.is_whole()) {
		// We need to resize the input matrix and the result matrix.
		size_t orig_num_cols = mat.get_num_cols();
		local_matrix_store::exposed_area orig_in = mat.get_exposed_area();
		local_matrix_store::exposed_area orig_res = results.get_exposed_area();
		local_matrix_store &mutable_mat = const_cast<local_matrix_store &>(mat);
		bool success = true;
		for (size_t col_idx = 0; col_idx < orig_num_cols; col_idx += part_len) {
			size_t llen = std::min(orig_num_cols - col_idx, part_len);
			success = mutable_mat.resize(orig_in.local_start_row,
					orig_in.local_start_col + col_idx, mat.get_num_rows(), llen);
			if (success)
				success = results.resize(orig_res.local_start_row,
						orig_res.local_start_col + col_idx,
						results.get_num_rows(), llen);
			if (!success) {
				assert(col_idx == 0);
				break;
			}
			// We need to reset the flags when we compute on smaller parts.
			size_t orig_num = agg_flags.size();
			agg_flags.clear();
			agg_flags.resize(orig_num, false);
			bool ret = _groupby(labels, mat, op, margin, results, agg_flags);
			if (!ret) {
				mutable_mat.restore_size(orig_in);
				results.restore_size(orig_res);
				return false;
			}
		}
		mutable_mat.restore_size(orig_in);
		results.restore_size(orig_res);
		if (!success)
			return _groupby(labels, mat, op, margin, results, agg_flags);
		else
			return true;
	}
	// Group by cols on a tall matrix.
	else if (dim == part_dim_t::PART_DIM1 && margin == matrix_margin::MAR_COL
			&& mat.get_num_rows() > part_len
			// the input matrix hasn't been resized.
			&& mat.is_whole()) {
		size_t orig_num_rows = mat.get_num_rows();
		local_matrix_store::exposed_area orig_in = mat.get_exposed_area();
		local_matrix_store::exposed_area orig_res = results.get_exposed_area();
		local_matrix_store &mutable_mat = const_cast<local_matrix_store &>(mat);
		bool success = true;
		// We need to resize the input matrix and the result matrix.
		for (size_t row_idx = 0; row_idx < orig_num_rows; row_idx += part_len) {
			size_t llen = std::min(orig_num_rows - row_idx, part_len);
			success = mutable_mat.resize(orig_in.local_start_row + row_idx,
					orig_in.local_start_col, llen, mat.get_num_cols());
			if (success)
				success = results.resize(orig_res.local_start_row + row_idx,
						orig_res.local_start_col, llen, results.get_num_cols());
			if (!success) {
				assert(row_idx == 0);
				break;
			}
			// We need to reset the flags when we compute on smaller parts.
			size_t orig_num = agg_flags.size();
			agg_flags.clear();
			agg_flags.resize(orig_num, false);
			bool ret = _groupby(labels, mat, op, margin, results, agg_flags);
			if (!ret) {
				mutable_mat.restore_size(orig_in);
				results.restore_size(orig_res);
				return false;
			}
		}
		mutable_mat.restore_size(orig_in);
		results.restore_size(orig_res);
		if (!success)
			return _groupby(labels, mat, op, margin, results, agg_flags);
		else
			return true;
	}
	// Group by rows on a tall matrix.
	else if (mat.is_virtual() && dim == part_dim_t::PART_DIM1
			&& margin == matrix_margin::MAR_ROW && mat.get_num_rows() > part_len
			// the input matrix hasn't been resized.
			&& mat.is_whole()) {
		// Here we need to resize the input matrix and the label vector.
		size_t orig_num_rows = mat.get_num_rows();
		local_matrix_store::exposed_area orig_in = mat.get_exposed_area();
		local_matrix_store::exposed_area orig_labels = labels.get_exposed_area();
		local_matrix_store &mutable_mat = const_cast<local_matrix_store &>(mat);
		local_matrix_store &mutable_labels = const_cast<local_matrix_store &>(
				labels);
		bool success = true;
		for (size_t row_idx = 0; row_idx < orig_num_rows; row_idx += part_len) {
			size_t llen = std::min(orig_num_rows - row_idx, part_len);
			success = mutable_mat.resize(orig_in.local_start_row + row_idx,
					orig_in.local_start_col, llen, mat.get_num_cols());
			if (success)
				success = mutable_labels.resize(orig_labels.local_start_row + row_idx,
						orig_labels.local_start_col, llen, labels.get_num_cols());
			if (!success) {
				assert(row_idx == 0);
				break;
			}
			bool ret = _groupby(labels, mat, op, margin, results, agg_flags);
			if (!ret) {
				mutable_mat.restore_size(orig_in);
				mutable_labels.restore_size(orig_labels);
				return false;
			}
		}
		mutable_mat.restore_size(orig_in);
		mutable_labels.restore_size(orig_labels);
		if (!success)
			return _groupby(labels, mat, op, margin, results, agg_flags);
		else
			return true;
	}
	// Group by cols on a wide matrix.
	else if (mat.is_virtual() && dim == part_dim_t::PART_DIM2
			&& margin == matrix_margin::MAR_COL && mat.get_num_cols() > part_len
			// the input matrix hasn't been resized.
			&& mat.is_whole()) {
		// Here we need to resize the input matrix and the label vector.
		size_t orig_num_cols = mat.get_num_cols();
		local_matrix_store::exposed_area orig_in = mat.get_exposed_area();
		local_matrix_store::exposed_area orig_labels = labels.get_exposed_area();
		local_matrix_store &mutable_mat = const_cast<local_matrix_store &>(mat);
		local_matrix_store &mutable_labels = const_cast<local_matrix_store &>(
				labels);
		bool success = true;
		for (size_t col_idx = 0; col_idx < orig_num_cols; col_idx += part_len) {
			size_t llen = std::min(orig_num_cols - col_idx, part_len);
			success = mutable_mat.resize(orig_in.local_start_row,
					orig_in.local_start_col + col_idx, mat.get_num_rows(), llen);
			if (success)
				success = mutable_labels.resize(orig_labels.local_start_row,
						orig_labels.local_start_col + col_idx,
						labels.get_num_rows(), llen);
			if (!success) {
				assert(col_idx == 0);
				break;
			}
			bool ret = _groupby(labels, mat, op, margin, results, agg_flags);
			if (!ret) {
				mutable_mat.restore_size(orig_in);
				results.restore_size(orig_labels);
				return false;
			}
		}
		mutable_mat.restore_size(orig_in);
		mutable_labels.restore_size(orig_labels);
		if (!success)
			return _groupby(labels, mat, op, margin, results, agg_flags);
		else
			return true;
	}
	else
		return _groupby(labels, mat, op, margin, results, agg_flags);
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

///////////////////////// BLAS matrix multiplication //////////////////////////

template<class T>
void tall_gemm_col(const local_matrix_store &Astore, const T *Amat,
		const local_matrix_store &Bstore, const T *Bmat,
		local_matrix_store &out, T *res_mat)
{
	assert(0);
}

template<>
void tall_gemm_col<double>(const local_matrix_store &Astore, const double *Amat,
		const local_matrix_store &Bstore, const double *Bmat,
		local_matrix_store &out, double *res_mat)
{
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			Astore.get_num_rows(), Bstore.get_num_cols(),
			Astore.get_num_cols(), 1, Amat,
			Astore.get_num_rows(), Bmat, Bstore.get_num_rows(),
			0, res_mat, out.get_num_rows());
}

template<>
void tall_gemm_col<float>(const local_matrix_store &Astore, const float *Amat,
		const local_matrix_store &Bstore, const float *Bmat,
		local_matrix_store &out, float *res_mat)
{
	cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			Astore.get_num_rows(), Bstore.get_num_cols(),
			Astore.get_num_cols(), 1, Amat,
			Astore.get_num_rows(), Bmat, Bstore.get_num_rows(),
			0, res_mat, out.get_num_rows());
}

template<class T>
void tall_gemm_row(const local_matrix_store &Astore, const T *Amat,
		const local_matrix_store &Bstore, const T *Bmat,
		local_matrix_store &out, T *res_mat)
{
	assert(0);
}

template<>
void tall_gemm_row<double>(const local_matrix_store &Astore, const double *Amat,
		const local_matrix_store &Bstore, const double *Bmat,
		local_matrix_store &out, double *res_mat)
{
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			Astore.get_num_rows(), Bstore.get_num_cols(),
			Astore.get_num_cols(), 1, Amat,
			Astore.get_num_cols(), Bmat, Bstore.get_num_cols(),
			0, res_mat, out.get_num_cols());
}

template<>
void tall_gemm_row<float>(const local_matrix_store &Astore, const float *Amat,
		const local_matrix_store &Bstore, const float *Bmat,
		local_matrix_store &out, float *res_mat)
{
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			Astore.get_num_rows(), Bstore.get_num_cols(),
			Astore.get_num_cols(), 1, Amat,
			Astore.get_num_cols(), Bmat, Bstore.get_num_cols(),
			0, res_mat, out.get_num_cols());
}

static local_matrix_store::ptr alloc_mat_buf(size_t num_rows, size_t num_cols,
		const scalar_type &type, matrix_layout_t layout)
{
	if (layout == matrix_layout_t::L_COL)
		return local_matrix_store::ptr(new local_buf_col_matrix_store(0, 0,
					num_rows, num_cols, type, -1));
	else
		return local_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
					num_rows, num_cols, type, -1));
}

template<class T>
void _matrix_tall_multiply(const local_matrix_store &Astore,
		const local_matrix_store &Bstore, local_matrix_store &out,
		std::pair<local_matrix_store::ptr, local_matrix_store::ptr> &bufs)
{
	const T *Amat = (const T *) Astore.get_raw_arr();
	// If the A matrix isn't stored in contiguous memory,
	// and the matrix buffer doesn't exist or has a wrong size.
	if (Amat == NULL || Astore.store_layout() != out.store_layout()) {
		if (bufs.first == NULL
				|| bufs.first->get_num_rows() != Astore.get_num_rows())
			// Let's make sure all matrices have the same data layout as
			// the result matrix.
			bufs.first = alloc_mat_buf(Astore.get_num_rows(),
					Astore.get_num_cols(), Astore.get_type(), out.store_layout());
		bufs.first->copy_from(Astore);
		Amat = reinterpret_cast<const T *>(bufs.first->get_raw_arr());
	}

	T *res_mat = (T *) out.get_raw_arr();
	// If the out matrix isn't stored in contiguous memory,
	// and the matrix buffer doesn't exist or has a wrong size.
	if (res_mat == NULL && (bufs.second == NULL
				|| bufs.second->get_num_rows() != out.get_num_rows()))
		bufs.second = alloc_mat_buf(out.get_num_rows(), out.get_num_cols(),
				out.get_type(), out.store_layout());
	// If the out matrix isn't stored in contiguous memory,
	if (res_mat == NULL)
		res_mat = reinterpret_cast<T *>(bufs.second->get_raw_arr());

	const T *Bmat = reinterpret_cast<const T *>(Bstore.get_raw_arr());
	assert(Bstore.store_layout() == out.store_layout());
	assert(Amat);
	assert(Bmat);
	assert(res_mat);

	if (out.store_layout() == matrix_layout_t::L_COL)
		tall_gemm_col<T>(Astore, Amat, Bstore, Bmat, out, res_mat);
	else
		tall_gemm_row<T>(Astore, Amat, Bstore, Bmat, out, res_mat);

	if (out.get_raw_arr() == NULL)
		out.copy_from(*bufs.second);
	Astore.complete();
	Bstore.complete();
}

void matrix_tall_multiply(const local_matrix_store &Astore,
		const local_matrix_store &Bstore, local_matrix_store &out,
		std::pair<local_matrix_store::ptr, local_matrix_store::ptr> &bufs)
{
	const size_t part_len = get_part_dim_len(Astore, part_dim_t::PART_DIM1);
	assert(Astore.get_type() == Bstore.get_type());
	assert(Astore.get_type() == out.get_type());
	// As long as Astore is taller than we expect, we partition it to compute
	// matrix multiplication. This may cause a little extra overhead if Astore
	// is row-major and isn't a virtual matrix.
	if (Astore.get_num_rows() > part_len
			// Both input matrices haven't been resized.
			&& Astore.is_whole() && Bstore.is_whole()) {
		size_t orig_num_rows = Astore.get_num_rows();
		local_matrix_store::exposed_area orig_A = Astore.get_exposed_area();
		local_matrix_store::exposed_area orig_out = out.get_exposed_area();
		local_matrix_store &mutableA = const_cast<local_matrix_store &>(Astore);
		bool success = true;
		for (size_t row_idx = 0; row_idx < orig_num_rows; row_idx += part_len) {
			size_t llen = std::min(orig_num_rows - row_idx, part_len);
			success = mutableA.resize(orig_A.local_start_row + row_idx,
					orig_A.local_start_col, llen, Astore.get_num_cols());
			if (success)
				success = out.resize(orig_out.local_start_row + row_idx,
						orig_out.local_start_col, llen, out.get_num_cols());
			if (!success)
				break;

			if (Astore.get_type() == get_scalar_type<double>())
				_matrix_tall_multiply<double>(Astore, Bstore, out, bufs);
			else
				_matrix_tall_multiply<float>(Astore, Bstore, out, bufs);
		}
		mutableA.restore_size(orig_A);
		out.restore_size(orig_out);
		if (!success) {
			if (Astore.get_type() == get_scalar_type<double>())
				_matrix_tall_multiply<double>(Astore, Bstore, out, bufs);
			else
				_matrix_tall_multiply<float>(Astore, Bstore, out, bufs);
		}
	}
	else {
		if (Astore.get_type() == get_scalar_type<double>())
			_matrix_tall_multiply<double>(Astore, Bstore, out, bufs);
		else
			_matrix_tall_multiply<float>(Astore, Bstore, out, bufs);
	}
}

void matrix_wide_multiply(const local_matrix_store &left,
		const local_matrix_store &right, local_matrix_store &out,
		std::pair<local_matrix_store::ptr, local_matrix_store::ptr> &bufs)
{
	// TODO we should split the matrix
}

void materialize_tall(
		const std::vector<detail::local_matrix_store::const_ptr> &ins)
{
	size_t orig_num_rows = ins[0]->get_num_rows();
	// We need all tall matrices have the same number of rows.
	for (size_t i = 1; i < ins.size(); i++)
		assert(ins[i]->get_num_rows() == orig_num_rows);

	// Most of the matrices materialized here are sink matrices.
	// When materializing these matrices, the corresponding computation
	// will futher break the matrices to increase CPU cache utilization.
	// For the best performance, we should let the computation to choose
	// the partition size. For example, matrix multiplication should
	// use a larger partition size to get performance from BLAS.
	for (size_t i = 0; i < ins.size(); i++)
		ins[i]->materialize_self();
}

void materialize_wide(
		const std::vector<detail::local_matrix_store::const_ptr> &ins)
{
	size_t orig_num_cols = ins[0]->get_num_cols();
	// We need all wide matrices have the same number of cols.
	for (size_t i = 1; i < ins.size(); i++)
		assert(ins[i]->get_num_cols() == orig_num_cols);

	// For the reason as above.
	for (size_t i = 0; i < ins.size(); i++)
		ins[i]->materialize_self();
}

size_t get_part_dim_len(const local_matrix_store &mat, part_dim_t dim)
{
	// In this case, the returned value shouldn't matter.
	if (dim == part_dim_t::PART_NONE)
		return std::max(mat.get_num_rows(), mat.get_num_cols());

	size_t other_dim
		= dim == part_dim_t::PART_DIM1 ? mat.get_num_cols() : mat.get_num_rows();

	size_t entry_size = mat.get_entry_size();
	size_t max_part_dim = L1_SIZE / other_dim / entry_size;
	// We use this minimal size for in the long dimension because normally,
	// the maximal size for the short dimension is 32.
	size_t part_dim = 128;
	// The maximal size for the long dimension is 1024.
	while (part_dim < max_part_dim && part_dim < 1024)
		part_dim *= 2;
	return part_dim;
}

size_t get_long_dim_len(const local_matrix_store &mat)
{
	return get_part_dim_len(mat,
			mat.is_wide() ? part_dim_t::PART_DIM2 : part_dim_t::PART_DIM1);
}

size_t get_long_dim_len(const local_matrix_store &mat1,
		const local_matrix_store &mat2)
{
	size_t short_dim1 = std::min(mat1.get_num_rows(), mat1.get_num_cols());
	size_t short_dim2 = std::min(mat2.get_num_rows(), mat2.get_num_cols());
	// We use the matrix with the larger length in the short dimension
	// to determine the length in the long dimension.
	if (short_dim1 < short_dim2)
		return get_long_dim_len(mat2);
	else
		return get_long_dim_len(mat1);
}

}

}
