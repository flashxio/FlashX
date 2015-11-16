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

#include <cblas.h>

#include "local_matrix_store.h"
#include "bulk_operate.h"
#include "dense_matrix.h"
#include "local_vec_store.h"

namespace fm
{

namespace detail
{

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
	if (get_raw_arr())
		memset(get_raw_arr(), 0, nrow * ncol * get_entry_size());
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
	for (size_t i = 0; i < nrow; i++)
		op.set(get_row(i), ncol, get_global_start_row() + i,
				get_global_start_col());
}

bool local_row_matrix_store::copy_from(const local_matrix_store &store)
{
	assert(this->get_num_rows() == store.get_num_rows());
	assert(this->get_num_cols() == store.get_num_cols());
	assert(store.get_type() == this->get_type());

	if (store.store_layout() == matrix_layout_t::L_ROW) {
		size_t ncol = get_num_cols();
		size_t nrow = get_num_rows();
		// If the store has data stored contiguously.
		if (get_raw_arr() && store.get_raw_arr())
			memcpy(get_raw_arr(), store.get_raw_arr(),
					nrow * ncol * get_entry_size());
		else {
			const local_row_matrix_store &row_store
				= static_cast<const local_row_matrix_store &>(store);
			for (size_t i = 0; i < nrow; i++)
				memcpy(get_row(i), row_store.get_row(i),
						ncol * get_entry_size());
		}
	}
	else {
		// TODO I should handle wide matrix and tall matrix differently
		// to minimize virtual function calls.
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
	return true;
}

void local_col_matrix_store::reset_data()
{
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	// If the store has data stored contiguously.
	if (get_raw_arr())
		memset(get_raw_arr(), 0, nrow * ncol * get_entry_size());
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
	for (size_t i = 0; i < ncol; i++)
		op.set(get_col(i), nrow, get_global_start_row(),
				get_global_start_col() + i);
}

bool local_col_matrix_store::copy_from(const local_matrix_store &store)
{
	assert(this->get_num_rows() == store.get_num_rows());
	assert(this->get_num_cols() == store.get_num_cols());
	assert(store.get_type() == this->get_type());

	if (store.store_layout() == matrix_layout_t::L_COL) {
		size_t ncol = get_num_cols();
		size_t nrow = get_num_rows();
		// If the store has data stored contiguously.
		if (get_raw_arr() && store.get_raw_arr())
			memcpy(get_raw_arr(), store.get_raw_arr(),
					nrow * ncol * get_entry_size());
		else {
			const local_col_matrix_store &col_store
				= static_cast<const local_col_matrix_store &>(store);
			for (size_t i = 0; i < ncol; i++)
				memcpy(get_col(i), col_store.get_col(i), nrow * get_entry_size());
		}
	}
	else {
		// TODO I should handle wide matrix and tall matrix differently
		// to minimize virtual function calls.
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
	return true;
}

namespace
{

const size_t SUB_CHUNK_SIZE = 1024;

/*
 * This class contains the information of a submatrix in the original matrix.
 * This is used to improve CPU cache hits.
 */
class sub_matrix_info
{
	size_t start_row;
	size_t start_col;
	size_t nrow;
	size_t ncol;
public:
	sub_matrix_info(size_t start_row, size_t nrow, size_t start_col,
			size_t ncol) {
		this->start_row = start_row;
		this->start_col = start_col;
		this->nrow = nrow;
		this->ncol = ncol;
	}

	size_t get_num_rows() const {
		return nrow;
	}

	size_t get_num_cols() const {
		return ncol;
	}

	size_t get_start_row() const {
		return start_row;
	}

	size_t get_start_col() const {
		return start_col;
	}
};

/*
 * This class contains the information of a submatrix in the a column-wise matrix.
 * It is mainly used in inner product.
 */
class sub_col_matrix_info: public sub_matrix_info
{
	const local_col_matrix_store &m;
public:
	sub_col_matrix_info(size_t start_row, size_t nrow, size_t start_col,
			size_t ncol, const local_col_matrix_store &_m): sub_matrix_info(
				start_row, nrow, start_col, ncol), m(_m) {
		assert(start_row + nrow <= m.get_num_rows());
		assert(start_col + ncol <= m.get_num_cols());
	}

	const char *get_col(size_t col) const {
		return m.get_col(get_start_col() + col)
			+ get_start_row() * m.get_entry_size();
	}
};

/*
 * This class contains the information of a submatrix in the a row-wise matrix.
 * It is mainly used in inner product.
 */
class sub_row_matrix_info: public sub_matrix_info
{
	const local_row_matrix_store &m;
public:
	sub_row_matrix_info(size_t start_row, size_t nrow, size_t start_col,
			size_t ncol, const local_row_matrix_store &_m): sub_matrix_info(
				start_row, nrow, start_col, ncol), m(_m) {
		assert(start_row + nrow <= m.get_num_rows());
		assert(start_col + ncol <= m.get_num_cols());
	}

	const char *get_row(size_t row) const {
		return m.get_row(get_start_row() + row) + get_start_col() * m.get_entry_size();
	}
};

/*
 * In this case, the left matrix is row-major and wide. The right matrix is tall and its data
 * is stored column-wise.
 */
void inner_prod_row_wide(const local_row_matrix_store &m1,
		const local_col_matrix_store &m2, const bulk_operate &left_op,
		const bulk_operate &right_op, local_row_matrix_store &res)
{
	size_t ncol = m1.get_num_cols();
	size_t nrow = m1.get_num_rows();
	std::unique_ptr<char[]> tmp_res(
			new char[SUB_CHUNK_SIZE * left_op.output_entry_size()]);
	std::unique_ptr<char[]> tmp_res2(
			new char[res.get_num_cols() * res.get_entry_size()]);
	for (size_t k = 0; k < ncol; k += SUB_CHUNK_SIZE) {
		size_t sub_ncol = std::min(SUB_CHUNK_SIZE, ncol - k);
		// We further split the input matrices into parts.
		sub_row_matrix_info sub_left(0, nrow, k, sub_ncol, m1);
		sub_col_matrix_info sub_right(k, sub_ncol, 0, m2.get_num_cols(), m2);
		for (size_t i = 0; i < sub_left.get_num_rows(); i++) {
			for (size_t j = 0; j < sub_right.get_num_cols(); j++) {
				left_op.runAA(sub_ncol, sub_left.get_row(i),
						sub_right.get_col(j), tmp_res.get());
				right_op.runAgg(sub_ncol, tmp_res.get(), NULL,
						tmp_res2.get() + res.get_entry_size() * j);
			}
			// This is fine because we assume the input type of the right operator
			// should be the same as the type of the output matrix.
			right_op.runAA(sub_right.get_num_cols(), tmp_res2.get(),
					res.get_row(i), res.get_row(i));
		}
	}
}

void inner_prod_col_wide(const local_col_matrix_store &m1,
		const local_col_matrix_store &m2, const bulk_operate &left_op,
		const bulk_operate &right_op, local_row_matrix_store &res)
{
	size_t ncol = m1.get_num_cols();
	size_t nrow = m1.get_num_rows();
	size_t left_entry_size = left_op.left_entry_size();
	std::unique_ptr<char[]> tmp_row(
			new char[SUB_CHUNK_SIZE * left_op.left_entry_size()]);
	std::unique_ptr<char[]> tmp_res(
			new char[SUB_CHUNK_SIZE * left_op.output_entry_size()]);
	std::unique_ptr<char[]> tmp_res2(
			new char[res.get_num_cols() * res.get_entry_size()]);
	std::vector<const char *> sub_left_row;
	const scatter_gather &sg = m1.get_type().get_sg();
	for (size_t k = 0; k < ncol; k += SUB_CHUNK_SIZE) {
		size_t sub_ncol = std::min(SUB_CHUNK_SIZE, ncol - k);

		// Get the pointers of the columns in the left matrix.
		sub_left_row.resize(sub_ncol);
		for (size_t i = 0; i < sub_ncol; i++)
			sub_left_row[i] = m1.get_col(i + k);

		sub_col_matrix_info sub_right(k, sub_ncol, 0, m2.get_num_cols(), m2);
		for (size_t i = 0; i < nrow; i++) {
			// Get the pointers to the elements in a row of the left matrix.
			if (i > 0) {
				for (size_t j = 0; j < sub_ncol; j++)
					sub_left_row[j] += left_entry_size;
			}
			sg.gather(sub_left_row, tmp_row.get());

			for (size_t j = 0; j < sub_right.get_num_cols(); j++) {
				left_op.runAA(sub_ncol, tmp_row.get(),
						sub_right.get_col(j), tmp_res.get());
				right_op.runAgg(sub_ncol, tmp_res.get(), NULL,
						tmp_res2.get() + res.get_entry_size() * j);
			}
			// This is fine because we assume the input type of the right operator
			// should be the same as the type of the output matrix.
			right_op.runAA(sub_right.get_num_cols(), tmp_res2.get(),
					res.get_row(i), res.get_row(i));
		}
	}
}

/*
 * In this case, the left matrix is tall and is stored in row major. I assume
 * the right matrix is small and is stored in column major. We don't need to consider
 * the case that the right matrix is wide because the product would
 * be too large to be stored in any storage media.
 */
void inner_prod_row_tall(const local_row_matrix_store &m1,
		const local_col_matrix_store &m2, const bulk_operate &left_op,
		const bulk_operate &right_op, local_row_matrix_store &res)
{
	size_t ncol = m1.get_num_cols();
	size_t nrow = m1.get_num_rows();
	char *tmp_res = (char *) malloc(ncol * res.get_entry_size());
	for (size_t i = 0; i < nrow; i++) {
		for (size_t j = 0; j < m2.get_num_cols(); j++) {
			left_op.runAA(ncol, m1.get_row(i), m2.get_col(j), tmp_res);
			right_op.runAgg(ncol, tmp_res, NULL, res.get(i, j));
		}
	}
	free(tmp_res);
}

/*
 * In this case, the left matrix is tall and stored in column major. I assume
 * the right matrix is small and don't care its format.
 */
void inner_prod_col_tall(const local_col_matrix_store &m1,
		const local_matrix_store &m2, const bulk_operate &left_op,
		const bulk_operate &right_op, local_col_matrix_store &res)
{
	size_t ncol = m1.get_num_cols();
	size_t nrow = m1.get_num_rows();
	char *tmp_res = (char *) malloc(SUB_CHUNK_SIZE * res.get_entry_size());
	// We further break the local matrix into small matrices to increase
	// CPU cache hits.
	for (size_t k = 0; k < nrow; k += SUB_CHUNK_SIZE) {
		sub_col_matrix_info subm(k, std::min(SUB_CHUNK_SIZE, nrow - k),
				0, ncol, m1);
		for (size_t i = 0; i < ncol; i++) {
			for (size_t j = 0; j < m2.get_num_cols(); j++) {
				left_op.runAE(subm.get_num_rows(), subm.get_col(i),
						m2.get(i, j), tmp_res);
				char *store_col = res.get_col(j) + k * res.get_entry_size();
				right_op.runAA(subm.get_num_rows(), tmp_res, store_col,
						store_col);
			}
		}
	}
	free(tmp_res);
}

}

void aggregate(const local_matrix_store &store, const bulk_operate &op,
		agg_margin margin, local_vec_store &res)
{
	size_t output_size = op.output_entry_size();
	size_t ncol = store.get_num_cols();
	size_t nrow = store.get_num_rows();
	if (margin == agg_margin::BOTH) {
		assert(res.get_length() == 1);
		// If the store has data stored contiguously.
		if (store.get_raw_arr())
			op.runAgg(ncol * nrow, store.get_raw_arr(), NULL, res.get_raw_arr());
		// For row-major matrix.
		else if (store.store_layout() == matrix_layout_t::L_ROW) {
			const local_row_matrix_store &row_store
				= static_cast<const local_row_matrix_store &>(store);
			// We need to initialize the result first.
			op.runAgg(ncol, row_store.get_row(0), NULL, res.get_raw_arr());
			for (size_t i = 1; i < nrow; i++)
				op.runAgg(ncol, row_store.get_row(i), res.get_raw_arr(),
						res.get_raw_arr());
		}
		else {
			assert(store.store_layout() == matrix_layout_t::L_COL);
			const local_col_matrix_store &col_store
				= static_cast<const local_col_matrix_store &>(store);
			// We need to initialize the result first.
			op.runAgg(nrow, col_store.get_col(0), NULL, res.get_raw_arr());
			for (size_t i = 1; i < ncol; i++)
				op.runAgg(nrow, col_store.get_col(i), res.get_raw_arr(),
						res.get_raw_arr());
		}
	}
	else if (margin == agg_margin::MAR_ROW) {
		local_matrix_store::const_ptr buf_mat;
		const local_row_matrix_store *row_store;
		if (store.store_layout() == matrix_layout_t::L_COL) {
			buf_mat = store.conv2(matrix_layout_t::L_ROW);
			assert(buf_mat);
			row_store = static_cast<const local_row_matrix_store *>(
					buf_mat.get());
		}
		else
			row_store = static_cast<const local_row_matrix_store *>(&store);
		assert(res.get_length() == store.get_num_rows());
		for (size_t i = 0; i < nrow; i++)
			op.runAgg(ncol, row_store->get_row(i), NULL, res.get(i));
	}
	else if (margin == agg_margin::MAR_COL) {
		local_matrix_store::const_ptr buf_mat;
		const local_col_matrix_store *col_store;
		if (store.store_layout() == matrix_layout_t::L_ROW) {
			buf_mat = store.conv2(matrix_layout_t::L_COL);
			col_store = static_cast<const local_col_matrix_store *>(
					buf_mat.get());
		}
		else
			col_store = static_cast<const local_col_matrix_store *>(&store);
		assert(res.get_length() == store.get_num_cols());
		for (size_t i = 0; i < ncol; i++)
			op.runAgg(nrow, col_store->get_col(i), NULL, res.get(i));
	}
	else {
		// This shouldn't happen.
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"aggregate on an unknown margin %1%") % margin;
		assert(0);
	}
}

void mapply2(const local_matrix_store &m1, const local_matrix_store &m2,
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

void sapply(const local_matrix_store &store, const bulk_uoperate &op,
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

void apply(int margin, const arr_apply_operate &op,
		const local_matrix_store &in_mat, local_matrix_store &out_mat)
{
	assert(margin == apply_margin::MAR_ROW || margin == apply_margin::MAR_COL);
	// In these two cases, we need to convert the matrix store layout
	// before we can apply the function to the matrix.
	local_matrix_store::const_ptr buf_mat;
	if (in_mat.store_layout() == matrix_layout_t::L_COL
			&& margin == apply_margin::MAR_ROW) {
		buf_mat = in_mat.conv2(matrix_layout_t::L_ROW);
		assert(buf_mat);
	}
	else if (in_mat.store_layout() == matrix_layout_t::L_ROW
			&& margin == apply_margin::MAR_COL) {
		buf_mat = in_mat.conv2(matrix_layout_t::L_COL);
		assert(buf_mat);
	}

	const local_matrix_store *this_mat;
	if (buf_mat)
		this_mat = buf_mat.get();
	else
		this_mat = &in_mat;

	if (margin == apply_margin::MAR_ROW) {
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

void inner_prod(const local_matrix_store &m1, const local_matrix_store &m2,
		const bulk_operate &left_op, const bulk_operate &right_op,
		local_matrix_store &res)
{
	if (m1.store_layout() == matrix_layout_t::L_ROW) {
		local_matrix_store::ptr new_m2;
		const local_matrix_store *col_m2 = &m2;
		if (m2.store_layout() == matrix_layout_t::L_ROW) {
			new_m2 = m2.conv2(matrix_layout_t::L_COL);
			col_m2 = new_m2.get();
		}
		assert(col_m2->store_layout() == matrix_layout_t::L_COL);
		assert(res.store_layout() == matrix_layout_t::L_ROW);
		if (m1.is_wide())
			inner_prod_row_wide(static_cast<const local_row_matrix_store &>(m1),
					static_cast<const local_col_matrix_store &>(*col_m2), left_op,
					right_op, static_cast<local_row_matrix_store &>(res));
		else
			inner_prod_row_tall(static_cast<const local_row_matrix_store &>(m1),
					static_cast<const local_col_matrix_store &>(*col_m2), left_op,
					right_op, static_cast<local_row_matrix_store &>(res));
	}
	else {
		if (m1.is_wide()) {
			local_matrix_store::ptr new_m2;
			const local_matrix_store *col_m2 = &m2;
			// TODO we can handle the case of left wide col-major matrix and
			// right tall row-major matrix more efficiently, so we don't need
			// to convert the right matrix.
			if (m2.store_layout() == matrix_layout_t::L_ROW) {
				new_m2 = m2.conv2(matrix_layout_t::L_COL);
				col_m2 = new_m2.get();
			}
			assert(col_m2->store_layout() == matrix_layout_t::L_COL);
			assert(res.store_layout() == matrix_layout_t::L_ROW);
			inner_prod_col_wide(static_cast<const local_col_matrix_store &>(m1),
					static_cast<const local_col_matrix_store &>(*col_m2), left_op,
					right_op, static_cast<local_row_matrix_store &>(res));
		}
		else {
			assert(res.store_layout() == matrix_layout_t::L_COL);
			inner_prod_col_tall(static_cast<const local_col_matrix_store &>(m1),
					m2, left_op, right_op,
					static_cast<local_col_matrix_store &>(res));
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
};

double_multiply_operate dm_op;

}

void scale_cols(const local_matrix_store &store, const local_vec_store &vals,
		local_matrix_store &res)
{
	assert(res.store_layout() == store.store_layout());
	size_t ncol = store.get_num_cols();
	size_t nrow = store.get_num_rows();
	const bulk_operate *op = &store.get_type().get_basic_ops().get_multiply();
	if (store.get_type() == get_scalar_type<double>())
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

void scale_rows(const local_matrix_store &store, const local_vec_store &vals,
		local_matrix_store &res)
{
	assert(res.store_layout() == store.store_layout());
	size_t ncol = store.get_num_cols();
	size_t nrow = store.get_num_rows();
	const bulk_operate *op = &store.get_type().get_basic_ops().get_multiply();
	if (store.get_type() == get_scalar_type<double>())
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
