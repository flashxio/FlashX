#ifndef __LOCAL_MATRIX_STORE_H__
#define __LOCAL_MATRIX_STORE_H__

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

#include <memory>

#include "matrix_store.h"
#include "raw_data_array.h"

namespace fm
{

class bulk_operate;
class bulk_uoperate;

namespace detail
{

/*
 * This is a base class that represents part of a matrix.
 * The original matrix can be an SMP matrix, a NUMA matrix,
 * an external-memory matrix or a distributed-memory matrix.
 * It is guaranteed that this matrix stores the portion of the original
 * matrix in the local memory.
 * This class has very similar interface as mem_matrix_store, but it's
 * used for different purpose.
 */
class local_matrix_store: public matrix_store
{
	off_t global_start_row;
	off_t global_start_col;
public:
	typedef std::shared_ptr<local_matrix_store> ptr;
	typedef std::shared_ptr<const local_matrix_store> const_ptr;

	local_matrix_store(off_t global_start_row, off_t global_start_col,
			size_t nrow, size_t ncol, const scalar_type &type): matrix_store(
				nrow, ncol, true, type) {
		this->global_start_row = global_start_row;
		this->global_start_col = global_start_col;
	}

	off_t get_global_start_row() const {
		return global_start_row;
	}
	off_t get_global_start_col() const {
		return global_start_col;
	}

	virtual bool read_only() const = 0;
	virtual const char *get_raw_arr() const = 0;
	virtual char *get_raw_arr() = 0;
	virtual const char *get(size_t row, size_t col) const = 0;
	virtual char *get(size_t row, size_t col) = 0;

	virtual matrix_store::const_ptr transpose() const {
		assert(0);
		return matrix_store::const_ptr();
	}

	template<class Type>
	Type get(size_t row, size_t col) const {
		return *(const Type *) get(row, col);
	}

	template<class Type>
	void set(size_t row, size_t col, Type val) {
		*(Type *) get(row, col) = val;
	}
};

class local_col_matrix_store: public local_matrix_store
{
public:
	local_col_matrix_store(off_t global_start_row, off_t global_start_col,
			size_t nrow, size_t ncol, const scalar_type &type): local_matrix_store(
				global_start_row, global_start_col, nrow, ncol, type) {
	}

	virtual void reset_data();
	virtual void set_data(const set_operate &op);

	virtual const char *get_col(size_t col) const = 0;
	virtual char *get_col(size_t col) = 0;

	virtual const char *get(size_t row, size_t col) const {
		return get_col(col) + row * get_entry_size();
	}

	virtual char *get(size_t row, size_t col) {
		return get_col(col) + row * get_entry_size();
	}

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_COL;
	}
};

class local_row_matrix_store: public local_matrix_store
{
public:
	local_row_matrix_store(off_t global_start_row, off_t global_start_col,
			size_t nrow, size_t ncol, const scalar_type &type): local_matrix_store(
				global_start_row, global_start_col, nrow, ncol, type) {
	}

	virtual void reset_data();
	virtual void set_data(const set_operate &op);

	virtual const char *get_row(size_t row) const = 0;
	virtual char *get_row(size_t row) = 0;

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_ROW;
	}

	virtual const char *get(size_t row, size_t col) const {
		return get_row(row) + col * get_entry_size();
	}
	virtual char *get(size_t row, size_t col) {
		return get_row(row) + col * get_entry_size();
	}
};

/*
 * A matrix owns data to store portion of data in a column-major matrix.
 */
class local_buf_col_matrix_store: public local_col_matrix_store
{
	raw_data_array data;
public:
	local_buf_col_matrix_store(off_t global_start_row, off_t global_start_col,
			size_t nrow, size_t ncol, const scalar_type &type): local_col_matrix_store(
				global_start_row, global_start_col, nrow, ncol, type) {
		if (nrow * ncol > 0)
			data = raw_data_array(nrow * ncol * type.get_size());
	}

	virtual bool read_only() const {
		return false;
	}

	virtual const char *get_raw_arr() const {
		return data.get_raw();
	}

	virtual char *get_raw_arr() {
		return data.get_raw();
	}

	virtual const char *get_col(size_t col) const {
		return data.get_raw() + col * get_num_rows() * get_entry_size();
	}
	virtual char *get_col(size_t col) {
		return data.get_raw() + col * get_num_rows() * get_entry_size();
	}
};

/*
 * A matrix owns data to store portion of data in a row-major matrix.
 */
class local_buf_row_matrix_store: public local_row_matrix_store
{
	raw_data_array data;
public:
	local_buf_row_matrix_store(off_t global_start_row, off_t global_start_col,
			size_t nrow, size_t ncol, const scalar_type &type): local_row_matrix_store(
				global_start_row, global_start_col, nrow, ncol, type) {
		if (nrow * ncol > 0)
			data = raw_data_array(nrow * ncol * type.get_size());
	}

	virtual bool read_only() const {
		return false;
	}

	virtual const char *get_raw_arr() const {
		return data.get_raw();
	}

	virtual char *get_raw_arr() {
		return data.get_raw();
	}

	virtual const char *get_row(size_t row) const {
		return data.get_raw() + row * get_num_cols() * get_entry_size();
	}

	virtual char *get_row(size_t row) {
		return data.get_raw() + row * get_num_cols() * get_entry_size();
	}
};

/*
 * A matrix that references portion of data in another column-major matrix.
 * The referenced data is stored contiguously.
 */
class local_ref_contig_col_matrix_store: public local_col_matrix_store
{
	char *data;
public:
	local_ref_contig_col_matrix_store(off_t global_start_row, off_t global_start_col,
			char *data, size_t nrow, size_t ncol,
			const scalar_type &type): local_col_matrix_store(global_start_row,
				global_start_col, nrow, ncol, type) {
		this->data = data;
	}

	virtual bool read_only() const {
		return false;
	}

	virtual const char *get_raw_arr() const {
		return data;
	}

	virtual char *get_raw_arr() {
		return data;
	}

	virtual const char *get_col(size_t col) const {
		return data + col * get_num_rows() * get_entry_size();
	}
	virtual char *get_col(size_t col) {
		return data + col * get_num_rows() * get_entry_size();
	}
};

/*
 * A matrix that references portion of data in another row-major matrix.
 * The referenced data is stored contiguously.
 */
class local_ref_contig_row_matrix_store: public local_row_matrix_store
{
	char *data;
public:
	local_ref_contig_row_matrix_store(off_t global_start_row, off_t global_start_col,
			char *data, size_t nrow, size_t ncol,
			const scalar_type &type): local_row_matrix_store(global_start_row,
				global_start_col, nrow, ncol, type) {
		this->data = data;
	}

	virtual bool read_only() const {
		return false;
	}

	virtual const char *get_raw_arr() const {
		return data;
	}

	virtual char *get_raw_arr() {
		return data;
	}

	virtual const char *get_row(size_t row) const {
		return data + row * get_num_cols() * get_entry_size();
	}

	virtual char *get_row(size_t row) {
		return data + row * get_num_cols() * get_entry_size();
	}
};

/*
 * A matrix that references portion of data in another column-major matrix.
 * The referenced data isn't guaranteed to be stored contiguously.
 */
class local_ref_col_matrix_store: public local_col_matrix_store
{
	std::vector<char *> cols;
public:
	local_ref_col_matrix_store(off_t global_start_row, off_t global_start_col,
			const std::vector<char *> &cols, size_t nrow,
			const scalar_type &type): local_col_matrix_store(global_start_row,
				global_start_col, nrow, cols.size(), type) {
		this->cols = cols;
	}

	virtual bool read_only() const {
		return false;
	}

	virtual const char *get_raw_arr() const {
		return NULL;
	}

	virtual char *get_raw_arr() {
		return NULL;
	}

	virtual const char *get_col(size_t col) const {
		return cols[col];
	}
	virtual char *get_col(size_t col) {
		return cols[col];
	}
};

/*
 * A matrix that references portion of data in another row-major matrix.
 * The referenced data isn't guaranteed to be stored contiguously.
 */
class local_ref_row_matrix_store: public local_row_matrix_store
{
	std::vector<char *> rows;
public:
	local_ref_row_matrix_store(off_t global_start_row, off_t global_start_col,
			const std::vector<char *> &rows, size_t ncol,
			const scalar_type &type): local_row_matrix_store(global_start_row,
				global_start_col, rows.size(), ncol, type) {
		this->rows = rows;
	}

	virtual bool read_only() const {
		return false;
	}

	virtual const char *get_raw_arr() const {
		return NULL;
	}

	virtual char *get_raw_arr() {
		return NULL;
	}

	virtual const char *get_row(size_t row) const {
		return rows[row];
	}

	virtual char *get_row(size_t row) {
		return rows[row];
	}
};

/*
 * A matrix that references portion of const data in another column-major matrix.
 * The referenced data is stored contiguously.
 */
class local_cref_contig_col_matrix_store: public local_col_matrix_store
{
	const char *data;
public:
	local_cref_contig_col_matrix_store(off_t global_start_row, off_t global_start_col,
			const char *data, size_t nrow, size_t ncol,
			const scalar_type &type): local_col_matrix_store(global_start_row,
				global_start_col, nrow, ncol, type) {
		this->data = data;
	}

	virtual bool read_only() const {
		return true;
	}

	virtual const char *get_raw_arr() const {
		return data;
	}
	virtual const char *get_col(size_t col) const {
		return data + col * get_num_rows() * get_entry_size();
	}

	virtual char *get_raw_arr() {
		assert(0);
		return NULL;
	}
	virtual char *get_col(size_t col) {
		assert(0);
		return NULL;
	}
};

/*
 * A matrix that references portion of const data in another row-major matrix.
 * The referenced data is stored contiguously.
 */
class local_cref_contig_row_matrix_store: public local_row_matrix_store
{
	const char *data;
public:
	local_cref_contig_row_matrix_store(off_t global_start_row, off_t global_start_col,
			const char *data, size_t nrow, size_t ncol,
			const scalar_type &type): local_row_matrix_store(global_start_row,
				global_start_col, nrow, ncol, type) {
		this->data = data;
	}

	virtual bool read_only() const {
		return true;
	}

	virtual const char *get_raw_arr() const {
		return data;
	}
	virtual const char *get_row(size_t row) const {
		return data + row * get_num_cols() * get_entry_size();
	}

	virtual char *get_raw_arr() {
		assert(0);
		return NULL;
	}
	virtual char *get_row(size_t row) {
		assert(0);
		return NULL;
	}
};

/*
 * A matrix that references portion of data in another column-major matrix.
 * The referenced data isn't guaranteed to be stored contiguously.
 */
class local_cref_col_matrix_store: public local_col_matrix_store
{
	std::vector<const char *> cols;
public:
	local_cref_col_matrix_store(off_t global_start_row, off_t global_start_col,
			const std::vector<const char *> &cols, size_t nrow,
			const scalar_type &type): local_col_matrix_store(global_start_row,
				global_start_col, nrow, cols.size(), type) {
		this->cols = cols;
	}

	virtual bool read_only() const {
		return true;
	}

	virtual const char *get_col(size_t col) const {
		return cols[col];
	}

	virtual const char *get_raw_arr() const {
		return NULL;
	}
	virtual char *get_raw_arr() {
		assert(0);
		return NULL;
	}
	virtual char *get_col(size_t col) {
		assert(0);
		return NULL;
	}
};

/*
 * A matrix that references portion of data in another row-major matrix.
 * The referenced data isn't guaranteed to be stored contiguously.
 */
class local_cref_row_matrix_store: public local_row_matrix_store
{
	std::vector<const char *> rows;
public:
	local_cref_row_matrix_store(off_t global_start_row, off_t global_start_col,
			const std::vector<const char *> &rows, size_t ncol,
			const scalar_type &type): local_row_matrix_store(global_start_row,
				global_start_col, rows.size(), ncol, type) {
		this->rows = rows;
	}

	virtual bool read_only() const {
		return true;
	}

	virtual const char *get_row(size_t row) const {
		return rows[row];
	}

	virtual const char *get_raw_arr() const {
		return NULL;
	}
	virtual char *get_raw_arr() {
		assert(0);
		return NULL;
	}
	virtual char *get_row(size_t row) {
		assert(0);
		return NULL;
	}
};

/*
 * These are the general operations on the local matrix store.
 */
void aggregate(const local_matrix_store &store, const bulk_operate &op, char *res);
void mapply2(const local_matrix_store &m1, const local_matrix_store &m2,
			const bulk_operate &op, local_matrix_store &res);
void sapply(const local_matrix_store &store, const bulk_uoperate &op,
		local_matrix_store &res);
void inner_prod(const local_matrix_store &m1, const local_matrix_store &m2,
		const bulk_operate &left_op, const bulk_operate &right_op,
		local_matrix_store &res);

}

}

#endif
