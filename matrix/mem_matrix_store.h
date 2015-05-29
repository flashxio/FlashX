#ifndef __MEM_MATRIX_STORE_H__
#define __MEM_MATRIX_STORE_H__

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
#include "raw_data_array.h"

namespace fm
{

namespace detail
{

class local_matrix_store;
class vec_store;

/*
 * This is the base class that represents an in-memory complete matrix.
 * It is used for SMP.
 */
class mem_matrix_store: public matrix_store
{
protected:
	/*
	 * We partition a matrix for parallel.
	 */
	static const size_t CHUNK_SIZE;

	bool write_header(FILE *f) const;
public:
	typedef std::shared_ptr<mem_matrix_store> ptr;
	typedef std::shared_ptr<const mem_matrix_store> const_ptr;

	static ptr load(const std::string &file_name);
	static ptr cast(matrix_store::ptr store);
	static const_ptr cast(matrix_store::const_ptr store);

	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, int num_nodes);

	mem_matrix_store(size_t nrow, size_t ncol,
			const scalar_type &type): matrix_store(nrow, ncol, true, type) {
	}

	virtual void reset_data();
	virtual void set_data(const set_operate &op);
	virtual matrix_store::ptr conv2(matrix_layout_t layout) const;

	virtual int get_num_nodes() const {
		return -1;
	}

	virtual const char *get(size_t row, size_t col) const = 0;
	virtual char *get(size_t row, size_t col) = 0;

	virtual matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const = 0;
	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const = 0;
	virtual std::shared_ptr<const vec_store> get_col_vec(off_t idx) const = 0;
	virtual std::shared_ptr<const vec_store> get_row_vec(off_t idx) const = 0;

	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const {
		assert(0);
		return std::shared_ptr<const local_matrix_store>();
	}
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) {
		assert(0);
		return std::shared_ptr<local_matrix_store>();
	}
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t id) const = 0;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t id) = 0;

	virtual bool write2file(const std::string &file_name) const = 0;

	virtual std::pair<size_t, size_t> get_portion_size() const {
		if (is_wide())
			return std::pair<size_t, size_t>(get_num_rows(), CHUNK_SIZE);
		else
			return std::pair<size_t, size_t>(CHUNK_SIZE, get_num_cols());
	}

	template<class T>
	T get(size_t row, size_t col) const {
		return *(const T *) get(row, col);
	}
	template<class T>
	void set(size_t row, size_t col, T val) {
		*(T *) get(row, col) = val;
	}
};

/*
 * This represents a column-major matrix. All columns are stored
 * in contiguous memory.
 */
class mem_col_matrix_store: public mem_matrix_store
{
	raw_data_array data;

	mem_col_matrix_store(size_t nrow, size_t ncol,
			const scalar_type &type): mem_matrix_store(nrow, ncol, type) {
		if (nrow * ncol > 0)
			data = raw_data_array(nrow * ncol * type.get_size());
	}
protected:
	mem_col_matrix_store(size_t nrow, size_t ncol, const scalar_type &type,
			const raw_data_array &data): mem_matrix_store(nrow, ncol,
				type) {
		this->data = data;
	}
public:
	typedef std::shared_ptr<mem_col_matrix_store> ptr;
	typedef std::shared_ptr<const mem_col_matrix_store> const_ptr;

	static ptr create(const raw_data_array &data, size_t nrow, size_t ncol,
			const scalar_type &type) {
		return ptr(new mem_col_matrix_store(nrow, ncol, type, data));
	}

	static ptr create(size_t nrow, size_t ncol, const scalar_type &type) {
		return ptr(new mem_col_matrix_store(nrow, ncol, type));
	}

	static const_ptr cast(matrix_store::const_ptr store);
	static ptr cast(matrix_store::ptr store);

	const raw_data_array &get_data() const {
		return data;
	}

	virtual const char *get_col(size_t col) const {
		return data.get_raw() + col * get_num_rows() * get_entry_size();
	}
	virtual char *get_col(size_t col) {
		return data.get_raw() + col * get_num_rows() * get_entry_size();
	}

	virtual const char *get(size_t row, size_t col) const {
		return get_col(col) + row * get_entry_size();
	}
	virtual char *get(size_t row, size_t col) {
		return get_col(col) + row * get_entry_size();
	}

	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols);
	virtual std::shared_ptr<const local_matrix_store> get_portion(size_t id) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_store::const_ptr get_cols(const std::vector<off_t> &idxs) const;
	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const {
		return matrix_store::const_ptr();
	}

	virtual std::shared_ptr<const vec_store> get_col_vec(off_t idx) const;
	virtual std::shared_ptr<const vec_store> get_row_vec(off_t idx) const;

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_COL;
	}
	virtual bool write2file(const std::string &file_name) const;
};

/*
 * This represents a row-major matrix. All rows are stored in contiguous
 * memory.
 */
class mem_row_matrix_store: public mem_matrix_store
{
	raw_data_array data;

	mem_row_matrix_store(size_t nrow, size_t ncol,
			const scalar_type &type): mem_matrix_store(nrow, ncol, type) {
		if (nrow * ncol > 0)
			data = raw_data_array(nrow * ncol * type.get_size());
	}
protected:
	mem_row_matrix_store(size_t nrow, size_t ncol, const scalar_type &type,
			const raw_data_array &data): mem_matrix_store(nrow, ncol,
				type) {
		this->data = data;
	}
public:
	typedef std::shared_ptr<mem_row_matrix_store> ptr;
	typedef std::shared_ptr<const mem_row_matrix_store> const_ptr;

	static ptr create(const raw_data_array &data, size_t nrow, size_t ncol,
			const scalar_type &type) {
		return ptr(new mem_row_matrix_store(nrow, ncol, type, data));
	}

	static ptr create(size_t nrow, size_t ncol, const scalar_type &type) {
		return ptr(new mem_row_matrix_store(nrow, ncol, type));
	}

	static ptr cast(matrix_store::ptr store);
	static const_ptr cast(matrix_store::const_ptr store);

	const raw_data_array &get_data() const {
		return data;
	}

	virtual const char *get_row(size_t row) const {
		return data.get_raw() + row * get_num_cols() * get_entry_size();
	}

	virtual char *get_row(size_t row) {
		return data.get_raw() + row * get_num_cols() * get_entry_size();
	}

	virtual const char *get(size_t row, size_t col) const {
		return get_row(row) + col * get_entry_size();
	}
	virtual char *get(size_t row, size_t col) {
		return get_row(row) + col * get_entry_size();
	}

	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols);
	virtual std::shared_ptr<const local_matrix_store> get_portion(size_t id) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const {
		return matrix_store::const_ptr();
	}
	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const;
	virtual std::shared_ptr<const vec_store> get_col_vec(off_t idx) const;
	virtual std::shared_ptr<const vec_store> get_row_vec(off_t idx) const;

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_ROW;
	}
	virtual bool write2file(const std::string &file_name) const;
};

/*
 * This matrix contains a few columns of a column-major matrix.
 * The columns in this matrix isn't necessarily stored contiguously.
 */
class mem_sub_col_matrix_store: public mem_col_matrix_store
{
	std::vector<off_t> orig_col_idxs;

	mem_sub_col_matrix_store(const mem_col_matrix_store &store,
			const std::vector<off_t> &col_idxs): mem_col_matrix_store(
				store.get_num_rows(), col_idxs.size(), store.get_type(),
				store.get_data()) {
		this->orig_col_idxs = col_idxs;
	}
	mem_sub_col_matrix_store(const raw_data_array &data,
			const std::vector<off_t> &col_idxs, size_t nrow,
			const scalar_type &type): mem_col_matrix_store(nrow, col_idxs.size(),
				type, data) {
		this->orig_col_idxs = col_idxs;
	}
public:
	static ptr create(const raw_data_array &data,
			const std::vector<off_t> &col_idxs, size_t nrow,
			const scalar_type &type) {
		return ptr(new mem_sub_col_matrix_store(data, col_idxs, nrow, type));
	}

	/*
	 * The column indexes are the absolute index on the original column matrix.
	 */
	static ptr create(const mem_col_matrix_store &store,
			const std::vector<off_t> &abs_col_idxs) {
		return ptr(new mem_sub_col_matrix_store(store, abs_col_idxs));
	}

	virtual char *get_col(size_t col) {
		return mem_col_matrix_store::get_col(orig_col_idxs[col]);
	}

	virtual const char *get_col(size_t col) const {
		return mem_col_matrix_store::get_col(orig_col_idxs[col]);
	}

	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const {
		// TODO
		assert(0);
		return std::shared_ptr<const local_matrix_store>();
	}
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) {
		// TODO
		assert(0);
		return std::shared_ptr<local_matrix_store>();
	}
	virtual std::shared_ptr<const local_matrix_store> get_portion(size_t id) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_store::const_ptr get_cols(const std::vector<off_t> &idxs) const;
	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const {
		return matrix_store::const_ptr();
	}

	virtual std::shared_ptr<const vec_store> get_col_vec(off_t idx) const;
	virtual std::shared_ptr<const vec_store> get_row_vec(off_t idx) const;
};

/*
 * This matrix contains a few rows of a row-major matrix.
 * The rows in this matrix isn't necessarily stored contiguously.
 */
class mem_sub_row_matrix_store: public mem_row_matrix_store
{
	std::vector<off_t> orig_row_idxs;

	mem_sub_row_matrix_store(const mem_row_matrix_store &store,
			const std::vector<off_t> &row_idxs): mem_row_matrix_store(
				row_idxs.size(), store.get_num_cols(), store.get_type(),
				store.get_data()) {
		this->orig_row_idxs = row_idxs;
	}
	mem_sub_row_matrix_store(const raw_data_array &data,
			const std::vector<off_t> &row_idxs, size_t ncol,
			const scalar_type &type): mem_row_matrix_store(row_idxs.size(),
				ncol, type, data) {
		this->orig_row_idxs = row_idxs;
	}
public:
	static ptr create(const raw_data_array &data,
			const std::vector<off_t> &row_idxs, size_t ncol,
			const scalar_type &type) {
		return ptr(new mem_sub_row_matrix_store(data, row_idxs, ncol, type));
	}
	/*
	 * The row indexes are the absolute index on the original row matrix.
	 */
	static ptr create(const mem_row_matrix_store &store,
			const std::vector<off_t> &abs_row_idxs) {
		return ptr(new mem_sub_row_matrix_store(store, abs_row_idxs));
	}

	virtual char *get_row(size_t row) {
		return mem_row_matrix_store::get_row(orig_row_idxs[row]);
	}

	virtual const char *get_row(size_t row) const {
		return mem_row_matrix_store::get_row(orig_row_idxs[row]);
	}

	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const {
		// TODO
		assert(0);
		return std::shared_ptr<const local_matrix_store>();
	}
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) {
		// TODO
		assert(0);
		return std::shared_ptr<local_matrix_store>();
	}
	virtual std::shared_ptr<const local_matrix_store> get_portion(size_t id) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const {
		return matrix_store::const_ptr();
	}
	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const;
	virtual std::shared_ptr<const vec_store> get_col_vec(off_t idx) const;
	virtual std::shared_ptr<const vec_store> get_row_vec(off_t idx) const;
};

}

}

#endif
