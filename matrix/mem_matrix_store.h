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

#include <boost/format.hpp>

#include "comm_exception.h"

#include "matrix_store.h"
#include "raw_data_array.h"
#include "data_io.h"

namespace fm
{

namespace detail
{

class local_matrix_store;
class vec_store;

/*
 * This is the base class that represents an in-memory complete matrix.
 */
class mem_matrix_store: public matrix_store
{
	const size_t mat_id;
protected:
	bool write_header(FILE *f) const;
	size_t get_mat_id() const {
		return mat_id;
	}
public:
	typedef std::shared_ptr<mem_matrix_store> ptr;
	typedef std::shared_ptr<const mem_matrix_store> const_ptr;
	/*
	 * We partition a matrix for parallel.
	 */
	static const size_t CHUNK_SIZE;

	/*
	 * This is to load a dense matrix from a file that contains a matrix header.
	 */
	static const_ptr load(const std::string &file_name, int num_nodes);
	static const_ptr load(file_io::ptr io, int num_nodes);

	/*
	 * This is to load a dense matrix from a file that doesn't contain
	 * a matrix header.
	 */
	static const_ptr load_raw(const std::string &file_name, size_t nrow, size_t ncol,
			matrix_layout_t layout, const scalar_type &type, int num_nodes);
	static const_ptr load_raw(file_io::ptr io, size_t nrow, size_t ncol,
			matrix_layout_t layout, const scalar_type &type, int num_nodes);

	static ptr cast(matrix_store::ptr store);
	static const_ptr cast(matrix_store::const_ptr store);

	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, int num_nodes);

	mem_matrix_store(size_t nrow, size_t ncol, const scalar_type &type);

	virtual bool write2file(const std::string &file_name,
			bool text = false, std::string sep = " ") const;

	/*
	 * This function symmetrizes a matrix.
	 * It only works for a SMP square matrix.
	 */
	bool symmetrize(bool upper2lower);

	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		std::unordered_map<size_t, size_t> ret;
		// TODO right now we only indicate the matrix. We set the number of
		// bytes to 0
		// We should also use data_id instead of mat_id.
		ret.insert(std::pair<size_t, size_t>(mat_id, 0));
		return ret;
	}
	virtual std::string get_name() const;
	virtual bool share_data(const matrix_store &store) const;

	virtual const char *get(size_t row, size_t col) const = 0;
	virtual char *get(size_t row, size_t col) = 0;
	virtual const char *get_row(size_t row) const {
		return NULL;
	}
	virtual char *get_row(size_t row) {
		return NULL;
	}
	virtual const char *get_col(size_t col) const {
		return NULL;
	}
	virtual char *get_col(size_t col) {
		return NULL;
	}
	virtual const char *get_raw_arr() const {
		return NULL;
	}
	virtual char *get_raw_arr() {
		return NULL;
	}

	virtual async_cres_t get_portion_async(size_t start_row, size_t start_col,
			size_t num_rows, size_t num_cols,
			std::shared_ptr<portion_compute> compute) const {
		return async_cres_t(true,
				get_portion(start_row, start_col, num_rows, num_cols));
	}
	virtual void write_portion_async(
			std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col);

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
 * mem_row_matrix_store and mem_col_matrix_store are two simple implementations
 * of in-memory matrix stores. They store all data in a single piece of
 * contiguous memory. They are mainly used for storing small matrices.
 */

/*
 * This represents a column-major matrix. All columns are stored
 * in contiguous memory.
 */
class mem_col_matrix_store: public mem_matrix_store
{
	simple_raw_array data;

	mem_col_matrix_store(size_t nrow, size_t ncol,
			const scalar_type &type): mem_matrix_store(nrow, ncol, type) {
		if (nrow * ncol > 0)
			data = simple_raw_array(nrow * ncol * type.get_size(), -1);
	}
protected:
	mem_col_matrix_store(size_t nrow, size_t ncol, const scalar_type &type,
			const simple_raw_array &data): mem_matrix_store(nrow, ncol,
				type) {
		this->data = data;
	}
public:
	typedef std::shared_ptr<mem_col_matrix_store> ptr;
	typedef std::shared_ptr<const mem_col_matrix_store> const_ptr;

	static ptr create(const simple_raw_array &data, size_t nrow, size_t ncol,
			const scalar_type &type) {
		return ptr(new mem_col_matrix_store(nrow, ncol, type, data));
	}

	static ptr create(size_t nrow, size_t ncol, const scalar_type &type) {
		return ptr(new mem_col_matrix_store(nrow, ncol, type));
	}

	static const_ptr cast(matrix_store::const_ptr store);
	static ptr cast(matrix_store::ptr store);

	virtual std::shared_ptr<const vec_store> conv2vec() const;

	const simple_raw_array &get_data() const {
		return data;
	}
	virtual const char *get_raw_arr() const {
		return data.get_raw();
	}
	virtual char *get_raw_arr() {
		return data.get_raw();
	}

	virtual size_t get_data_id() const;

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
	virtual int get_portion_node_id(size_t id) const {
		return -1;
	}

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_store::const_ptr get_cols(const std::vector<off_t> &idxs) const;
	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const;

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_COL;
	}

	virtual bool resize(size_t num_rows, size_t num_cols);
};

/*
 * This represents a row-major matrix. All rows are stored in contiguous
 * memory.
 */
class mem_row_matrix_store: public mem_matrix_store
{
	simple_raw_array data;

	mem_row_matrix_store(size_t nrow, size_t ncol,
			const scalar_type &type): mem_matrix_store(nrow, ncol, type) {
		if (nrow * ncol > 0)
			data = simple_raw_array(nrow * ncol * type.get_size(), -1);
	}
protected:
	mem_row_matrix_store(size_t nrow, size_t ncol, const scalar_type &type,
			const simple_raw_array &data): mem_matrix_store(nrow, ncol,
				type) {
		this->data = data;
	}
public:
	typedef std::shared_ptr<mem_row_matrix_store> ptr;
	typedef std::shared_ptr<const mem_row_matrix_store> const_ptr;

	static ptr create(const simple_raw_array &data, size_t nrow, size_t ncol,
			const scalar_type &type) {
		return ptr(new mem_row_matrix_store(nrow, ncol, type, data));
	}

	static ptr create(size_t nrow, size_t ncol, const scalar_type &type) {
		return ptr(new mem_row_matrix_store(nrow, ncol, type));
	}

	static ptr cast(matrix_store::ptr store);
	static const_ptr cast(matrix_store::const_ptr store);

	virtual std::shared_ptr<const vec_store> conv2vec() const;

	const simple_raw_array &get_data() const {
		return data;
	}
	virtual const char *get_raw_arr() const {
		return data.get_raw();
	}
	virtual char *get_raw_arr() {
		return data.get_raw();
	}

	virtual size_t get_data_id() const;

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
	virtual int get_portion_node_id(size_t id) const {
		return -1;
	}

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const;

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_ROW;
	}

	virtual bool resize(size_t num_rows, size_t num_cols);
};

/*
 * This matrix contains a few columns of a column-major matrix.
 * The columns in this matrix isn't necessarily stored contiguously.
 */
class mem_sub_col_matrix_store: public mem_col_matrix_store
{
	std::shared_ptr<const std::vector<off_t> > orig_col_idxs;

	mem_sub_col_matrix_store(const mem_col_matrix_store &store,
			std::shared_ptr<const std::vector<off_t> > col_idxs): mem_col_matrix_store(
				store.get_num_rows(), col_idxs->size(), store.get_type(),
				store.get_data()) {
		this->orig_col_idxs = col_idxs;
	}
	mem_sub_col_matrix_store(const simple_raw_array &data,
			std::shared_ptr<const std::vector<off_t> > col_idxs, size_t nrow,
			const scalar_type &type): mem_col_matrix_store(nrow, col_idxs->size(),
				type, data) {
		this->orig_col_idxs = col_idxs;
	}
public:
	static ptr create(const simple_raw_array &data,
			std::shared_ptr<const std::vector<off_t> > col_idxs, size_t nrow,
			const scalar_type &type) {
		return ptr(new mem_sub_col_matrix_store(data, col_idxs, nrow, type));
	}

	/*
	 * The column indexes are the absolute index on the original column matrix.
	 */
	static ptr create(const mem_col_matrix_store &store,
			const std::vector<off_t> &abs_col_idxs) {
		std::shared_ptr<std::vector<off_t> > idxs(new std::vector<off_t>());
		*idxs = abs_col_idxs;
		return ptr(new mem_sub_col_matrix_store(store, idxs));
	}

	virtual size_t get_data_id() const;

	virtual char *get_col(size_t col) {
		return mem_col_matrix_store::get_col(orig_col_idxs->at(col));
	}

	virtual const char *get_col(size_t col) const {
		return mem_col_matrix_store::get_col(orig_col_idxs->at(col));
	}

	// Its data is likely not in contiguous memory.
	virtual const char *get_raw_arr() const {
		if (get_num_cols() == 1)
			return get_col(0);
		else
			return NULL;
	}
	virtual char *get_raw_arr() {
		if (get_num_cols() == 1)
			return get_col(0);
		else
			return NULL;
	}

	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols);

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_store::const_ptr get_cols(const std::vector<off_t> &idxs) const;

	virtual bool resize(size_t num_rows, size_t num_cols) {
		throw unsupported_exception("sub col matrix doesn't support resize");
	}
};

/*
 * This matrix contains a few rows of a row-major matrix.
 * The rows in this matrix isn't necessarily stored contiguously.
 */
class mem_sub_row_matrix_store: public mem_row_matrix_store
{
	std::shared_ptr<const std::vector<off_t> > orig_row_idxs;

	mem_sub_row_matrix_store(const mem_row_matrix_store &store,
			std::shared_ptr<const std::vector<off_t> > row_idxs): mem_row_matrix_store(
				row_idxs->size(), store.get_num_cols(), store.get_type(),
				store.get_data()) {
		this->orig_row_idxs = row_idxs;
	}
	mem_sub_row_matrix_store(const simple_raw_array &data,
			std::shared_ptr<const std::vector<off_t> > row_idxs, size_t ncol,
			const scalar_type &type): mem_row_matrix_store(row_idxs->size(),
				ncol, type, data) {
		this->orig_row_idxs = row_idxs;
	}
public:
	static ptr create(const simple_raw_array &data,
			std::shared_ptr<const std::vector<off_t> > row_idxs, size_t ncol,
			const scalar_type &type) {
		return ptr(new mem_sub_row_matrix_store(data, row_idxs, ncol, type));
	}
	/*
	 * The row indexes are the absolute index on the original row matrix.
	 */
	static ptr create(const mem_row_matrix_store &store,
			const std::vector<off_t> &abs_row_idxs) {
		std::shared_ptr<std::vector<off_t> > idxs(new std::vector<off_t>());
		*idxs = abs_row_idxs;
		return ptr(new mem_sub_row_matrix_store(store, idxs));
	}

	virtual size_t get_data_id() const;

	virtual char *get_row(size_t row) {
		return mem_row_matrix_store::get_row(orig_row_idxs->at(row));
	}

	virtual const char *get_row(size_t row) const {
		return mem_row_matrix_store::get_row(orig_row_idxs->at(row));
	}

	// Its data is likely not in contiguous memory.
	virtual const char *get_raw_arr() const {
		if (get_num_rows() == 1)
			return get_row(0);
		else
			return NULL;
	}
	virtual char *get_raw_arr() {
		if (get_num_rows() == 1)
			return get_row(0);
		else
			return NULL;
	}

	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols);

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const;

	virtual bool resize(size_t num_rows, size_t num_cols) {
		throw unsupported_exception("sub row matrix doesn't support resize");
	}
};

class mem_matrix_stream: public matrix_stream
{
	mem_matrix_store::ptr store;

	mem_matrix_stream(mem_matrix_store::ptr store) {
		this->store = store;
	}
public:
	static ptr create(mem_matrix_store::ptr store) {
		return ptr(new mem_matrix_stream(store));
	}

	virtual void write_async(std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col) {
		store->write_portion_async(portion, start_row, start_col);
	}
	virtual bool is_complete() const {
		// We write data to in-mem matrix directly, so it's always complete.
		return true;
	}
	virtual const matrix_store &get_mat() const {
		return *store;
	}
};

}

}

#endif
