#ifndef __NUMA_DENSE_MATRIX_H__
#define __NUMA_DENSE_MATRIX_H__

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

#include "comm_exception.h"

#include "raw_data_array.h"
#include "NUMA_mapper.h"
#include "NUMA_vector.h"
#include "mem_matrix_store.h"

namespace fm
{

namespace detail
{

class NUMA_matrix_store: public mem_matrix_store
{
protected:
	NUMA_matrix_store(size_t nrow, size_t ncol,
			const scalar_type &type): mem_matrix_store(nrow, ncol, type) {
	}
public:
	typedef std::shared_ptr<NUMA_matrix_store> ptr;
	typedef std::shared_ptr<const NUMA_matrix_store> const_ptr;

	static ptr cast(matrix_store::ptr store);
	static const_ptr cast(matrix_store::const_ptr store);

	static ptr create(size_t nrow, size_t ncol, int num_nodes,
			matrix_layout_t layout, const scalar_type &type);

	virtual matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const {
		assert(0);
		return matrix_store::const_ptr();
	}
	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const {
		assert(0);
		return matrix_store::const_ptr();
	}
	virtual bool write2file(const std::string &file_name) const {
		assert(0);
		return false;
	}
};

class NUMA_row_matrix_store: public NUMA_matrix_store
{
protected:
	NUMA_row_matrix_store(size_t nrow, size_t ncol,
			const scalar_type &type): NUMA_matrix_store(nrow, ncol, type) {
	}
public:
	typedef std::shared_ptr<NUMA_row_matrix_store> ptr;

	static ptr cast(matrix_store::ptr store) {
		// TODO check store type.
		return std::static_pointer_cast<NUMA_row_matrix_store>(store);
	}

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_ROW;
	}
};

class NUMA_col_matrix_store: public NUMA_matrix_store
{
protected:
	NUMA_col_matrix_store(size_t nrow, size_t ncol,
			const scalar_type &type): NUMA_matrix_store(nrow, ncol, type) {
	}
public:
	typedef std::shared_ptr<NUMA_col_matrix_store> ptr;

	static ptr cast(matrix_store::ptr store) {
		// TODO check store type.
		return std::static_pointer_cast<NUMA_col_matrix_store>(store);
	}

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_COL;
	}
};

class NUMA_col_wide_matrix_store;
class NUMA_row_wide_matrix_store;

/*
 * This class defines an in-memory dense matrix whose data is organized
 * in row major and has many more rows than columns. The matrix is optimized
 * for a NUMA machine. Rows in the matrix are distributed across multiple
 * NUMA node. Multiple adjacent rows are stored in contiguous memory and
 * a single row is guaranteed to be stored together.
 */
class NUMA_row_tall_matrix_store: public NUMA_row_matrix_store
{
	// This is to map rows to different NUMA nodes.
	detail::NUMA_mapper mapper;
	std::vector<detail::raw_data_array> data;

	// The copy constructor performs shallow copy.
	NUMA_row_tall_matrix_store(
			const NUMA_row_tall_matrix_store &mat): NUMA_row_matrix_store(
			mat.get_num_rows(), mat.get_num_cols(),
			mat.get_type()), mapper(mat.get_num_nodes()) {
		this->data = mat.data;
	}

	NUMA_row_tall_matrix_store(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type);
public:
	typedef std::shared_ptr<NUMA_row_tall_matrix_store> ptr;

	static ptr create(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type) {
		return ptr(new NUMA_row_tall_matrix_store(nrow, ncol, num_nodes, type));
	}

	int get_num_nodes() const {
		return data.size();
	}

	const char *get_row(size_t row_idx) const;
	char *get_row(size_t row_idx);
	using NUMA_matrix_store::get_rows;
	const char *get_rows(size_t row_start, size_t row_end) const;
	char *get_rows(size_t row_start, size_t row_end);

	const char *get(size_t row_idx, size_t col_idx) const {
		const char *row = get_row(row_idx);
		return row + col_idx * get_entry_size();
	}

	char *get(size_t row_idx, size_t col_idx) {
		char *row = get_row(row_idx);
		return row + col_idx * get_entry_size();
	}


	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols);
	virtual std::shared_ptr<const local_matrix_store> get_portion(size_t id) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);

	virtual vec_store::const_ptr get_row_vec(off_t row) const {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't get a row from a NUMA tall row matrix";
		return vec_store::const_ptr();
	}
	virtual vec_store::const_ptr get_col_vec(off_t idx) const {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't get a column from a NUMA tall row matrix";
		return vec_store::const_ptr();
	}
	virtual matrix_store::const_ptr append_cols(
			const std::vector<matrix_store::const_ptr> &mats) const {
		throw unsupported_exception(
				"can't add columns to a row-major matrix");
	}

	virtual matrix_store::const_ptr transpose() const;
	virtual bool write2file(const std::string &file_name) const;

	friend class NUMA_col_wide_matrix_store;
};

/*
 * This class defines an in-memory dense matrix whose data is organized
 * in column major and has many more rows than columns. The matrix is optimized
 * for a NUMA machine. Each column in the matrix are distributed across multiple
 * NUMA nodes.
 */
class NUMA_col_tall_matrix_store: public NUMA_col_matrix_store
{
	std::vector<NUMA_vec_store::ptr> data;

	NUMA_col_tall_matrix_store(
			const std::vector<NUMA_vec_store::ptr> &cols): NUMA_col_matrix_store(
				cols.front()->get_length(), cols.size(),
				cols.front()->get_type()) {
		this->data = cols;
	}

	// The copy constructor performs shallow copy.
	NUMA_col_tall_matrix_store(
			const NUMA_col_tall_matrix_store &mat): NUMA_col_matrix_store(
			mat.get_num_rows(), mat.get_num_cols(), mat.get_type()) {
		this->data = mat.data;
	}

	NUMA_col_tall_matrix_store(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type);
public:
	typedef std::shared_ptr<NUMA_col_tall_matrix_store> ptr;

	static ptr create(const std::vector<NUMA_vec_store::ptr> &cols) {
		return ptr(new NUMA_col_tall_matrix_store(cols));
	}

	static ptr create(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type) {
		return ptr(new NUMA_col_tall_matrix_store(nrow, ncol, num_nodes, type));
	}

	int get_num_nodes() const {
		return data[0]->get_num_nodes();
	}

	char *get(size_t row_idx, size_t col_idx) {
		return data[col_idx]->get(row_idx);
	}

	const char *get(size_t row_idx, size_t col_idx) const {
		return data[col_idx]->get(row_idx);
	}


	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols);
	virtual std::shared_ptr<const local_matrix_store> get_portion(size_t id) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);

	virtual vec_store::ptr get_col_vec(off_t col) {
		return data[col];
	}
	virtual vec_store::const_ptr get_col_vec(off_t col) const {
		return data[col];
	}
	virtual vec_store::const_ptr get_row_vec(off_t idx) const {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't get a row from a NUMA tall column matrix";
		return vec_store::const_ptr();
	}
	virtual matrix_store::const_ptr append_cols(
			const std::vector<matrix_store::const_ptr> &mats) const;

	virtual matrix_store::const_ptr get_cols(const std::vector<off_t> &idxs) const;

	virtual matrix_store::const_ptr transpose() const;
	virtual bool write2file(const std::string &file_name) const;

	friend class NUMA_row_wide_matrix_store;
};

class NUMA_row_wide_matrix_store: public NUMA_row_matrix_store
{
	NUMA_col_tall_matrix_store store;

	NUMA_row_wide_matrix_store(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type): NUMA_row_matrix_store(nrow, ncol,
				type), store(ncol, nrow, num_nodes, type) {
	}

	/*
	 * The constructed matrix store will be the transpose of _store.
	 */
	NUMA_row_wide_matrix_store(
			const NUMA_col_tall_matrix_store &_store): NUMA_row_matrix_store(
				_store.get_num_cols(), _store.get_num_rows(),
				_store.get_type()), store(_store) {
	}
public:
	typedef std::shared_ptr<NUMA_row_wide_matrix_store> ptr;

	static ptr create(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type) {
		return ptr(new NUMA_row_wide_matrix_store(nrow, ncol, num_nodes, type));
	}

	static ptr create_transpose(const NUMA_col_tall_matrix_store &store) {
		return ptr(new NUMA_row_wide_matrix_store(store));
	}

	int get_num_nodes() const {
		return store.get_num_nodes();
	}

	char *get(size_t row_idx, size_t col_idx) {
		return store.get(col_idx, row_idx);
	}

	const char *get(size_t row_idx, size_t col_idx) const {
		return store.get(col_idx, row_idx);
	}

	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols);
	virtual std::shared_ptr<const local_matrix_store> get_portion(size_t id) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);

	virtual vec_store::const_ptr get_row_vec(off_t row) const {
		return store.get_col_vec(row);
	}
	virtual vec_store::const_ptr get_col_vec(off_t idx) const {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't get a column from a NUMA wide row matrix";
		return vec_store::const_ptr();
	}
	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const {
		return store.get_cols(idxs);
	}
	virtual matrix_store::const_ptr append_cols(
			const std::vector<matrix_store::const_ptr> &mats) const {
		throw unsupported_exception(
				"can't add columns to a row-major matrix");
	}

	virtual matrix_store::const_ptr transpose() const;
};

class NUMA_col_wide_matrix_store: public NUMA_col_matrix_store
{
	NUMA_row_tall_matrix_store store;

	NUMA_col_wide_matrix_store(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type): NUMA_col_matrix_store(nrow, ncol,
				type), store(ncol, nrow, num_nodes, type) {
	}

	/*
	 * The constructed matrix store will be the transpose of _store.
	 */
	NUMA_col_wide_matrix_store(
			const NUMA_row_tall_matrix_store &_store): NUMA_col_matrix_store(
				_store.get_num_cols(), _store.get_num_rows(),
				_store.get_type()), store(_store) {
	}
public:
	typedef std::shared_ptr<NUMA_col_wide_matrix_store> ptr;

	static ptr create(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type) {
		return ptr(new NUMA_col_wide_matrix_store(nrow, ncol, num_nodes, type));
	}

	static ptr create_transpose(const NUMA_row_tall_matrix_store &store) {
		return ptr(new NUMA_col_wide_matrix_store(store));
	}

	int get_num_nodes() const {
		return store.get_num_nodes();
	}

	const char *get_col(size_t col_idx) const {
		return store.get_row(col_idx);
	}

	char *get_col(size_t col_idx) {
		return store.get_row(col_idx);
	}

	using NUMA_matrix_store::get_cols;
	const char *get_cols(size_t col_start, size_t col_end) const {
		return store.get_rows(col_start, col_end);
	}

	char *get_cols(size_t col_start, size_t col_end) {
		return store.get_rows(col_start, col_end);
	}

	const char *get(size_t row_idx, size_t col_idx) const {
		return store.get(col_idx, row_idx);
	}

	char *get(size_t row_idx, size_t col_idx) {
		return store.get(col_idx, row_idx);
	}

	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols);
	virtual std::shared_ptr<const local_matrix_store> get_portion(size_t id) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);

	virtual vec_store::const_ptr get_col_vec(off_t col) const {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't get a column from a NUMA wide col matrix";
		return vec_store::const_ptr();
	}
	virtual vec_store::const_ptr get_row_vec(off_t row) const {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't get a row from a NUMA wide col matrix";
		return vec_store::const_ptr();
	}
	virtual matrix_store::const_ptr append_cols(
			const std::vector<matrix_store::const_ptr> &mats) const {
		throw unsupported_exception(
				"can't add columns to a wide column-major matrix");
	}

	virtual matrix_store::const_ptr transpose() const;
};

}

}

#endif
