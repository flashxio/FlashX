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

/*
 * NUMA_matrix_store is optimized for storing large matrices. It mainly targets
 * large NUMA machines, but it should also be used to for large matrices in
 * a SMP machine. It has the same format as the external-memory dense matrix.
 * It chunks a dense matrix into portions horizontally for tall matrices and
 * vertically for wide matrices and stores all elements in a portion
 * contiguously. In this way, we can allocate memory in chunks to store
 * a matrix and memory chunks can be reused when the matrix is destroyed.
 * The reason that we want to reuse memory chunks is that it is expensive
 * to allocate large memory chunks. When we allocate large memory chunks
 * and populate them with pages, we keep them and reuse them to avoid
 * the overhead of populating memory with pages.
 */

class NUMA_matrix_store: public mem_matrix_store
{
protected:
	const data_id_t::ptr data_id;
	NUMA_matrix_store(size_t nrow, size_t ncol, const scalar_type &type,
			data_id_t::ptr id = data_id_t::create(mat_counter++)): mem_matrix_store(
			nrow, ncol, type), data_id(id) {
	}
public:
	typedef std::shared_ptr<NUMA_matrix_store> ptr;
	typedef std::shared_ptr<const NUMA_matrix_store> const_ptr;

	static ptr cast(matrix_store::ptr store);
	static const_ptr cast(matrix_store::const_ptr store);

	static ptr create(size_t nrow, size_t ncol, int num_nodes,
			matrix_layout_t layout, const scalar_type &type);

	virtual void inc_dag_ref(size_t id) {
		data_id->inc_ref(id);
	}
	virtual void reset_dag_ref() {
		data_id->reset_ref();
	}
	virtual size_t get_dag_ref() const {
		return data_id->get_ref();
	}

	virtual size_t get_data_id() const {
		return data_id->get_id();
	}
	virtual bool share_data(const matrix_store &store) const {
		return matrix_store::share_data(store);
	}

	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		std::unordered_map<size_t, size_t> ret;
		ret.insert(std::pair<size_t, size_t>(data_id->get_id(),
					get_num_rows() * get_num_cols()));
		return ret;
	}
	virtual void write_portion_async(
			std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col);
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
class NUMA_row_tall_matrix_store: public NUMA_matrix_store
{
	// This is to map rows to different NUMA nodes.
	NUMA_mapper mapper;
	std::vector<detail::chunked_raw_array> data;

	// The copy constructor performs shallow copy.
	NUMA_row_tall_matrix_store(const NUMA_row_tall_matrix_store &mat);

	NUMA_row_tall_matrix_store(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type);
	NUMA_row_tall_matrix_store(const std::vector<detail::chunked_raw_array> &data,
			size_t nrow, size_t ncol, const NUMA_mapper &mapper,
			const scalar_type &type);

	bool get_portion_check(size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;

public:
	typedef std::shared_ptr<NUMA_row_tall_matrix_store> ptr;

	static ptr create(const std::vector<detail::chunked_raw_array> &data,
			size_t nrow, size_t ncol, const NUMA_mapper &mapper,
			const scalar_type &type) {
		return ptr(new NUMA_row_tall_matrix_store(data, nrow, ncol,
					mapper, type));
	}

	static ptr create(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type) {
		return ptr(new NUMA_row_tall_matrix_store(nrow, ncol, num_nodes, type));
	}

	virtual vec_store::const_ptr conv2vec() const {
		return NUMA_vec_store::create(get_num_rows() * get_num_cols(), get_type(),
				data, mapper);
	}

	int get_num_nodes() const {
		return data.size();
	}

	const char *get_row(size_t row_idx) const;
	char *get_row(size_t row_idx);
	using NUMA_matrix_store::get_rows;

	const char *get(size_t row_idx, size_t col_idx) const {
		const char *row = get_row(row_idx);
		return row + col_idx * get_entry_size();
	}

	char *get(size_t row_idx, size_t col_idx) {
		char *row = get_row(row_idx);
		return row + col_idx * get_entry_size();
	}


	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const;
	virtual matrix_store::const_ptr get_cols(const std::vector<off_t> &idxs) const;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols);
	virtual std::shared_ptr<const local_matrix_store> get_portion(size_t id) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);
	virtual int get_portion_node_id(size_t id) const;

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_ROW;
	}

	virtual bool resize(size_t num_rows, size_t num_cols);

	friend class NUMA_col_wide_matrix_store;
};

/*
 * This class defines an in-memory dense matrix whose data is organized
 * in column major and has many more rows than columns. The matrix is optimized
 * for a NUMA machine. Each column in the matrix are distributed across multiple
 * NUMA nodes.
 */
class NUMA_col_tall_matrix_store: public NUMA_matrix_store
{
	// This is to map rows to different NUMA nodes.
	NUMA_mapper mapper;
	std::vector<detail::chunked_raw_array> data;

	// The copy constructor performs shallow copy.
	NUMA_col_tall_matrix_store(
			const NUMA_col_tall_matrix_store &mat): NUMA_matrix_store(
			mat.get_num_rows(), mat.get_num_cols(), mat.get_type(),
			mat.data_id), mapper(mat.mapper) {
		this->data = mat.data;
	}

	NUMA_col_tall_matrix_store(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type);
	NUMA_col_tall_matrix_store(
			const std::vector<detail::chunked_raw_array>& data,
			size_t nrow, size_t ncol, const NUMA_mapper &_mapper,
			const scalar_type &type): NUMA_matrix_store(nrow, ncol, type),
			mapper(_mapper) {
		this->data = data;
	}
public:
	typedef std::shared_ptr<NUMA_col_tall_matrix_store> ptr;

	static ptr create(const std::vector<detail::chunked_raw_array>& data,
			size_t nrow, size_t ncol, const NUMA_mapper &mapper,
			const scalar_type &type) {
		return ptr(new NUMA_col_tall_matrix_store(data, nrow, ncol,
					mapper, type));
	}

	static ptr create(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type) {
		return ptr(new NUMA_col_tall_matrix_store(nrow, ncol, num_nodes, type));
	}

	virtual vec_store::const_ptr conv2vec() const {
		return NUMA_vec_store::create(get_num_rows() * get_num_cols(), get_type(),
				data, mapper);
	}

	int get_num_nodes() const {
		return data.size();
	}

	char *get(size_t row_idx, size_t col_idx);
	const char *get(size_t row_idx, size_t col_idx) const;

	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const;
	virtual matrix_store::const_ptr get_cols(const std::vector<off_t> &idxs) const;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols);
	virtual std::shared_ptr<const local_matrix_store> get_portion(size_t id) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);
	virtual int get_portion_node_id(size_t id) const;

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_COL;
	}

	virtual bool resize(size_t num_rows, size_t num_cols);

	friend class NUMA_row_wide_matrix_store;
};

class NUMA_row_wide_matrix_store: public NUMA_matrix_store
{
	NUMA_col_tall_matrix_store store;

	NUMA_row_wide_matrix_store(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type): NUMA_matrix_store(nrow, ncol, type),
			store(ncol, nrow, num_nodes, type) {
	}

	/*
	 * The constructed matrix store will be the transpose of _store.
	 */
	NUMA_row_wide_matrix_store(
			const NUMA_col_tall_matrix_store &_store): NUMA_matrix_store(
				_store.get_num_cols(), _store.get_num_rows(),
				_store.get_type(), _store.data_id), store(_store) {
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

	virtual vec_store::const_ptr conv2vec() const {
		return store.conv2vec();
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
	virtual int get_portion_node_id(size_t id) const {
		return store.get_portion_node_id(id);
	}

	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const {
		auto cols = store.get_cols(idxs);
		if (cols == NULL)
			return matrix_store::const_ptr();
		else
			return cols->transpose();
	}

	virtual matrix_store::const_ptr transpose() const;

	virtual bool resize(size_t num_rows, size_t num_cols) {
		bool ret = store.resize(num_cols, num_rows);
		if (!ret)
			return false;
		return matrix_store::resize(num_rows, num_cols);
	}

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_ROW;
	}
};

class NUMA_col_wide_matrix_store: public NUMA_matrix_store
{
	NUMA_row_tall_matrix_store store;

	NUMA_col_wide_matrix_store(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type): NUMA_matrix_store(nrow, ncol, type), store(
				ncol, nrow, num_nodes, type) {
	}

	/*
	 * The constructed matrix store will be the transpose of _store.
	 */
	NUMA_col_wide_matrix_store(
			const NUMA_row_tall_matrix_store &_store): NUMA_matrix_store(
				_store.get_num_cols(), _store.get_num_rows(),
				_store.get_type(), _store.data_id), store(_store) {
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

	virtual vec_store::const_ptr conv2vec() const {
		return store.conv2vec();
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

	const char *get(size_t row_idx, size_t col_idx) const {
		return store.get(col_idx, row_idx);
	}

	char *get(size_t row_idx, size_t col_idx) {
		return store.get(col_idx, row_idx);
	}

	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const {
		auto cols = store.get_cols(idxs);
		if (cols == NULL)
			return matrix_store::const_ptr();
		else
			return cols->transpose();
	}

	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols);
	virtual std::shared_ptr<const local_matrix_store> get_portion(size_t id) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);
	virtual int get_portion_node_id(size_t id) const {
		return store.get_portion_node_id(id);
	}

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_COL;
	}

	virtual bool resize(size_t num_rows, size_t num_cols) {
		bool ret = store.resize(num_cols, num_rows);
		if (!ret)
			return false;
		return matrix_store::resize(num_rows, num_cols);
	}

	virtual matrix_store::const_ptr transpose() const;
};

}

}

#endif
