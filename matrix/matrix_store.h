#ifndef __MATRIX_STORE_H__
#define __MATRIX_STORE_H__

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
#include <atomic>
#include <unordered_map>

#include "safs_file.h"

#include "matrix_header.h"
#include "generic_type.h"

namespace fm
{

class set_operate;

namespace detail
{

class portion_compute;
class local_matrix_store;
class vec_store;

typedef std::pair<bool, std::shared_ptr<const local_matrix_store> > async_cres_t;

const size_t INVALID_MAT_ID = std::numeric_limits<size_t>::max();

class data_id_t
{
	size_t id;
	/*
	 * This counts the number of other virtual matrices that take this matrix
	 * as the direct input in the DAG.
	 * Multiple matrix_store may reference to the same matrix or the same matrix
	 * data, so we should keep the reference in this data structure instead of
	 * in matrix_store directly.
	 * The reference set contains data id of other virtual matrices. The reference
	 * counts how many times this matrix is needed to create the data of other
	 * matrices.
	 */
	std::set<size_t> refs;

	data_id_t(size_t id) {
		this->id = id;
	}
public:
	typedef std::shared_ptr<data_id_t> ptr;

	static ptr create(size_t id) {
		return ptr(new data_id_t(id));
	}

	size_t get_id() const {
		return id;
	}

	void inc_ref(size_t data_id) {
		if (data_id != INVALID_MAT_ID)
			refs.insert(data_id);
	}
	void reset_ref() {
		refs.clear();
	}
	size_t get_ref() const {
		return refs.size();
	}
};

class matrix_store
{
	size_t nrow;
	size_t ncol;
	bool in_mem;
	// This is kind of redundant because we can always get the entry size
	// from the type. However, getting the entry size of the type requires
	// to call a virtual method. Storing the entry size here can avoid
	// the function call. It doesn't increase the size of the data structure
	// due to the data alignment by the compiler.
	int entry_size;
	// The type is a reference. It makes the dense matrix object uncopiable.
	// Maybe this is what we want.
	const scalar_type &type;
	/*
	 * This indicates whether or not we cache a portion in each worker thread.
	 * By default, this is enabled.
	 */
	bool cache_portion;
public:
	typedef std::shared_ptr<matrix_store> ptr;
	typedef std::shared_ptr<const matrix_store> const_ptr;

	static std::atomic<size_t> mat_counter;
	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, int num_nodes, bool in_mem,
			safs::safs_file_group::ptr group = NULL);

	matrix_store(size_t nrow, size_t ncol, bool in_mem,
			const scalar_type &_type);

	template<class T>
	void init_randu(T min, T max) {
		std::shared_ptr<const set_operate> op = create_urand_init<T>(min, max);
		set_data(*op);
	}

	template<class T>
	void init_randn(T mean, T var) {
		std::shared_ptr<const set_operate> op = create_nrand_init<T>(mean, var);
		set_data(*op);
	}

	virtual ~matrix_store() {
	}

	virtual bool resize(size_t num_rows, size_t num_cols) {
		this->nrow = num_rows;
		this->ncol = num_cols;
		return true;
	}

	size_t get_num_rows() const {
		return nrow;
	}

	size_t get_num_cols() const {
		return ncol;
	}

	size_t get_entry_size() const {
		return entry_size;
	}

	const scalar_type &get_type() const {
		return type;
	}

	bool is_in_mem() const {
		return in_mem;
	}

	// A matrix can represent a sparse matrix.
	// By default, we assume all matrices are dense.
	virtual bool is_sparse() const {
		return false;
	}

	virtual int get_num_nodes() const {
		return -1;
	}

	virtual void inc_dag_ref(size_t data_id) {
	}
	virtual void reset_dag_ref() {
	}
	virtual size_t get_dag_ref() const {
		return 0;
	}

	/*
	 * The shape of a matrix: a tall matrix or a wide matrix.
	 * We care about the shape of a large matrix. We deal with a tall matrix
	 * different from a wide matrix.
	 */
	virtual bool is_wide() const {
		return get_num_cols() > get_num_rows();
	}

	/*
	 * Data Id is used to identify the data in a matrix.
	 * When a matrix is transposed or move to a different storage memory or
	 * converted into a different data layout, it should have the same data
	 * Id.
	 */
	virtual size_t get_data_id() const = 0;

	/*
	 * Test if this matrix shares the same data as the other matrix.
	 */
	virtual bool share_data(const matrix_store &store) const;

	/*
	 * This method gets underlying materialized matrix IDs and the number of
	 * elements in each of these materialized matrices.
	 */
	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const = 0;
	virtual std::string get_name() const = 0;

	virtual matrix_layout_t store_layout() const = 0;

	virtual void reset_data();
	virtual void set_data(const set_operate &op);

	virtual matrix_store::const_ptr transpose() const = 0;

	/*
	 * When matrix data is move to faster memory, data is moved in one chunk
	 * at a time. The chunk size is defined by a specific implementation of
	 * matrix store. Each chunk is assigned with an identifier, which is
	 * defined sequentially.
	 */
	size_t get_num_portions() const;
	virtual std::pair<size_t, size_t> get_portion_size() const = 0;
	/*
	 * This method gets a portion of data from the matrix asynchronously.
	 * We only need one version for asynchronous access to portions because
	 * the returned portion from this method is a copy from the original
	 * matrix and modification on the portion doesn't change the data in
	 * the original matrix.
	 *
	 * When a local matrix store is returned, it's not guaranteed that the
	 * data in the local matrix store is valid. A status in the returned value
	 * indicates whether the data is valid. If the data is invalid when it's
	 * returned, the computation passed to
	 * this method is invoked when the portion of data is loaded
	 * in memory. During the time between returning from the method and
	 * the portion of data becomes available, it's users' responsibility
	 * of keep the local matrix store alive.
	 */
	virtual async_cres_t get_portion_async(size_t start_row, size_t start_col,
			size_t num_rows, size_t num_cols,
			std::shared_ptr<portion_compute> compute) const = 0;
	/*
	 * These versions fetches the portion of data. It's guaranteed that
	 * the data in the returned local matrix store is valid.
	 */
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const = 0;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) = 0;
	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id);
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t id) const;
	/*
	 * This gets the node Id of the specified portion.
	 */
	virtual int get_portion_node_id(size_t id) const = 0;
	virtual void write_portion_async(
			std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col) = 0;

	virtual std::shared_ptr<const vec_store> conv2vec() const {
		return std::shared_ptr<const vec_store>();
	}

	/*
	 * Get rows/columns within [start, end).
	 */
	virtual matrix_store::const_ptr get_cols(off_t start, off_t end) const;
	virtual matrix_store::const_ptr get_rows(off_t start, off_t end) const;
	virtual matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const;
	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const {
		return matrix_store::const_ptr();
	}

	virtual bool is_virtual() const {
		return false;
	}
	virtual bool is_sink() const {
		return false;
	}
	virtual void materialize_self() const {
	}

	/*
	 * This allow users to enable/disable caching data portions.
	 * It's used by EM matrix and mapply virtual matrix.
	 */
	virtual void set_cache_portion(bool cache_portion) {
		this->cache_portion = cache_portion;
	}
	virtual bool is_cache_portion() const {
		return cache_portion;
	}
	virtual void set_prefetches(size_t num, std::pair<size_t, size_t> range) {
	}
};

/*
 * This matrix stream is designed to deal with the case that we know exactly
 * where data should be written but the data may come from multiple threads
 * in a non-sequential order.
 */
class matrix_stream
{
public:
	typedef std::shared_ptr<matrix_stream> ptr;

	static ptr create(matrix_store::ptr store);

	virtual void write_async(std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col) = 0;
	virtual bool is_complete() const = 0;
	virtual const matrix_store &get_mat() const = 0;
};

/*
 * This class is designed for the case that we don't know where exactly data
 * should be written to, but we know relative location of each write.
 * Each write is identified by sequence numbers.
 */
class matrix_append
{
	matrix_store::ptr res;

	// This is used to fill the slot in the queue that doesn't have the data
	// from a portion.
	std::shared_ptr<const local_matrix_store> empty_portion;

	spin_lock lock;
	// A queue that contains portions appended out of order. The location
	// of a portion in the queue is determined by the difference between
	// its sequence id and `last_append'.
	std::deque<std::shared_ptr<const local_matrix_store>> q;
	// The sequence number of the last append.
	off_t last_append;
	// The total number of elements written to the result matrix.
	std::atomic<size_t> written_eles;

	matrix_append(matrix_store::ptr store);
public:
	typedef std::shared_ptr<matrix_append> ptr;

	~matrix_append();

	static ptr create(matrix_store::ptr store) {
		return ptr(new matrix_append(store));
	}

	void write_async(std::shared_ptr<const local_matrix_store> portion,
			off_t seq_id);

	void flush();

	size_t get_written_eles() const {
		return written_eles;
	}
};

}

}

#endif
