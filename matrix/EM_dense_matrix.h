#ifndef __EM_DENSE_MATRIX_H__
#define __EM_DENSE_MATRIX_H__

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
#include <boost/format.hpp>

#include "log.h"

#include "bulk_operate.h"
#include "matrix_store.h"
#include "EM_object.h"
#include "mem_worker_thread.h"

namespace fm
{

namespace detail
{

class local_matrix_store;
class mem_matrix_store;
class EM_matrix_stream;

class EM_matrix_store: public matrix_store, public EM_object
{
	/*
	 * The difference between the two identifiers are:
	 * `mat_id' identifies the matrix data structure. Whenever the matrix
	 * is shallow copied or transposed, `mat_id' changes.
	 * `data_id' identifies the content in a matrix.
	 * So when a matrix is transposed or shallow copied, it should share
	 * the same data id.
	 */
	const size_t mat_id;
	const data_id_t::ptr data_id;

	matrix_layout_t layout;
	file_holder::ptr holder;
	io_set::ptr ios;

	std::shared_ptr<EM_matrix_stream> stream;

	/*
	 * These two fields are used for sub matrix.
	 * They indicates the actual number of rows and columns stored on disks.
	 * In contrast, get_num_rows() and get_num_cols() are #rows and columns
	 * exposed to users.
	 */
	size_t orig_num_rows;
	size_t orig_num_cols;

	size_t num_prefetches;
	std::pair<size_t, size_t> prefetch_range;

	size_t get_orig_num_rows() const {
		return orig_num_rows;
	}

	size_t get_orig_num_cols() const {
		return orig_num_cols;
	}

	EM_matrix_store(file_holder::ptr holder, io_set::ptr ios, size_t nrow,
			size_t ncol, matrix_layout_t layout, const scalar_type &type,
			safs::safs_file_group::ptr group);
	EM_matrix_store(file_holder::ptr holder, io_set::ptr ios, size_t nrow,
			size_t ncol, size_t orig_nrow, size_t orig_ncol,
			matrix_layout_t layout, const scalar_type &type,
			data_id_t::ptr data_id = data_id_t::create(mat_counter++));
	void _write_portion_async(
			std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col);
	virtual std::pair<size_t, size_t> get_orig_portion_size() const;
public:
	static const size_t CHUNK_SIZE;

	typedef std::shared_ptr<EM_matrix_store> ptr;
	typedef std::shared_ptr<const EM_matrix_store> const_ptr;

	static ptr create(const std::string &mat_file);
	static ptr create(file_holder::ptr holder, io_set::ptr ios, size_t nrow,
			size_t ncol, matrix_layout_t layout, const scalar_type &type,
			safs::safs_file_group::ptr group = NULL) {
		return ptr(new EM_matrix_store(holder, ios, nrow, ncol, nrow, ncol,
					layout, type));
	}

	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, safs::safs_file_group::ptr group = NULL);

	/*
	 * This is to load the matrix data from an external file and create
	 * an EM matrix on SAFS.
	 */
	static ptr load(const std::string &ext_file, size_t nrow, size_t ncol,
			matrix_layout_t layout, const scalar_type &type,
			safs::safs_file_group::ptr group = NULL);

	static ptr cast(matrix_store::ptr store) {
		return std::dynamic_pointer_cast<EM_matrix_store>(store);
	}

	static const_ptr cast(matrix_store::const_ptr store) {
		return std::dynamic_pointer_cast<const EM_matrix_store>(store);
	}

	EM_matrix_store::const_ptr shallow_copy() const {
		return EM_matrix_store::const_ptr(new EM_matrix_store(holder, ios,
					get_num_rows(), get_num_cols(), orig_num_rows, orig_num_cols,
					store_layout(), get_type(), data_id));
	}

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

	virtual void set_prefetches(size_t num, std::pair<size_t, size_t> range) {
		this->num_prefetches = num;
		this->prefetch_range = range;
	}

	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		std::unordered_map<size_t, size_t> ret;
		ret.insert(std::pair<size_t, size_t>(data_id->get_id(),
					get_num_rows() * get_num_cols()));
		return ret;
	}
	virtual std::string get_name() const {
		return (boost::format("EM_mat-%1%(%2%,%3%)") % mat_id % get_num_rows()
			% get_num_cols()).str();
	}

	virtual matrix_layout_t store_layout() const {
		return layout;
	}

	virtual matrix_store::const_ptr transpose() const;

	virtual std::vector<safs::io_interface::ptr> create_ios() const;

	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols);
	virtual int get_portion_node_id(size_t id) const {
		return -1;
	}

	virtual std::pair<size_t, size_t> get_portion_size() const;
	virtual async_cres_t get_portion_async(size_t start_row, size_t start_col,
			size_t num_rows, size_t num_cols,
			portion_compute::ptr compute) const;
	virtual void write_portion_async(
			std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col);
	/*
	 * This starts a stream that collects users writes and writes them to disks
	 * with a large I/O. This should be called in the main thread.
	 */
	void start_stream();
	/*
	 * This ends a stream and makes sure all users' data is written to disks.
	 * This should be called in the main thread.
	 */
	void end_stream();
	void wait4complete() const;
	void wait4complete(size_t num_ios) const;

	virtual std::shared_ptr<const vec_store> conv2vec() const;

	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const;

	/*
	 * Set this matrix persistent in SAFS, so that even if there isn't
	 * a reference to the matrix, its data still stored in SAFS.
	 * This method isn't thread-safe.
	 */
	bool set_persistent(const std::string &name) const;
	/*
	 * Unset the persistency of the matrix in SAFS, so that the matrix file
	 * is deleted after all references to the matrix are gone.
	 */
	void unset_persistent() const;

	friend class EM_matrix_stream;
};

/*
 * This class helps to stream data to the output dense matrix on disks.
 * This assumes that the data comes from multiple threads without a specific
 * order, but they are written to the locations close to each other on disks.
 * Once we order the incoming data, we can stream data to disks sequentially.
 * This data structure is shared by multiple threads.
 */
class EM_matrix_stream: public matrix_stream
{
	EM_matrix_store::ptr mat;

	/*
	 * This stores data for multiple portions of a matrix.
	 */
	class filled_portion
	{
		std::shared_ptr<local_matrix_store> data;
		std::vector<std::shared_ptr<local_matrix_store> > portions;
		// The number of filled elements.
		std::atomic<size_t> num_filled;
		// The total number of elements to fill.
		size_t tot_num_eles;

		size_t portion_size;
		bool is_wide;
	public:
		typedef std::shared_ptr<filled_portion> ptr;

		filled_portion(local_raw_array &arr, off_t start_row, off_t start_col,
				size_t num_rows, size_t num_cols, const scalar_type &type,
				matrix_layout_t layout, size_t portion_size, bool is_wide);
		std::shared_ptr<local_matrix_store> create_portion(char *arr,
				off_t start_row, off_t start_col, size_t num_rows,
				size_t num_cols, matrix_layout_t layout, const scalar_type &type);
		std::shared_ptr<local_matrix_store> get_portion(off_t start_row,
				off_t start_col);

		// If the write completely fill the local buffer, it returns true.
		bool write(std::shared_ptr<const local_matrix_store> portion,
				off_t global_start_row, off_t global_start_col);

		std::shared_ptr<const local_matrix_store> get_whole_portion() const {
			return data;
		}
	};

	pthread_spinlock_t lock;
	// This keeps the buffers with partial data in EM matrix portions.
	// If an EM matrix portion is complete, the portion is flushed to disks
	// and it is deleted from the hashtable.
	std::unordered_map<off_t, filled_portion::ptr> portion_bufs;

	// The number of portions we need to write to disks together.
	// It determines the minimal I/O size.
	size_t min_io_portions;

	EM_matrix_stream(EM_matrix_store::ptr mat);
public:
	typedef std::shared_ptr<EM_matrix_stream> ptr;

	static ptr create(EM_matrix_store::ptr mat) {
		return ptr(new EM_matrix_stream(mat));
	}

	virtual ~EM_matrix_stream();

	virtual void write_async(std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col);
	virtual bool is_complete() const;

	virtual const matrix_store &get_mat() const {
		return *mat;
	}
};

}

}

#endif
