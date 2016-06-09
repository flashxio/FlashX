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
	const size_t data_id;

	matrix_layout_t layout;
	file_holder::ptr holder;
	io_set::ptr ios;

	/*
	 * This indicates whether or not we cache a portion in each worker thread.
	 * By default, this is enabled.
	 */
	bool cache_portion;

	/*
	 * These two fields are used for sub matrix.
	 * They indicates the actual number of rows and columns stored on disks.
	 * In contrast, get_num_rows() and get_num_cols() are #rows and columns
	 * exposed to users.
	 */
	size_t orig_num_rows;
	size_t orig_num_cols;

	size_t get_orig_num_rows() const {
		return orig_num_rows;
	}

	size_t get_orig_num_cols() const {
		return orig_num_cols;
	}

	EM_matrix_store(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, safs::safs_file_group::ptr group);
	EM_matrix_store(file_holder::ptr holder, io_set::ptr ios, size_t nrow,
			size_t ncol, size_t orig_nrow, size_t orig_ncol,
			matrix_layout_t layout, const scalar_type &type, size_t _data_id);
public:
	static const size_t CHUNK_SIZE;

	typedef std::shared_ptr<EM_matrix_store> ptr;
	typedef std::shared_ptr<const EM_matrix_store> const_ptr;

	static ptr create(const std::string &mat_file);

	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout,
			const scalar_type &type, safs::safs_file_group::ptr group = NULL) {
		return ptr(new EM_matrix_store(nrow, ncol, layout, type, group));
	}

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

	virtual void set_cache_portion(bool cache_portion) {
		this->cache_portion = cache_portion;
	}

	bool is_cache_portion() const {
		return cache_portion;
	}

	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		std::unordered_map<size_t, size_t> ret;
		ret.insert(std::pair<size_t, size_t>(data_id,
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
	virtual async_res_t get_portion_async(size_t start_row, size_t start_col,
			size_t num_rows, size_t num_cols,
			portion_compute::ptr compute);
	virtual void write_portion_async(
			std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col);
	void wait4complete();

	virtual matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const;
	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const;
	virtual std::shared_ptr<const vec_store> get_col_vec(off_t idx) const;
	virtual std::shared_ptr<const vec_store> get_row_vec(off_t idx) const;

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
};

/*
 * This class helps to stream data to the output dense matrix on disks.
 * This assumes that the data comes from multiple threads without a specific
 * order, but they are written to the locations close to each other on disks.
 * Once we order the incoming data, we can stream data to disks sequentially.
 * This data structure is shared by multiple threads.
 */
class EM_matrix_stream
{
	matrix_store::ptr mat;

	/*
	 * This stores data in a portion of a matrix.
	 */
	class filled_portion
	{
		std::shared_ptr<local_matrix_store> data;
		// The number of filled elements.
		std::atomic<size_t> num_filled;
	public:
		typedef std::shared_ptr<filled_portion> ptr;

		filled_portion(std::shared_ptr<local_matrix_store> data) {
			this->data = data;
			num_filled = 0;
		}
		// If the write completely fill the local buffer, it returns true.
		bool write(std::shared_ptr<const local_matrix_store> portion,
				off_t global_start_row, off_t global_start_col);

		std::shared_ptr<const local_matrix_store> get_whole_portion() const {
			return data;
		}
	};

	class portion_queue
	{
		// All local matrices have the same size as the portion size of
		// the EM matrix.
		std::map<off_t, std::shared_ptr<const local_matrix_store> > bufs;
		// The last portion that has been flushed to disks.
		size_t num_flushed;

		const size_t portion_size;
		const bool is_wide;
		const matrix_layout_t out_layout;
	public:
		portion_queue(size_t _portion_size, bool _is_wide,
				matrix_layout_t layout): portion_size(_portion_size), is_wide(
					_is_wide), out_layout(layout) {
			num_flushed = 0;
		}

		void add(std::shared_ptr<const local_matrix_store> lmat,
				off_t start_row, off_t start_col);
		std::shared_ptr<const local_matrix_store> pop_contig(size_t num_portions);

		size_t get_num_flushed() const {
			return num_flushed;
		}
		size_t get_buf_size() const {
			return bufs.size();
		}
	};

	pthread_spinlock_t lock;
	// This keeps the buffers with partial data in EM matrix portions.
	// If an EM matrix portion is complete, the portion is flushed to disks
	// and it is deleted from the hashtable.
	std::unordered_map<off_t, filled_portion::ptr> portion_bufs;

	std::shared_ptr<portion_queue> portion_q;
	// The number of portions we need to write to disks together.
	// It determines the minimal I/O size.
	size_t min_io_portions;

	EM_matrix_stream(matrix_store::ptr mat);
	void write_portion_async(std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col);
public:
	typedef std::shared_ptr<EM_matrix_stream> ptr;

	static ptr create(matrix_store::ptr mat) {
		return ptr(new EM_matrix_stream(mat));
	}

	~EM_matrix_stream();

	void write_async(std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col);
	void flush();
	bool is_complete() const;
};

}

}

#endif
