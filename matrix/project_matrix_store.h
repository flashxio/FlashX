#ifndef __FM_PROJECT_MATRIX_H__
#define __FM_PROJECT_MATRIX_H__

/*
 * Copyright 2016 Open Connectome Project (http://openconnecto.me)
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

#include "mem_matrix_store.h"
#include "local_matrix_store.h"

namespace fm
{

namespace detail
{

class sparse_project_matrix_store: public mem_matrix_store
{
public:
	struct nnz_idx {
		off_t row_idx;
		off_t col_idx;
		bool operator==(const nnz_idx &e) const {
			return row_idx == e.row_idx && col_idx == e.col_idx;
		}
	};

	/*
	 * This orders the nnz in row-major order in a tall matrix.
	 */
	struct tall_row_first_comp {
		bool operator()(const nnz_idx &e1, const nnz_idx &e2) const {
			if (e1.row_idx < e2.row_idx)
				return true;
			else if (e1.row_idx > e2.row_idx)
				return false;
			else
				return e1.col_idx < e2.col_idx;
		}
	};

	/*
	 * This orders the nnz in col-major order in a tall matrix.
	 * The nnz is first partitioned in chunks.
	 */
	struct tall_col_first_comp {
		bool operator()(const nnz_idx &e1, const nnz_idx &e2) const {
			off_t portion1 = e1.row_idx / mem_matrix_store::CHUNK_SIZE;
			off_t portion2 = e2.row_idx / mem_matrix_store::CHUNK_SIZE;
			if (portion1 < portion2)
				return true;
			else if (portion1 > portion2)
				return false;
			else if (e1.col_idx < e2.col_idx)
				return true;
			else if (e1.col_idx > e2.col_idx)
				return false;
			else
				return e1.row_idx < e2.row_idx;
		}
	};

	/*
	 * This orders the nnz in row-major order in a wide matrix.
	 * The nnz is first partitioned in chunks.
	 */
	struct wide_row_first_comp {
		bool operator()(const nnz_idx &e1, const nnz_idx &e2) const {
			off_t portion1 = e1.col_idx / mem_matrix_store::CHUNK_SIZE;
			off_t portion2 = e2.col_idx / mem_matrix_store::CHUNK_SIZE;
			if (portion1 < portion2)
				return true;
			else if (portion1 > portion2)
				return false;
			else if (e1.row_idx < e2.row_idx)
				return true;
			else if (e1.row_idx > e2.row_idx)
				return false;
			else
				return e1.col_idx < e2.col_idx;
		}
	};

	/*
	 * This orders the nnz in col-major order in a wide matrix.
	 */
	struct wide_col_first_comp {
		bool operator()(const nnz_idx &e1, const nnz_idx &e2) const {
			if (e1.col_idx < e2.col_idx)
				return true;
			else if (e1.col_idx > e2.col_idx)
				return false;
			else
				return e1.row_idx < e2.row_idx;
		}
	};

private:
	std::vector<nnz_idx> nnz_idxs;
	mem_col_matrix_store::const_ptr vals;
	matrix_layout_t layout;

	sparse_project_matrix_store(size_t nrow, size_t ncol,
			matrix_layout_t layout, const scalar_type &type);
public:
	typedef std::shared_ptr<sparse_project_matrix_store> ptr;
	typedef std::shared_ptr<const sparse_project_matrix_store> const_ptr;

	static ptr create_sparse_rand(size_t nrow, size_t ncol,
			matrix_layout_t layout, const scalar_type &type, double density);

	virtual char *get(size_t row, size_t col) {
		return NULL;
	}

	virtual bool write2file(const std::string &file_name) const {
		throw unsupported_exception(
				"don't support writing a sparse matrix to a file");
	}
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) {
		throw unsupported_exception("a sparse matrix is read-only");
	}
	virtual int get_portion_node_id(size_t id) const {
		return -1;
	}

	virtual bool is_sparse() const {
		return true;
	}

	virtual const char *get(size_t row, size_t col) const;

	virtual std::string get_name() const {
		return (boost::format("sparse_mat-%1%(%2%,%3%)") % get_mat_id()
				% get_num_rows() % get_num_cols()).str();
	}

	virtual matrix_layout_t store_layout() const {
		return layout;
	}

	virtual matrix_store::const_ptr transpose() const;

	virtual async_cres_t get_portion_async(size_t start_row, size_t start_col,
			size_t num_rows, size_t num_cols,
			std::shared_ptr<portion_compute> compute) const {
		// The matrix is always in-memory.
		return async_cres_t(true, get_portion(start_row, start_col, num_rows,
					num_cols));
	}
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;

	virtual void write_portion_async(
			std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col);

	size_t get_nnz() const {
		return nnz_idxs.size();
	}

	matrix_store::const_ptr conv_dense() const;
};

class local_sparse_matrix_store: public local_matrix_store
{
	struct col_comp {
		bool operator()(const sparse_project_matrix_store::nnz_idx &e1,
				const sparse_project_matrix_store::nnz_idx &e2) const {
			return e1.col_idx < e2.col_idx;
		}
	};
	struct row_comp {
		bool operator()(const sparse_project_matrix_store::nnz_idx &e1,
				const sparse_project_matrix_store::nnz_idx &e2) const {
			return e1.row_idx < e2.row_idx;
		}
	};

	std::vector<sparse_project_matrix_store::nnz_idx> local_idxs;
	local_matrix_store::const_ptr vals;
	matrix_layout_t layout;
public:
	typedef std::shared_ptr<local_sparse_matrix_store> ptr;
	typedef std::shared_ptr<const local_sparse_matrix_store> const_ptr;

	local_sparse_matrix_store(off_t global_start_row, off_t global_start_col,
			size_t nrow, size_t ncol,
			const std::vector<sparse_project_matrix_store::nnz_idx> &local_idxs,
			local_matrix_store::const_ptr vals, matrix_layout_t layout);

	virtual bool read_only() const {
		return true;
	}
	virtual char *get_raw_arr() {
		return NULL;
	}
	virtual const char *get_raw_arr() const {
		return NULL;
	}
	virtual char *get(size_t row, size_t col) {
		return NULL;
	}
	virtual bool copy_from(const local_matrix_store &store) {
		return false;
	}
	virtual matrix_store::ptr transpose() {
		return matrix_store::ptr();
	}
	virtual size_t get_all_rows(std::vector<char *> &rows) {
		throw unsupported_exception(
				"don't support getting all rows of a local sparse matrix");
	}
	virtual size_t get_all_cols(std::vector<char *> &cols) {
		throw unsupported_exception(
				"don't support getting all cols of a local sparse matrix");
	}
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) {
		throw unsupported_exception("a sparse matrix is read-only");
	}

	virtual int get_portion_node_id(size_t id) const {
		return -1;
	}

	virtual matrix_layout_t store_layout() const {
		return layout;
	}
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols);
	virtual local_matrix_store::ptr conv2(matrix_layout_t layout) const;
	virtual size_t get_all_rows(std::vector<const char *> &rows) const;
	virtual size_t get_all_cols(std::vector<const char *> &cols) const;
	virtual const char *get(size_t row, size_t col) const;
	virtual matrix_store::const_ptr transpose() const;
	const char *get_col_nnz(off_t col_idx, std::vector<off_t> &row_idxs) const;
	const char *get_row_nnz(off_t col_idx, std::vector<off_t> &row_idxs) const;
	size_t get_nnz() const {
		return local_idxs.size();
	}
};

}

}

#endif
