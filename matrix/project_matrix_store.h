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

		nnz_idx() {
			row_idx = 0;
			col_idx = 0;
		}

		nnz_idx(off_t row_idx, off_t col_idx) {
			this->row_idx = row_idx;
			this->col_idx = col_idx;
		}

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

class lsparse_col_matrix_store: public lvirtual_col_matrix_store
{
	struct col_first_comp {
		bool operator()(const sparse_project_matrix_store::nnz_idx &e1,
				const sparse_project_matrix_store::nnz_idx &e2) const {
			if (e1.col_idx < e2.col_idx)
				return true;
			else if (e1.col_idx > e2.col_idx)
				return false;
			else
				return e1.row_idx < e2.row_idx;
		}
	};

	std::vector<sparse_project_matrix_store::nnz_idx> local_idxs;
	local_matrix_store::const_ptr vals;

	local_col_matrix_store::const_ptr buf;
public:
	typedef std::shared_ptr<lsparse_col_matrix_store> ptr;
	typedef std::shared_ptr<const lsparse_col_matrix_store> const_ptr;

	lsparse_col_matrix_store(off_t start_row, off_t start_col, size_t nrow, size_t ncol,
			const std::vector<sparse_project_matrix_store::nnz_idx> &local_idxs,
			local_matrix_store::const_ptr vals): lvirtual_col_matrix_store(
				start_row, start_col, nrow, ncol, vals->get_type(),
				vals->get_node_id()) {
		this->local_idxs = local_idxs;
		this->vals = vals;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		off_t gstart_row = get_global_start_row() + local_start_row;
		off_t gstart_col = get_global_start_col() + local_start_col;
		// If the buffered matrix doesn't have the data for the specified
		// location.
		if (buf && (buf->get_global_start_row() != gstart_row
				|| buf->get_global_start_col() != gstart_col
				|| buf->get_num_rows() != local_num_rows
				|| buf->get_num_cols() != local_num_cols))
			buf = NULL;
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		// If the buffered matrix doesn't have the data for the entire local
		// matrix.
		if (buf && (buf->get_global_start_row() != get_global_start_row()
				|| buf->get_global_start_col() != get_global_start_col()
				|| buf->get_num_rows() != get_num_rows()
				|| buf->get_num_cols() != get_num_cols()))
			buf = NULL;
		local_matrix_store::reset_size();
	}

	using lvirtual_col_matrix_store::get_raw_arr;
	virtual const char *get_raw_arr() const {
		materialize_self();
		return buf->get_raw_arr();
	}

	using lvirtual_col_matrix_store::transpose;
	virtual matrix_store::const_ptr transpose() const;

	using lvirtual_col_matrix_store::get_col;
	virtual const char *get_col(size_t col) const {
		materialize_self();
		return buf->get_col(col);
	}

	virtual local_matrix_store::const_ptr get_portion(
			size_t local_start_row, size_t local_start_col, size_t num_rows,
			size_t num_cols) const;

	virtual void materialize_self() const;

	const char *get_col_nnz(off_t col_idx, std::vector<off_t> &row_idxs) const;
	size_t get_nnz() const {
		return local_idxs.size();
	}
};

class lsparse_row_matrix_store: public lvirtual_row_matrix_store
{
	struct row_first_comp {
		bool operator()(const sparse_project_matrix_store::nnz_idx &e1,
				const sparse_project_matrix_store::nnz_idx &e2) const {
			if (e1.row_idx < e2.row_idx)
				return true;
			else if (e1.row_idx > e2.row_idx)
				return false;
			else
				return e1.col_idx < e2.col_idx;
		}
	};

	std::vector<sparse_project_matrix_store::nnz_idx> local_idxs;
	local_matrix_store::const_ptr vals;

	local_row_matrix_store::const_ptr buf;
public:
	typedef std::shared_ptr<lsparse_row_matrix_store> ptr;
	typedef std::shared_ptr<const lsparse_row_matrix_store> const_ptr;

	lsparse_row_matrix_store(off_t start_row, off_t start_col, size_t nrow, size_t ncol,
			const std::vector<sparse_project_matrix_store::nnz_idx> &local_idxs,
			local_matrix_store::const_ptr vals): lvirtual_row_matrix_store(
				start_row, start_col, nrow, ncol, vals->get_type(),
				vals->get_node_id()) {
		this->local_idxs = local_idxs;
		this->vals = vals;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		off_t gstart_row = get_global_start_row() + local_start_row;
		off_t gstart_col = get_global_start_col() + local_start_col;
		// If the buffered matrix doesn't have the data for the specified
		// location.
		if (buf && (buf->get_global_start_row() != gstart_row
				|| buf->get_global_start_col() != gstart_col
				|| buf->get_num_rows() != local_num_rows
				|| buf->get_num_cols() != local_num_cols))
			buf = NULL;
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		// If the buffered matrix doesn't have the data for the entire local
		// matrix.
		if (buf && (buf->get_global_start_row() != get_global_start_row()
				|| buf->get_global_start_col() != get_global_start_col()
				|| buf->get_num_rows() != get_num_rows()
				|| buf->get_num_cols() != get_num_cols()))
			buf = NULL;
		local_matrix_store::reset_size();
	}

	using lvirtual_row_matrix_store::get_raw_arr;
	virtual const char *get_raw_arr() const {
		materialize_self();
		return buf->get_raw_arr();
	}

	using lvirtual_row_matrix_store::transpose;
	virtual matrix_store::const_ptr transpose() const;

	using lvirtual_row_matrix_store::get_row;
	virtual const char *get_row(size_t row) const {
		materialize_self();
		return buf->get_row(row);
	}

	virtual local_matrix_store::const_ptr get_portion(
			size_t local_start_row, size_t local_start_col, size_t num_rows,
			size_t num_cols) const;
	virtual void materialize_self() const;

	size_t get_nnz() const {
		return local_idxs.size();
	}

	const char *get_row_nnz(off_t row_idx, std::vector<off_t> &col_idxs) const;
	// This returns rows in the specified range.
	// The rows must be stored contiguously. Otherwise, it'll return NULL.
	// The local store can be resized. The indexes of non-zero entries are
	// relative in the resized store.
	const char *get_rows_nnz(off_t start_row, off_t end_row,
			std::vector<sparse_project_matrix_store::nnz_idx> &idxs) const;
};

}

}

#endif
