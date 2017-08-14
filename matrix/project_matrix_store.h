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

#include "virtual_matrix_store.h"
#include "local_matrix_store.h"

namespace fm
{

namespace detail
{

class sparse_project_matrix_store: public virtual_matrix_store
{
public:
	struct nz_idx {
		off_t row_idx;
		off_t col_idx;

		nz_idx() {
			row_idx = 0;
			col_idx = 0;
		}

		nz_idx(off_t row_idx, off_t col_idx) {
			this->row_idx = row_idx;
			this->col_idx = col_idx;
		}

		bool operator==(const nz_idx &e) const {
			return row_idx == e.row_idx && col_idx == e.col_idx;
		}
	};

	/*
	 * This orders the nnz in row-major order in a tall matrix.
	 */
	struct tall_row_first_comp {
		bool operator()(const nz_idx &e1, const nz_idx &e2) const {
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
		bool operator()(const nz_idx &e1, const nz_idx &e2) const {
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
		bool operator()(const nz_idx &e1, const nz_idx &e2) const {
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
		bool operator()(const nz_idx &e1, const nz_idx &e2) const {
			if (e1.col_idx < e2.col_idx)
				return true;
			else if (e1.col_idx > e2.col_idx)
				return false;
			else
				return e1.row_idx < e2.row_idx;
		}
	};

private:
	const size_t mat_id;
	std::vector<off_t> portion_offs;
	// This contains the global row and col indeces for non-zero entries.
	std::vector<nz_idx> nz_idxs;
	mem_col_matrix_store::ptr vals;
	matrix_layout_t layout;

	sparse_project_matrix_store(size_t nrow, size_t ncol,
			matrix_layout_t layout, const scalar_type &type);
	sparse_project_matrix_store(size_t nrow, size_t ncol,
			matrix_layout_t layout, const scalar_type &type, double density);
	bool is_entire_portion(size_t start_row, size_t start_col,
			size_t num_rows, size_t num_cols) const;
public:
	typedef std::shared_ptr<sparse_project_matrix_store> ptr;
	typedef std::shared_ptr<const sparse_project_matrix_store> const_ptr;

	static ptr create_sparse_rand(size_t nrow, size_t ncol,
			matrix_layout_t layout, const scalar_type &type, double density);

	virtual bool has_materialized() const {
		return false;
	}

	virtual size_t get_data_id() const {
		// TODO I should return a real data id.
		return INVALID_MAT_ID;
	}

	virtual char *get(size_t row, size_t col) {
		return NULL;
	}

	virtual int get_portion_node_id(size_t id) const {
		return -1;
	}

	virtual bool is_sparse() const {
		return true;
	}

	virtual const char *get(size_t row, size_t col) const;

	virtual std::string get_name() const {
		return (boost::format("sparse_mat-%1%(%2%,%3%)") % mat_id
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

	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const;
	virtual matrix_store::const_ptr materialize(bool in_mem,
			int num_nodes) const;
	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		std::unordered_map<size_t, size_t> ret;
		// TODO right now we only indicate the matrix. We set the number of
		// bytes to 0
		// We should also use data_id instead of mat_id.
		ret.insert(std::pair<size_t, size_t>(mat_id, 0));
		return ret;
	}
	virtual std::pair<size_t, size_t> get_portion_size() const;

	size_t get_nnz() const {
		return nz_idxs.size();
	}
};

class lsparse_col_matrix_store: public lvirtual_col_matrix_store
{
	struct col_first_comp {
		bool operator()(const sparse_project_matrix_store::nz_idx &e1,
				const sparse_project_matrix_store::nz_idx &e2) const {
			if (e1.col_idx < e2.col_idx)
				return true;
			else if (e1.col_idx > e2.col_idx)
				return false;
			else
				return e1.row_idx < e2.row_idx;
		}
	};

	// This contains the local row and col indices of non-zero entries
	// in the matrix.
	std::vector<sparse_project_matrix_store::nz_idx> local_idxs;
	local_matrix_store::const_ptr vals;

	local_col_matrix_store::const_ptr buf;
public:
	typedef std::shared_ptr<lsparse_col_matrix_store> ptr;
	typedef std::shared_ptr<const lsparse_col_matrix_store> const_ptr;

	lsparse_col_matrix_store(off_t start_row, off_t start_col, size_t nrow, size_t ncol,
			const std::vector<sparse_project_matrix_store::nz_idx> &local_idxs,
			local_matrix_store::const_ptr vals): lvirtual_col_matrix_store(
				start_row, start_col, nrow, ncol, vals->get_type(),
				vals->get_node_id()) {
		this->local_idxs = local_idxs;
		this->vals = vals;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		bool ret = local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
		if (!ret)
			return false;
		// If the buffered matrix doesn't have the data for the specified
		// location.
		if (buf && (buf->get_global_start_row() != get_global_start_row()
				|| buf->get_global_start_col() != get_global_start_col()
				|| buf->get_num_rows() != local_num_rows
				|| buf->get_num_cols() != local_num_cols))
			buf = NULL;
		return true;
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
	virtual local_matrix_store::ptr conv2(matrix_layout_t layout) const;

	// get the non-zero entries in a specified column.
	// `row_idxs' contains the local row indices.
	// If the local portion is reiszed, the row indices are relative to
	// the start row of the resized portion.
	const char *get_col_nnz(off_t col_idx, std::vector<off_t> &row_idxs) const;

	// This is mainly for testing.
	// We don't need to do anything after the local portion is resized.
	size_t get_nnz() const {
		return local_idxs.size();
	}
};

class lsparse_row_matrix_store: public lvirtual_row_matrix_store
{
	struct row_first_comp {
		bool operator()(const sparse_project_matrix_store::nz_idx &e1,
				const sparse_project_matrix_store::nz_idx &e2) const {
			if (e1.row_idx < e2.row_idx)
				return true;
			else if (e1.row_idx > e2.row_idx)
				return false;
			else
				return e1.col_idx < e2.col_idx;
		}
	};

	std::vector<sparse_project_matrix_store::nz_idx> local_idxs;
	local_matrix_store::const_ptr vals;

	local_row_matrix_store::const_ptr buf;
public:
	typedef std::shared_ptr<lsparse_row_matrix_store> ptr;
	typedef std::shared_ptr<const lsparse_row_matrix_store> const_ptr;

	lsparse_row_matrix_store(off_t start_row, off_t start_col, size_t nrow, size_t ncol,
			const std::vector<sparse_project_matrix_store::nz_idx> &local_idxs,
			local_matrix_store::const_ptr vals): lvirtual_row_matrix_store(
				start_row, start_col, nrow, ncol, vals->get_type(),
				vals->get_node_id()) {
		this->local_idxs = local_idxs;
		this->vals = vals;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		bool ret = local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
		if (!ret)
			return false;
		// If the buffered matrix doesn't have the data for the specified
		// location.
		if (buf && (buf->get_global_start_row() != get_global_start_row()
				|| buf->get_global_start_col() != get_global_start_col()
				|| buf->get_num_rows() != local_num_rows
				|| buf->get_num_cols() != local_num_cols))
			buf = NULL;
		return true;
	}
	virtual void reset_size() {
		local_matrix_store::reset_size();
		// If the buffered matrix doesn't have the data for the entire local
		// matrix.
		if (buf && (buf->get_global_start_row() != get_global_start_row()
				|| buf->get_global_start_col() != get_global_start_col()
				|| buf->get_num_rows() != get_num_rows()
				|| buf->get_num_cols() != get_num_cols()))
			buf = NULL;
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
	virtual local_matrix_store::ptr conv2(matrix_layout_t layout) const;

	// This is mainly for testing.
	// We don't need to do anything after the local portion is resized.
	size_t get_nnz() const {
		return local_idxs.size();
	}

	// get the non-zero entries in a specified row.
	// `col_idxs' contains the local col indices.
	// If the local portion is reiszed, the col indices are relative to
	// the start col of the resized portion.
	const char *get_row_nnz(off_t row_idx, std::vector<off_t> &col_idxs) const;
	// This returns rows in the specified range.
	// The rows must be stored contiguously. Otherwise, it'll return NULL.
	// The local store can be resized. The indexes of non-zero entries are
	// relative in the resized store.
	const char *get_rows_nnz(off_t start_row, off_t end_row,
			std::vector<sparse_project_matrix_store::nz_idx> &idxs) const;
};

}

}

#endif
