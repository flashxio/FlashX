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

#include "EM_vector.h"

namespace fm
{

class mem_dense_matrix;
class mem_col_dense_matrix;
class bulk_operate;

class submatrix_compute
{
	size_t start_row;
	size_t start_col;
public:
	typedef std::shared_ptr<submatrix_compute> ptr;
	submatrix_compute(size_t start_row, size_t start_col) {
		this->start_row = start_row;
		this->start_col = start_col;
	}

	virtual void run(const mem_dense_matrix &subm) = 0;
};

class EM_dense_matrix_accessor
{
public:
	typedef std::shared_ptr<EM_dense_matrix_accessor> ptr;

	virtual ~EM_dense_matrix_accessor() {
	}

	virtual bool fetch_submatrix(size_t start_row, size_t nrow,
			size_t start_col, size_t ncol,
			submatrix_compute::ptr compute) const = 0;
	virtual bool set_submatrix(size_t start_row, size_t start_col,
			std::shared_ptr<mem_dense_matrix> subm) = 0;
	virtual void wait4complete(int num) = 0;
	virtual void wait4all() = 0;
};

class EM_dense_matrix
{
public:
	typedef std::shared_ptr<EM_dense_matrix> ptr;

	virtual ~EM_dense_matrix() {
	}

	virtual EM_dense_matrix::ptr inner_prod(const mem_dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op) = 0;

	virtual void resize(size_t nrow, size_t ncol) = 0;
	virtual size_t get_num_rows() const = 0;
	virtual size_t get_num_cols() const = 0;
	virtual size_t get_entry_size() const = 0;
	virtual EM_dense_matrix_accessor::ptr create_accessor() = 0;
};

/*
 * In this matrix class, data is stored in columns.
 * This has to be a very narrow matrix, i.e., the column length must be much
 * larger than the row length.
 */
class EM_col_dense_matrix: public EM_dense_matrix
{
	// The number of elements.
	static const size_t COL_CHUNK_SIZE = 1024 * 1024;

	size_t entry_size;
	std::vector<EM_vector::ptr> cols;

	EM_col_dense_matrix(size_t entry_size) {
		this->entry_size = entry_size;
	}

	EM_col_dense_matrix(size_t nrow, size_t ncol, size_t entry_size) {
		this->entry_size = entry_size;
		this->resize(nrow, ncol);
	}
public:
	static ptr create(size_t entry_size) {
		return ptr(new EM_col_dense_matrix(entry_size));
	}

	static ptr create(size_t nrow, size_t ncol, size_t entry_size) {
		return ptr(new EM_col_dense_matrix(nrow, ncol, entry_size));
	}

	virtual void resize(size_t nrow, size_t ncol);

	virtual EM_dense_matrix::ptr inner_prod(const mem_dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op);

	virtual size_t get_num_rows() const {
		if (cols.empty())
			return 0;
		else
			return cols[0]->get_size();
	}

	virtual size_t get_num_cols() const {
		return cols.size();
	}

	virtual size_t get_entry_size() const {
		return entry_size;
	}

	EM_dense_matrix_accessor::ptr create_accessor();
};

class EM_col_matrix_accessor: public EM_dense_matrix_accessor
{
	EM_col_dense_matrix &m;
	std::vector<EM_vector_accessor::ptr> accessors;
public:
	typedef std::shared_ptr<EM_col_matrix_accessor> ptr;

	EM_col_matrix_accessor(EM_col_dense_matrix &_m,
			const std::vector<EM_vector::ptr> &cols): m(_m) {
		accessors.resize(cols.size());
		for (size_t i = 0; i < cols.size(); i++)
			accessors[i] = cols[i]->create_accessor();
	}

	bool fetch_submatrix(size_t start_row, size_t nrow,
			size_t start_col, size_t ncol, submatrix_compute::ptr compute) const;
	bool set_submatrix(size_t start_row, size_t start_col,
			std::shared_ptr<mem_dense_matrix> subm);
	virtual void wait4complete(int num);
	virtual void wait4all();
};

}

#endif
