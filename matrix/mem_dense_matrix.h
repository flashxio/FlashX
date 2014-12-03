#ifndef __MEM_DENSE_MATRIX_H__
#define __MEM_DENSE_MATRIX_H__

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

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <malloc.h>

#include "common.h"
#include "bulk_operate.h"

namespace fm
{

class bulk_operate;

enum matrix_layout_t
{
	L_COL,
	L_ROW,
};

class dense_matrix
{
public:
	virtual ~dense_matrix() {
	}
	virtual size_t get_entry_size() const = 0;
	virtual size_t get_num_rows() const = 0;
	virtual size_t get_num_cols() const = 0;
};

class mem_dense_matrix: public dense_matrix
{
	size_t nrow;
	size_t ncol;
	size_t entry_size;
protected:
	mem_dense_matrix(size_t nrow, size_t ncol, size_t entry_size) {
		this->nrow = nrow;
		this->ncol = ncol;
		this->entry_size = entry_size;
	}
public:
	typedef std::shared_ptr<mem_dense_matrix> ptr;

	virtual const char *get(size_t row, size_t col) const = 0;
	virtual matrix_layout_t store_layout() const = 0;

	size_t get_entry_size() const {
		return entry_size;
	}

	virtual size_t get_num_rows() const {
		return nrow;
	}

	virtual size_t get_num_cols() const {
		return ncol;
	}

	virtual mem_dense_matrix::ptr inner_prod(const mem_dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op) const = 0;
};

class mem_row_dense_matrix: public mem_dense_matrix
{
	char *data;

	mem_row_dense_matrix(size_t nrow, size_t ncol,
			size_t entry_size): mem_dense_matrix(nrow, ncol, entry_size) {
		if (nrow * ncol > 0) {
			data = (char *) memalign(PAGE_SIZE, nrow * ncol * entry_size);
			assert(data);
		}
	}
public:
	typedef std::shared_ptr<mem_row_dense_matrix> ptr;

	static ptr create(size_t nrow, size_t ncol, size_t entry_size) {
		return ptr(new mem_row_dense_matrix(nrow, ncol, entry_size));
	}

	~mem_row_dense_matrix() {
		free(data);
	}

	mem_dense_matrix::ptr inner_prod(const mem_dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op) const;

	char *get_row(size_t row) {
		return data + row * get_num_cols() * get_entry_size();
	}

	const char *get_row(size_t row) const {
		return data + row * get_num_cols() * get_entry_size();
	}

	char *get(size_t row, size_t col) {
		return get_row(row) + col * get_entry_size();
	}

	const char *get(size_t row, size_t col) const {
		return get_row(row) + col * get_entry_size();
	}

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_ROW;
	}
};

class mem_col_dense_matrix: public mem_dense_matrix
{
	// The number of rows.
	static const size_t SUB_CHUNK_SIZE = 1024;

	class sub_matrix {
		size_t start_row;
		size_t start_col;
		size_t nrow;
		size_t ncol;
		const mem_col_dense_matrix &m;
	public:
		sub_matrix(size_t start_row, size_t nrow, size_t start_col,
				size_t ncol, const mem_col_dense_matrix &_m): m(_m) {
			this->start_row = start_row;
			this->start_col = start_col;
			this->nrow = nrow;
			this->ncol = ncol;
			assert(start_row + nrow <= m.get_num_rows());
			assert(start_col + ncol <= m.get_num_cols());
		}

		size_t get_num_rows() const {
			return nrow;
		}

		size_t get_num_cols() const {
			return ncol;
		}

		const char *get_col(size_t col) const {
			return m.get_col(start_col + col) + start_row * m.get_entry_size();
		}
	};

	char *data;

	mem_col_dense_matrix(size_t nrow, size_t ncol,
			size_t entry_size): mem_dense_matrix(nrow, ncol, entry_size) {
		if (nrow * ncol > 0) {
			data = (char *) memalign(PAGE_SIZE, nrow * ncol * entry_size);
			assert(data);
		}
	}
public:
	typedef std::shared_ptr<mem_col_dense_matrix> ptr;

	static ptr create(size_t nrow, size_t ncol, size_t entry_size) {
		return ptr(new mem_col_dense_matrix(nrow, ncol, entry_size));
	}

	~mem_col_dense_matrix() {
		free(data);
	}

	mem_dense_matrix::ptr inner_prod(const mem_dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op) const;

	void set_col(char *buf, size_t size, size_t col) {
		assert(size == get_entry_size() * get_num_rows());
		memcpy(get_col(col), buf, size);
	}

	char *get_col(size_t col) {
		return data + col * get_num_rows() * get_entry_size();
	}

	const char *get_col(size_t col) const {
		return data + col * get_num_rows() * get_entry_size();
	}

	char *get(size_t row, size_t col) {
		return get_col(col) + row * get_entry_size();
	}

	const char *get(size_t row, size_t col) const {
		return get_col(col) + row * get_entry_size();
	}

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_COL;
	}
};

template<class EntryType>
class type_mem_dense_matrix
{
	mem_dense_matrix::ptr m;

	type_mem_dense_matrix(size_t nrow, size_t ncol, matrix_layout_t layout) {
		if (layout == matrix_layout_t::L_COL)
			m = mem_col_dense_matrix::create(nrow, ncol, sizeof(EntryType));
		else if (layout == matrix_layout_t::L_ROW)
			m = mem_row_dense_matrix::create(nrow, ncol, sizeof(EntryType));
		else
			assert(0);
	}

	type_mem_dense_matrix(mem_dense_matrix::ptr m) {
		this->m = m;
	}
public:
	typedef std::shared_ptr<type_mem_dense_matrix<EntryType> > ptr;

	static ptr create(size_t nrow, size_t ncol, matrix_layout_t layout) {
		return ptr(new type_mem_dense_matrix<EntryType>(nrow, ncol, layout));
	}

	static ptr create(mem_dense_matrix::ptr m) {
		assert(m->get_entry_size() == sizeof(EntryType));
		return ptr(new type_mem_dense_matrix<EntryType>(m));
	}

	size_t get_num_rows() const {
		return m->get_num_rows();
	}

	size_t get_num_cols() const {
		return m->get_num_cols();
	}

	void set(size_t row, size_t col, const EntryType &v) {
		*(EntryType *) m->get(row, col) = v;
	}

	EntryType get(size_t row, size_t col) const {
		return *(EntryType *) m->get(row, col);
	}

	const mem_dense_matrix &get_matrix() const {
		return *m;
	}
};

typedef type_mem_dense_matrix<int> I_mem_dense_matrix;
typedef type_mem_dense_matrix<double> D_mem_dense_matrix;

template<class LeftType, class RightType, class ResType>
mem_dense_matrix::ptr  multiply(mem_col_dense_matrix &m1, mem_dense_matrix &m2)
{
	basic_ops_impl<LeftType, RightType, ResType> ops;
	return m1.inner_prod(m2, ops.get_multiply(), ops.get_add());
}

template<class LeftType, class RightType, class ResType>
typename type_mem_dense_matrix<ResType>::ptr  multiply(
		type_mem_dense_matrix<LeftType> &m1,
		type_mem_dense_matrix<RightType> &m2)
{
	basic_ops_impl<LeftType, RightType, ResType> ops;
	return type_mem_dense_matrix<ResType>::create(
			m1.get_matrix().inner_prod(m2.get_matrix(), ops.get_multiply(),
			ops.get_add()));
}

}

#endif
