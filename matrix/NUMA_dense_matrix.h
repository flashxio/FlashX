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

#include "mem_dense_matrix.h"
#include "raw_data_array.h"
#include "NUMA_mapper.h"
#include "NUMA_vector.h"

namespace fm
{

/*
 * This class defines an in-memory dense matrix whose data is organized
 * in row major and has many more rows than columns. The matrix is optimized
 * for a NUMA machine. Rows in the matrix are distributed across multiple
 * NUMA node. Multiple adjacent rows are stored in contiguous memory and
 * a single row is guaranteed to be stored together.
 */
class NUMA_row_tall_dense_matrix: public dense_matrix
{
	// This is to map rows to different NUMA nodes.
	detail::NUMA_mapper mapper;
	std::vector<detail::raw_data_array> data;

	// The copy constructor performs shallow copy.
	NUMA_row_tall_dense_matrix(const NUMA_row_tall_dense_matrix &mat): dense_matrix(
			mat.get_num_rows(), mat.get_num_cols(), mat.get_type(),
			true), mapper(mat.get_num_nodes()) {
		this->data = mat.data;
	}

	NUMA_row_tall_dense_matrix(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type);
public:
	typedef std::shared_ptr<NUMA_row_tall_dense_matrix> ptr;

	static ptr create(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type) {
		return ptr(new NUMA_row_tall_dense_matrix(nrow, ncol, num_nodes, type));
	}

	static ptr cast(dense_matrix::ptr mat) {
		if (!mat->is_in_mem()) {
			BOOST_LOG_TRIVIAL(error) << "the matrix isn't in memory";
			return ptr();
		}
		if (mat->store_layout() != matrix_layout_t::L_ROW) {
			BOOST_LOG_TRIVIAL(error) << "the matrix isn't row major";
			return ptr();
		}
		// TODO how do we make sure it's a NUMA matrix?
		return std::static_pointer_cast<NUMA_row_tall_dense_matrix>(mat);
	}

	size_t get_num_nodes() const {
		return data.size();
	}

	const char *get_row(off_t row_idx) const;
	char *get_row(off_t row_idx);
	const char *get_rows(off_t row_start, off_t row_end) const;
	char *get_rows(off_t row_start, off_t row_end);

	template<class T>
	void set(size_t row_idx, size_t col_idx, T val) {
		char *row = get_row(row_idx);
		*(T *) (row + col_idx * get_entry_size()) = val;
	}

	template<class T>
	T get(size_t row_idx, size_t col_idx) const {
		const char *row = get_row(row_idx);
		return *(const T *) (row + col_idx * get_entry_size());
	}

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_ROW;
	}

	virtual void reset_data() {
		reset_arrays(data);
	}

	virtual void set_data(const set_operate &op);

	virtual bool copy_from(const dense_matrix &mat) {
		// TODO
		assert(0);
	}

	virtual bool write2file(const std::string &file_name) const {
		// TODO
		assert(0);
	}

	virtual dense_matrix::ptr shallow_copy() const {
		// TODO
		assert(0);
	}

	virtual dense_matrix::ptr deep_copy() const;

	virtual dense_matrix::ptr conv2(size_t nrow, size_t ncol, bool byrow) const {
		BOOST_LOG_TRIVIAL(error) << "NUMA matrix doesn't support conv2";
		return dense_matrix::ptr();
	}

	virtual dense_matrix::ptr transpose() const {
		// TODO
		assert(0);
	}

	virtual dense_matrix::ptr inner_prod(const dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op) const {
		// TODO
		assert(0);
	}

	virtual scalar_variable::ptr aggregate(const bulk_operate &op) const {
		// TODO
		assert(0);
	}

	virtual dense_matrix::ptr mapply2(const dense_matrix &m,
			const bulk_operate &op) const {
		// TODO
		assert(0);
	}

	virtual dense_matrix::ptr sapply(const bulk_uoperate &op) const {
		// TODO
		assert(0);
	}

	virtual dense_matrix::ptr apply(apply_margin margin,
			const arr_apply_operate &op) const {
		// TODO
		assert(0);
	}
};

/*
 * This class defines an in-memory dense matrix whose data is organized
 * in column major and has many more rows than columns. The matrix is optimized
 * for a NUMA machine. Each column in the matrix are distributed across multiple
 * NUMA nodes.
 */
class NUMA_col_tall_dense_matrix: public dense_matrix
{
	std::vector<NUMA_vector::ptr> data;

	NUMA_col_tall_dense_matrix(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type);
public:
	typedef std::shared_ptr<NUMA_col_tall_dense_matrix> ptr;

	static ptr create(size_t nrow, size_t ncol, int num_nodes,
			const scalar_type &type) {
		return ptr(new NUMA_col_tall_dense_matrix(nrow, ncol, num_nodes, type));
	}

	static ptr cast(dense_matrix::ptr mat) {
		if (!mat->is_in_mem()) {
			BOOST_LOG_TRIVIAL(error) << "the matrix isn't in memory";
			return ptr();
		}
		if (mat->store_layout() != matrix_layout_t::L_COL) {
			BOOST_LOG_TRIVIAL(error) << "the matrix isn't col major";
			return ptr();
		}
		// TODO how do we make sure it's a NUMA matrix?
		return std::static_pointer_cast<NUMA_col_tall_dense_matrix>(mat);
	}

	size_t get_num_nodes() const {
		return data[0]->get_num_nodes();
	}

	template<class T>
	void set(size_t row_idx, size_t col_idx, T val) {
		data[col_idx]->set<T>(row_idx, val);
	}

	template<class T>
	T get(size_t row_idx, size_t col_idx) const {
		data[col_idx]->get<T>(row_idx);
	}

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_COL;
	}

	virtual void reset_data() {
		for (size_t i = 0; i < data.size(); i++)
			data[i]->reset_data();
	}

	virtual void set_data(const set_operate &op);

	virtual bool copy_from(const dense_matrix &mat) {
		// TODO
		assert(0);
	}

	virtual bool write2file(const std::string &file_name) const {
		// TODO
		assert(0);
	}

	virtual dense_matrix::ptr shallow_copy() const {
		// TODO
		assert(0);
	}

	virtual dense_matrix::ptr deep_copy() const;

	virtual dense_matrix::ptr conv2(size_t nrow, size_t ncol, bool byrow) const {
		BOOST_LOG_TRIVIAL(error) << "NUMA matrix doesn't support conv2";
		return dense_matrix::ptr();
	}

	virtual dense_matrix::ptr transpose() const {
		// TODO
		assert(0);
	}

	virtual dense_matrix::ptr inner_prod(const dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op) const;

	virtual scalar_variable::ptr aggregate(const bulk_operate &op) const {
		// TODO
		assert(0);
	}

	virtual dense_matrix::ptr mapply2(const dense_matrix &m,
			const bulk_operate &op) const;

	virtual dense_matrix::ptr sapply(const bulk_uoperate &op) const {
		// TODO
		assert(0);
	}

	virtual dense_matrix::ptr apply(apply_margin margin,
			const arr_apply_operate &op) const {
		// TODO
		assert(0);
	}

};

}

#endif
