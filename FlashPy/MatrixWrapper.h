#ifndef __FLASHPY_MATRIX_WRAPPER_H__
#define __FLASHPY_MATRIX_WRAPPER_H__

/*
 * Copyright 2017 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashPy.
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

#include <numpy/ndarraytypes.h>

#include "dense_matrix.h"
#include "col_vec.h"
#include "factor.h"
#include "mem_matrix_store.h"

namespace flashpy
{

typedef int bulk_op_idx_t;
typedef int bulk_uop_idx_t;

static inline fm::bulk_operate::const_ptr get_op(const fm::scalar_type &type,
		bulk_op_idx_t idx)
{
	if (idx >= fm::basic_ops::op_idx::NUM_OPS)
		throw std::invalid_argument("invalid operation");
	auto op = type.get_basic_ops().get_op((fm::basic_ops::op_idx) idx);
	return fm::bulk_operate::conv2ptr(*op);
}

static inline fm::bulk_uoperate::const_ptr get_uop(const fm::scalar_type &type,
		bulk_uop_idx_t idx)
{
	if (idx >= fm::basic_uops::op_idx::NUM_OPS)
		throw std::invalid_argument("invalid operation");
	auto op = type.get_basic_uops().get_op((fm::basic_uops::op_idx) idx);
	return fm::bulk_uoperate::conv2ptr(*op);
}

static inline fm::agg_operate::const_ptr get_agg(const fm::scalar_type &type,
		bulk_op_idx_t idx)
{
	return fm::agg_operate::create(get_op(type, idx));
}

static inline std::shared_ptr<fm::col_vec> get_vec(fm::dense_matrix::ptr mat)
{
	auto ret = std::dynamic_pointer_cast<fm::col_vec>(mat);
	if (ret == NULL)
		throw std::invalid_argument("invalid matrix, wants a col vector");
	else
		return ret;
}

static inline std::shared_ptr<fm::factor_col_vector> get_factor(
		fm::dense_matrix::ptr mat)
{
	auto ret = std::dynamic_pointer_cast<fm::factor_col_vector>(mat);
	if (ret == NULL)
		throw std::invalid_argument("invalid matrix, want a factor col vector");
	else
		return ret;
}

const fm::scalar_type &convT_py2fm(const std::string &t);

/*
 * This is a wrapper for dense_matrix in FlashMatrix.
 * The reason that we create this class is that we want to minimize the number
 * of classes we expose to Cython, so that we only need to declare one C++
 * class.
 */
class matrix_wrapper
{
	fm::dense_matrix::ptr mat;

	void check_mat() const {
		if (!is_valid())
			throw std::invalid_argument("invalid matrix");
	}

	bool is_valid() const {
		return mat != NULL;
	}

	static fm::matrix_layout_t get_layout(const std::string &layout) {
		if (layout == "c")
			return fm::matrix_layout_t::L_COL;
		else if (layout == "r")
			return fm::matrix_layout_t::L_ROW;
		else
			throw std::invalid_argument("wrong layout");
	}

public:
	matrix_wrapper() {
	}

	matrix_wrapper(intptr_t data_addr, size_t length, const std::string &t);

	matrix_wrapper(intptr_t data_addr, size_t nrow, size_t ncol,
			const std::string &t, const std::string layout);

	matrix_wrapper(size_t nrow, size_t ncol, const std::string &t,
			const std::string layout) {
		mat = fm::dense_matrix::create(nrow, ncol, get_layout(layout),
				convT_py2fm(t));
	}

	matrix_wrapper(fm::dense_matrix::ptr mat) {
		this->mat = mat;
	}

	bool is_vector() const;

	std::string get_type_str() const;
	/*
	 * Get a Python type.
	 */
	enum NPY_TYPES get_type_py() const;

	template<class T>
	void init_seq(T start, T stride, size_t nrow, size_t ncol,
			std::string layout, bool byrow, int num_nodes, bool in_mem) {
		fm::matrix_layout_t layout_val = layout
			== "c" ? fm::matrix_layout_t::L_COL : fm::matrix_layout_t::L_ROW;
		mat = fm::dense_matrix::create_seq<T>(start, stride, nrow, ncol,
				layout_val, byrow, num_nodes, in_mem);
	}

	size_t get_num_rows() const {
		check_mat();
		return mat->get_num_rows();
	}

	size_t get_num_cols() const {
		check_mat();
		return mat->get_num_cols();
	}

	size_t get_entry_size() const {
		check_mat();
		return mat->get_type().get_size();
	}

	const char *get_raw_arr() const {
		check_mat();
		mat->move_store(true, -1);
		auto store = std::dynamic_pointer_cast<const fm::detail::mem_matrix_store>(
				mat->get_raw_store());
		return store->get_raw_arr();
	}

	char *get_raw_arr() {
		check_mat();
		mat->move_store(true, -1);
		auto store = std::dynamic_pointer_cast<const fm::detail::mem_matrix_store>(
				mat->get_raw_store());
		return (char *) store->get_raw_arr();
	}

	bool is_in_mem() const {
		check_mat();
		return mat->is_in_mem();
	}

	bool is_virtual() const {
		check_mat();
		return mat->is_virtual();
	}

	bool materialize_self() const {
		check_mat();
		return mat->materialize_self();
	}

	matrix_wrapper get_cols(const std::vector<off_t> &idxs) const {
		check_mat();
		return matrix_wrapper(mat->get_cols(idxs));
	}

	matrix_wrapper get_rows(const std::vector<off_t> &idxs) const {
		check_mat();
		return matrix_wrapper(mat->get_rows(idxs));
	}

	matrix_wrapper get_cols(matrix_wrapper idxs) const {
		check_mat();
		idxs.check_mat();
		return matrix_wrapper(mat->get_cols(get_vec(idxs.mat)));
	}

	matrix_wrapper get_rows(matrix_wrapper idxs) const {
		check_mat();
		idxs.check_mat();
		return matrix_wrapper(mat->get_rows(get_vec(idxs.mat)));
	}

	matrix_wrapper get_cols(size_t start, size_t end) const {
		check_mat();
		return matrix_wrapper(mat->get_cols(start, end));
	}

	matrix_wrapper get_rows(size_t start, size_t end) const {
		check_mat();
		return matrix_wrapper(mat->get_rows(start, end));
	}

	matrix_wrapper set_cols(const std::vector<off_t> &idxs,
			matrix_wrapper cols) {
		check_mat();
		return matrix_wrapper(mat->set_cols(idxs, cols.mat));
	}
	matrix_wrapper set_rows(const std::vector<off_t> &idxs,
			matrix_wrapper rows) {
		check_mat();
		return matrix_wrapper(mat->set_rows(idxs, rows.mat));
	}

	matrix_wrapper transpose() const {
		check_mat();
		return matrix_wrapper(mat->transpose());
	}

	matrix_wrapper conv_store(bool in_mem, int num_nodes) const {
		check_mat();
		return matrix_wrapper(mat->conv_store(in_mem, num_nodes));
	}

	matrix_wrapper inner_prod(matrix_wrapper m, bulk_op_idx_t left_op,
			bulk_op_idx_t right_op) const {
		check_mat();
		m.check_mat();
		return matrix_wrapper(mat->inner_prod(*m.mat,
					get_op(mat->get_type(), left_op),
					get_op(mat->get_type(), right_op)));
	}
	matrix_wrapper multiply(matrix_wrapper m) const {
		check_mat();
		m.check_mat();
		return matrix_wrapper(mat->multiply(*m.mat));
	}

	matrix_wrapper agg_row(bulk_op_idx_t op) const {
		check_mat();
		return matrix_wrapper(mat->aggregate(fm::matrix_margin::MAR_ROW,
					get_agg(mat->get_type(), op)));
	}
	matrix_wrapper agg_col(bulk_op_idx_t op) const {
		check_mat();
		return matrix_wrapper(mat->aggregate(fm::matrix_margin::MAR_COL,
					get_agg(mat->get_type(), op)));
	}

	matrix_wrapper groupby_row(matrix_wrapper labels, bulk_op_idx_t op) const {
		check_mat();
		labels.check_mat();
		return matrix_wrapper(mat->groupby_row(get_factor(labels.mat),
					get_agg(mat->get_type(), op)));
	}

	matrix_wrapper mapply_cols(matrix_wrapper vals, bulk_op_idx_t op) const {
		check_mat();
		vals.check_mat();
		return matrix_wrapper(mat->mapply_cols(get_vec(vals.mat),
					get_op(mat->get_type(), op)));
	}

	matrix_wrapper mapply_rows(matrix_wrapper vals, bulk_op_idx_t op) const {
		check_mat();
		vals.check_mat();
		return matrix_wrapper(mat->mapply_rows(get_vec(vals.mat),
					get_op(mat->get_type(), op)));
	}

	matrix_wrapper mapply2(matrix_wrapper m, bulk_op_idx_t op) const {
		check_mat();
		m.check_mat();
		auto res = mat->mapply2(*m.mat, get_op(mat->get_type(), op));
		if (m.is_vector() && is_vector())
			return matrix_wrapper(fm::col_vec::create(res));
		else
			return matrix_wrapper(res);
	}

	matrix_wrapper sapply(bulk_uop_idx_t op) const {
		check_mat();
		auto res = mat->sapply(get_uop(mat->get_type(), op));
		if (is_vector())
			return matrix_wrapper(fm::col_vec::create(res));
		else
			return matrix_wrapper(res);
	}
};

}

#endif
