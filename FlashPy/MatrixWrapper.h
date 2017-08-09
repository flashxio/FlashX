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
typedef int agg_op_idx_t;

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
		agg_op_idx_t idx)
{
	if (idx >= fm::agg_ops::op_idx::NUM_OPS)
		throw std::invalid_argument("invalid agg operation");
	return type.get_agg_ops().get_op((fm::agg_ops::op_idx) idx);
}

static inline std::shared_ptr<fm::col_vec> get_vec(fm::dense_matrix::ptr mat)
{
	auto ret = std::dynamic_pointer_cast<fm::col_vec>(mat);
	if (ret)
		return ret;
	ret = fm::col_vec::create(mat);
	return ret;
}

const fm::scalar_type &convT_py2fm(const std::string &t);

class scalar_wrapper
{
	fm::scalar_variable::ptr var;
public:
	scalar_wrapper() {
	}

	scalar_wrapper(fm::scalar_variable::ptr var) {
		this->var = var;
	}

	const fm::scalar_type &get_type() {
		return var->get_type();
	}

	const char *get_raw() const {
		return var->get_raw();
	}

	fm::scalar_variable::ptr get_var(const fm::scalar_type *type = NULL) const;
};

template<class T>
scalar_wrapper create_scalar_wrapper(T val)
{
	fm::scalar_variable::ptr var(new fm::scalar_variable_impl<T>(val));
	return scalar_wrapper(var);
}

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

	static fm::matrix_layout_t get_layout(const std::string &layout) {
		// The layout used by FORTRAN
		if (layout == "F")
			return fm::matrix_layout_t::L_COL;
		else if (layout == "C")
			return fm::matrix_layout_t::L_ROW;
		else
			throw std::invalid_argument("wrong layout");
	}

public:
	static matrix_wrapper cbind(const std::vector<matrix_wrapper> &mats);
	static matrix_wrapper rbind(const std::vector<matrix_wrapper> &mats);
	matrix_wrapper() {
	}

	matrix_wrapper(intptr_t data_addr, size_t length, const std::string &t);
	matrix_wrapper(intptr_t data_addr, size_t nrow, size_t ncol,
			const std::string &t, const std::string layout);
	matrix_wrapper(size_t length, std::string &t);
	matrix_wrapper(size_t nrow, size_t ncol, const std::string &t,
			const std::string layout);

	matrix_wrapper(fm::dense_matrix::ptr mat) {
		this->mat = mat;
	}

	bool is_valid() const {
		return mat != NULL;
	}

	bool is_vector() const;

	std::string get_type_str() const;
	/*
	 * Get a Python type.
	 */
	enum NPY_TYPES get_type_py() const;

	template<class T>
	void init_seq(T start, T stride, bool byrow) {
		check_mat();
		bool is_vec = is_vector();
		mat = fm::dense_matrix::create_seq<T>(start, stride,
				mat->get_num_rows(), mat->get_num_cols(), mat->store_layout(),
				byrow, mat->get_raw_store()->get_num_nodes(), mat->is_in_mem());
		if (is_vec)
			mat = fm::col_vec::create(mat);
	}

	void init_const_float(double val);
	void init_const_int(long val);

	void set_cached(bool);

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

	std::string get_layout() const {
		if (mat->store_layout() == fm::matrix_layout_t::L_COL)
			return "F";
		else
			return "C";
	}

	bool copy_rows_to(char *arr, size_t len) const;

	bool is_floating_point() const {
		return mat->get_type().is_floating_point();
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

	matrix_wrapper as_factor(int num_levels) const;
	matrix_wrapper as_vector() const;
	matrix_wrapper as_matrix() const;

	matrix_wrapper cast_ele_type(std::string dtyp) const;

	matrix_wrapper get_col(long idx) const;
	matrix_wrapper get_row(long idx) const;

	matrix_wrapper get_cols(const std::vector<off_t> &idxs) const;
	matrix_wrapper get_rows(const std::vector<off_t> &idxs) const;

	matrix_wrapper get_cols(matrix_wrapper idxs) const;
	matrix_wrapper get_rows(matrix_wrapper idxs) const;

	matrix_wrapper get_cols(size_t start, size_t end, long step) const;
	matrix_wrapper get_rows(size_t start, size_t end, long step) const;

	matrix_wrapper get_eles(matrix_wrapper idxs) const;

	matrix_wrapper set_cols(const std::vector<off_t> &idxs,
			matrix_wrapper cols);
	matrix_wrapper set_rows(const std::vector<off_t> &idxs,
			matrix_wrapper rows);

	matrix_wrapper set_cols(long start, long stop, long step, matrix_wrapper cols) {
		check_mat();
		return matrix_wrapper(mat->set_cols(start, stop, step, cols.mat));
	}

	matrix_wrapper set_rows(long start, long stop, long step, matrix_wrapper rows) {
		check_mat();
		return matrix_wrapper(mat->set_rows(start, stop, step, rows.mat));
	}

	matrix_wrapper transpose() const {
		check_mat();
		return matrix_wrapper(mat->transpose());
	}

	matrix_wrapper conv_layout(const std::string &layout) const {
		check_mat();
		return matrix_wrapper(mat->conv2(get_layout(layout)));
	}

	matrix_wrapper conv_store(bool in_mem, int num_nodes) const {
		check_mat();
		return matrix_wrapper(mat->conv_store(in_mem, num_nodes));
	}

	matrix_wrapper inner_prod(matrix_wrapper m, bulk_op_idx_t left_op,
			bulk_op_idx_t right_op) const;
	matrix_wrapper multiply(matrix_wrapper m) const;

	matrix_wrapper aggregate(agg_op_idx_t op) const {
		check_mat();
		return matrix_wrapper(mat->aggregate(fm::matrix_margin::BOTH,
					get_agg(mat->get_type(), op)));
	}

	matrix_wrapper agg_row(agg_op_idx_t op) const {
		check_mat();
		return matrix_wrapper(mat->aggregate(fm::matrix_margin::MAR_ROW,
					get_agg(mat->get_type(), op)));
	}
	matrix_wrapper agg_col(agg_op_idx_t op) const {
		check_mat();
		return matrix_wrapper(mat->aggregate(fm::matrix_margin::MAR_COL,
					get_agg(mat->get_type(), op)));
	}

	matrix_wrapper cum_row(agg_op_idx_t op) const;
	matrix_wrapper cum_col(agg_op_idx_t op) const;

	matrix_wrapper groupby_row(matrix_wrapper labels, agg_op_idx_t op) const;
	// TODO we should return data_frame later.
	std::pair<matrix_wrapper, matrix_wrapper> groupby(agg_op_idx_t op,
			bool with_val) const;

	matrix_wrapper mapply_cols(matrix_wrapper vals, bulk_op_idx_t op) const;
	matrix_wrapper mapply_rows(matrix_wrapper vals, bulk_op_idx_t op) const;
	matrix_wrapper mapply2(matrix_wrapper m, bulk_op_idx_t op) const;

	matrix_wrapper sapply(bulk_uop_idx_t op) const {
		check_mat();
		auto res = mat->sapply(get_uop(mat->get_type(), op));
		if (is_vector())
			return matrix_wrapper(fm::col_vec::create(res));
		else
			return matrix_wrapper(res);
	}

	matrix_wrapper apply_scalar(scalar_wrapper var, bulk_op_idx_t op) const;

	matrix_wrapper ifelse(matrix_wrapper x, matrix_wrapper y) const;
};

bool init_flashpy_c(const std::string &file);

}

#endif
