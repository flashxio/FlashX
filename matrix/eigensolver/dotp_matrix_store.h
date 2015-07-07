#ifndef __NORM2_MATRIX_STORE_H__
#define __NORM2_MATRIX_STORE_H__

/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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
#include "EM_object.h"

namespace fm
{

namespace detail
{

class sub_dotp_matrix_store: public virtual_matrix_store
{
	matrix_store::const_ptr orig_store;
	std::vector<off_t> idxs;
	std::vector<double> col_dot_prods;
public:
	sub_dotp_matrix_store(matrix_store::const_ptr orig_store,
			const std::vector<double> &dot_prods,
			const std::vector<off_t> &idxs): virtual_matrix_store(
				orig_store->get_num_rows(), idxs.size(),
				orig_store->is_in_mem(), orig_store->get_type()) {
		this->orig_store = orig_store;
		this->col_dot_prods = dot_prods;
		this->idxs = idxs;
		assert(dot_prods.size() == idxs.size());
		assert(orig_store->get_type() == get_scalar_type<double>());
	}

	matrix_store::const_ptr get_orig_store() const {
		return orig_store;
	}

	double get_col_dot(off_t idx) const {
		return col_dot_prods[idx];
	}

	std::string get_name() const {
		return (boost::format("sub_dotp_mat(%1%,%2%)")
				% orig_store->get_num_rows() % idxs.size()).str();
	}

	matrix_layout_t store_layout() const {
		return orig_store->store_layout();
	}

	matrix_store::const_ptr transpose() const {
		return matrix_store::const_ptr();
	}

	virtual local_matrix_store::const_ptr get_portion_async(size_t start_row,
			size_t start_col, size_t num_rows, size_t num_cols,
			portion_compute::ptr compute) const {
		return local_matrix_store::const_ptr();
	}
	virtual local_matrix_store::const_ptr get_portion(size_t start_row,
			size_t start_col, size_t num_rows, size_t num_cols) const {
		return local_matrix_store::const_ptr();
	}
	virtual std::pair<size_t, size_t> get_portion_size() const {
		return std::pair<size_t, size_t>(0, 0);
	}

	virtual matrix_store::ptr materialize() const {
		return matrix_store::ptr();
	}

};

/*
 * This matrix store helps maintain the dot products of each column
 * of the matrix. When we need to fetch a column, this class returns a special
 * virtual matrix that has the dot product of that column.
 */
class dotp_matrix_store: public virtual_matrix_store
{
	matrix_store::const_ptr orig_store;
	std::vector<double> col_dot_prods;

	dotp_matrix_store(matrix_store::const_ptr orig_store): virtual_matrix_store(
				orig_store->get_num_rows(), orig_store->get_num_cols(),
				orig_store->is_in_mem(), orig_store->get_type()) {
		this->orig_store = orig_store;

		dense_matrix::ptr mat = dense_matrix::create(orig_store);
		const bulk_uoperate *op = get_type().get_basic_uops().get_op(
				basic_uops::op_idx::SQ);
		dense_matrix::ptr sq_mat = mat->sapply(bulk_uoperate::conv2ptr(*op));
		vector::ptr sums = sq_mat->col_sum();

		const fm::detail::smp_vec_store &smp_res
			= dynamic_cast<const fm::detail::smp_vec_store &>(sums->get_data());
		col_dot_prods.resize(sums->get_length());
		for (size_t i = 0; i < sums->get_length(); i++)
			col_dot_prods[i] = smp_res.get<double>(i);
	}
public:
	typedef std::shared_ptr<dotp_matrix_store> ptr;

	static ptr create(matrix_store::const_ptr store) {
		return ptr(new dotp_matrix_store(store));
	}

	const std::vector<double> get_col_dot_prods() const {
		return col_dot_prods;
	}

	virtual matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const {
		std::vector<double> sub_dot_prods(idxs.size());
		for (size_t i = 0; i < idxs.size(); i++)
			sub_dot_prods[i] = col_dot_prods[idxs[i]];
		return matrix_store::const_ptr(new sub_dotp_matrix_store(orig_store,
					sub_dot_prods, idxs));
	}

	std::string get_name() const {
		return orig_store->get_name();
	}

	matrix_layout_t store_layout() const {
		return orig_store->store_layout();
	}

	matrix_store::const_ptr transpose() const {
		return orig_store->transpose();
	}

	virtual local_matrix_store::const_ptr get_portion_async(size_t start_row,
			size_t start_col, size_t num_rows, size_t num_cols,
			portion_compute::ptr compute) const {
		return orig_store->get_portion_async(start_row, start_col, num_rows,
				num_cols, compute);
	}
	virtual local_matrix_store::const_ptr get_portion(size_t start_row,
			size_t start_col, size_t num_rows, size_t num_cols) const {
		return orig_store->get_portion(start_row, start_col, num_rows, num_cols);
	}
	virtual std::pair<size_t, size_t> get_portion_size() const {
		return orig_store->get_portion_size();
	}

	virtual matrix_store::ptr materialize() const {
		return matrix_store::ptr(new dotp_matrix_store(*this));
	}

	virtual int get_num_nodes() const {
		return orig_store->get_num_nodes();
	}
};

}

}

#endif
