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

#include "collected_col_matrix_store.h"
#include "dense_matrix.h"

namespace fm
{

namespace eigen
{

namespace
{

std::string cat_mat_names(
		const std::vector<detail::matrix_store::const_ptr> &mats)
{
	std::string ret;
	for (size_t i = 0; i < mats.size() - 1; i++)
		ret += mats[i]->get_name() + ", ";
	ret += mats.back()->get_name();
	return ret;
}

class merge_cols_op: public detail::portion_mapply_op
{
public:
	merge_cols_op(size_t num_rows, size_t num_cols,
			const scalar_type &type): detail::portion_mapply_op(
				num_rows, num_cols, type) {
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	virtual portion_mapply_op::const_ptr transpose() const;
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return cat_mat_names(mats);
	}
};

class merge_rows_op: public detail::portion_mapply_op
{
public:
	merge_rows_op(size_t num_rows, size_t num_cols,
			const scalar_type &type): detail::portion_mapply_op(
				num_rows, num_cols, type) {
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	virtual portion_mapply_op::const_ptr transpose() const;
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return cat_mat_names(mats);
	}
};

detail::portion_mapply_op::const_ptr merge_rows_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new merge_cols_op(
				get_out_num_cols(), get_out_num_rows(), get_output_type()));
}

detail::portion_mapply_op::const_ptr merge_cols_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new merge_rows_op(
				get_out_num_cols(), get_out_num_rows(), get_output_type()));
}

void merge_rows_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(out.store_layout() == matrix_layout_t::L_ROW);
	detail::local_row_matrix_store &row_out
		= static_cast<detail::local_row_matrix_store &>(out);

	size_t out_idx = 0;
	for (size_t i = 0; i < ins.size(); i++) {
		assert(ins[i]->store_layout() == matrix_layout_t::L_ROW);
		const detail::local_row_matrix_store &row_in
			= static_cast<const detail::local_row_matrix_store &>(*ins[i]);
		for (size_t idx = 0; idx < row_in.get_num_rows(); idx++) {
			assert(out_idx < out.get_num_rows());
			memcpy(row_out.get_row(out_idx), row_in.get_row(idx),
					row_in.get_num_cols() * row_in.get_entry_size());
			out_idx++;
		}
	}
	assert(out_idx == out.get_num_rows());
}

void merge_cols_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(out.store_layout() == matrix_layout_t::L_COL);
	detail::local_col_matrix_store &col_out
		= static_cast<detail::local_col_matrix_store &>(out);

	size_t out_idx = 0;
	for (size_t i = 0; i < ins.size(); i++) {
		assert(ins[i]->store_layout() == matrix_layout_t::L_COL);
		const detail::local_col_matrix_store &col_in
			= static_cast<const detail::local_col_matrix_store &>(*ins[i]);
		for (size_t idx = 0; idx < col_in.get_num_cols(); idx++) {
			assert(out_idx < out.get_num_cols());
			memcpy(col_out.get_col(out_idx), col_in.get_col(idx),
					col_in.get_num_rows() * col_in.get_entry_size());
			out_idx++;
		}
	}
	assert(out_idx == out.get_num_cols());
}

}

collected_matrix_store::ptr collected_matrix_store::create(
		const std::vector<detail::matrix_store::const_ptr> &mats,
		size_t num_cols)
{
	assert(!mats.empty());
	size_t num_rows = mats[0]->get_num_rows();
	const scalar_type &type = mats[0]->get_type();
	for (size_t i = 0; i < mats.size(); i++) {
		if (num_rows != mats[i]->get_num_rows()) {
			BOOST_LOG_TRIVIAL(error)
				<< "collected cols have different numbers of rows";
			return collected_matrix_store::ptr();
		}
		if (type != mats[i]->get_type()) {
			BOOST_LOG_TRIVIAL(error) << "collected cols have different types";
			return collected_matrix_store::ptr();
		}
	}
	return ptr(new collected_matrix_store(mats, num_cols));
}

collected_matrix_store::collected_matrix_store(
		const std::vector<detail::matrix_store::const_ptr> &mats,
		size_t num_cols): virtual_matrix_store(mats[0]->get_num_rows(),
			num_cols, mats[0]->is_in_mem(), mats[0]->get_type())
{
	size_t num_rows = mats[0]->get_num_rows();
	const scalar_type &type = mats[0]->get_type();
	num_nodes = mats[0]->get_num_nodes();
	for (size_t i = 0; i < mats.size(); i++) {
		// The input matrix might be collected_matrix_store.
		// For the sake of performance, we want to collect the original
		// matrix stores.
		const collected_matrix_store *tmp
			= dynamic_cast<const collected_matrix_store *>(
					mats[i].get());
		if (tmp)
			orig_mats.insert(orig_mats.end(), tmp->orig_mats.begin(),
					tmp->orig_mats.end());
		else
			orig_mats.push_back(mats[i]);

		for (size_t j = 0; j < mats[i]->get_num_cols(); j++)
			rc_idxs.push_back(rc_idx(i, j));
	}
	assert(rc_idxs.size() <= num_cols);

	merged_mat = fm::detail::__mapply_portion_virtual(orig_mats,
			detail::portion_mapply_op::const_ptr(new merge_cols_op(
					num_rows, rc_idxs.size(), type)),
			mats[0]->store_layout());
}

collected_matrix_store::collected_matrix_store(size_t num_rows,
		size_t num_cols, bool in_mem, const scalar_type &type,
		size_t num_nodes): virtual_matrix_store(num_rows, num_cols, in_mem, type)
{
	this->num_nodes = num_nodes;
}

detail::matrix_store::const_ptr collected_matrix_store::get_cols(
		const std::vector<off_t> &idxs) const
{
	assert(std::is_sorted(idxs.begin(), idxs.end()));
	std::vector<detail::matrix_store::const_ptr> ret_mats;

	size_t col_idx = 0;
	size_t mat_idx = rc_idxs[col_idx].mat_idx;
	std::vector<off_t> inner_idxs;
	while (col_idx < idxs.size()) {
		// If the current still belongs to the same matrix, record the column
		// index and continue.
		if (rc_idxs[col_idx].mat_idx == mat_idx)
			inner_idxs.push_back(rc_idxs[col_idx].inner_idx);
		else {
			// We now get to another matrix.
			ret_mats.push_back(orig_mats[mat_idx]->get_cols(inner_idxs));

			mat_idx = rc_idxs[col_idx].mat_idx;
			inner_idxs.clear();
			inner_idxs.push_back(rc_idxs[col_idx].inner_idx);
		}

		col_idx++;
	}
	ret_mats.push_back(orig_mats[mat_idx]->get_cols(inner_idxs));

	size_t num_cols = 0;
	for (size_t i = 0; i < ret_mats.size(); i++)
		num_cols += ret_mats[i]->get_num_cols();
	assert(num_cols == idxs.size());

	return collected_matrix_store::create(ret_mats, num_cols);
}

detail::matrix_store::const_ptr collected_matrix_store::transpose() const
{
	collected_matrix_store *store = new collected_matrix_store(
			get_num_cols(), get_num_rows(), is_in_mem(), get_type(),
			get_num_nodes());
	store->merged_mat = merged_mat->transpose();
	store->orig_mats.resize(orig_mats.size());
	store->rc_idxs = rc_idxs;
	for (size_t i = 0; i < orig_mats.size(); i++)
		store->orig_mats[i] = orig_mats[i]->transpose();
	return detail::matrix_store::const_ptr(store);
}

}

}
