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
#include <unordered_set>

#include "combined_matrix_store.h"
#include "dense_matrix.h"
#include "local_matrix_store.h"

namespace fm
{

namespace detail
{

struct mat_info
{
	size_t nrow;
	size_t ncol;
	bool in_mem;
};

static mat_info get_matrix_info(const std::vector<matrix_store::const_ptr> &mats)
{
	mat_info info;
	bool is_wide = mats.front()->is_wide();
	info.nrow = mats.front()->get_num_rows();
	info.ncol = mats.front()->get_num_cols();
	info.in_mem = mats.front()->is_in_mem();

	for (size_t i = 1; i < mats.size(); i++) {
		if (is_wide)
			info.nrow += mats[i]->get_num_rows();
		else
			info.ncol += mats[i]->get_num_cols();
		info.in_mem = info.in_mem && mats[i]->is_in_mem();
	}
	return info;
}

combined_matrix_store::ptr combined_matrix_store::create(
		const std::vector<matrix_store::const_ptr> &mats, matrix_layout_t layout)
{
	if (mats.empty()) {
		BOOST_LOG_TRIVIAL(error) << "can't combine 0 matrices";
		return combined_matrix_store::ptr();
	}

	std::vector<matrix_store::const_ptr> tmp_mats;
	for (size_t i = 0; i < mats.size(); i++) {
		auto tmp = std::dynamic_pointer_cast<const combined_matrix_store>(mats[i]);
		if (tmp == NULL)
			tmp_mats.push_back(mats[i]);
		else {
			// We always flatten a combined matrix, so we don't need to check
			// the matrix stores inside a combined matrix.
			for (size_t j = 0; j < tmp->mats.size(); j++)
				tmp_mats.push_back(tmp->mats[j]);
		}
	}
	const scalar_type &type = tmp_mats.front()->get_type();
	bool is_wide = tmp_mats.front()->is_wide();
	for (size_t i = 1; i < tmp_mats.size(); i++) {
		if (tmp_mats[i]->get_type() != type) {
			BOOST_LOG_TRIVIAL(error)
				<< "can't combine matrices of different element types";
			return combined_matrix_store::ptr();
		}
		if (tmp_mats[i]->is_wide() != is_wide) {
			BOOST_LOG_TRIVIAL(error)
				<< "can't combine matrices of different row/col length";
			return combined_matrix_store::ptr();
		}
	}
	if (is_wide) {
		size_t num_cols = tmp_mats.front()->get_num_cols();
		for (size_t i = 1; i < tmp_mats.size(); i++) {
			if (tmp_mats[i]->get_num_cols() != num_cols) {
				BOOST_LOG_TRIVIAL(error)
					<< "can't combine matrices with different row lengths";
				return combined_matrix_store::ptr();
			}
		}
		return ptr(new combined_matrix_store(tmp_mats, layout));
	}
	else {
		size_t num_rows = tmp_mats.front()->get_num_rows();
		for (size_t i = 1; i < tmp_mats.size(); i++) {
			if (tmp_mats[i]->get_num_rows() != num_rows) {
				BOOST_LOG_TRIVIAL(error)
					<< "can't combine matrices with different col lengths";
				return combined_matrix_store::ptr();
			}
		}
		return ptr(new combined_matrix_store(tmp_mats, layout));
	}
}

namespace
{

class combine_tall_op: public portion_mapply_op
{
public:
	combine_tall_op(size_t num_rows, size_t num_cols,
			const scalar_type &t): portion_mapply_op(num_rows, num_cols, t) {
	}

	virtual portion_mapply_op::const_ptr transpose() const;
	virtual void run(const std::vector<local_matrix_store::const_ptr> &ins,
			local_matrix_store &out) const;
	virtual std::string to_string(
			const std::vector<matrix_store::const_ptr> &mats) const {
		// TODO
		return "";
	}
};

class combine_wide_op: public portion_mapply_op
{
public:
	combine_wide_op(size_t num_rows, size_t num_cols,
			const scalar_type &t): portion_mapply_op(num_rows, num_cols, t) {
	}

	virtual portion_mapply_op::const_ptr transpose() const;
	virtual void run(const std::vector<local_matrix_store::const_ptr> &ins,
			local_matrix_store &out) const;
	virtual std::string to_string(
			const std::vector<matrix_store::const_ptr> &mats) const {
		// TODO
		return "";
	}
};

portion_mapply_op::const_ptr combine_tall_op::transpose() const
{
	return portion_mapply_op::const_ptr(new combine_wide_op(get_out_num_cols(),
				get_out_num_rows(), get_output_type()));
}

void combine_tall_op::run(const std::vector<local_matrix_store::const_ptr> &ins,
		local_matrix_store &out) const
{
	size_t col_idx = 0;
	local_matrix_store::exposed_area area = out.get_exposed_area();
	off_t local_start_row = out.get_local_start_row();
	off_t local_start_col = out.get_local_start_col();
	for (size_t i = 0; i < ins.size(); i++) {
		out.resize(local_start_row, local_start_col + col_idx,
				out.get_num_rows(), ins[i]->get_num_cols());
		out.copy_from(*ins[i]);
		col_idx += ins[i]->get_num_cols();
	}
	out.restore_size(area);
	assert(col_idx == out.get_num_cols());
}

portion_mapply_op::const_ptr combine_wide_op::transpose() const
{
	return portion_mapply_op::const_ptr(new combine_tall_op(get_out_num_cols(),
				get_out_num_rows(), get_output_type()));
}

void combine_wide_op::run(const std::vector<local_matrix_store::const_ptr> &ins,
		local_matrix_store &out) const
{
	size_t row_idx = 0;
	local_matrix_store::exposed_area area = out.get_exposed_area();
	off_t local_start_row = out.get_local_start_row();
	off_t local_start_col = out.get_local_start_col();
	for (size_t i = 0; i < ins.size(); i++) {
		out.resize(local_start_row + row_idx, local_start_col,
				ins[i]->get_num_rows(), out.get_num_cols());
		// TODO I need to test the performance of memory copy
		// It might be slow.
		out.copy_from(*ins[i]);
		row_idx += ins[i]->get_num_rows();
	}
	out.restore_size(area);
	assert(row_idx == out.get_num_rows());
}

}

static portion_mapply_op::const_ptr get_combine_op(
		const std::vector<matrix_store::const_ptr> &mats)
{
	mat_info info = get_matrix_info(mats);
	const scalar_type &type = mats.front()->get_type();
	if (mats.front()->is_wide())
		return portion_mapply_op::const_ptr(new combine_wide_op(info.nrow,
					info.ncol, type));
	else
		return portion_mapply_op::const_ptr(new combine_tall_op(info.nrow,
					info.ncol, type));
}

combined_matrix_store::combined_matrix_store(
		const std::vector<matrix_store::const_ptr> &mats,
		matrix_layout_t layout): mapply_matrix_store(mats,
			get_combine_op(mats), layout)
{
	this->mats = mats;
}

combined_matrix_store::combined_matrix_store(
		const std::vector<matrix_store::const_ptr> &mats,
		matrix_layout_t layout, data_id_t::ptr data_id): mapply_matrix_store(
			mats, get_combine_op(mats), layout, data_id)
{
	this->mats = mats;
}

std::string combined_matrix_store::get_name() const
{
	std::string name = std::string("combine(") + mats[0]->get_name();
	for (size_t i = 1; i < mats.size(); i++)
		name += std::string(", ") + mats[i]->get_name();
	name += ")";
	return name;
}

matrix_store::const_ptr combined_matrix_store::transpose() const
{
	std::vector<matrix_store::const_ptr> tmp(mats.size());
	for (size_t i = 0; i < tmp.size(); i++)
		tmp[i] = mats[i]->transpose();
	matrix_layout_t layout;
	if (store_layout() == matrix_layout_t::L_ROW)
		layout = matrix_layout_t::L_COL;
	else
		layout = matrix_layout_t::L_ROW;
	return matrix_store::const_ptr(new combined_matrix_store(tmp, layout,
				get_id_ptr()));
}

/*
 * Find the required row from the vector of matrices.
 * It returns the matrix idx of the matrix with the row and
 * the local row idx in the matrix.
 */
static inline std::pair<off_t, off_t> find_mat(
		std::vector<matrix_store::const_ptr>::const_iterator start,
		std::vector<matrix_store::const_ptr>::const_iterator end,
		off_t start_row, off_t row)
{
	for (auto it = start; it != end; it++) {
		if (row >= start_row
				&& (size_t) row < start_row + (*it)->get_num_rows()) {
			return std::pair<off_t, off_t>(it - start, row - start_row);
		}
		start_row += (*it)->get_num_rows();
	}
	return std::pair<off_t, off_t>(-1, -1);
}

/*
 * Find the last location that has the same value as the one at `idx'.
 */
off_t find_last(const std::vector<off_t> &vals, off_t idx)
{
	for (size_t i = idx; i < vals.size(); i++) {
		if (vals[i] != vals[idx])
			return i - 1;
	}
	return vals.size() - 1;
}

static inline matrix_store::const_ptr _get_rows(matrix_store::const_ptr store,
		const std::vector<off_t> &idxs)
{
	matrix_store::const_ptr rows = store->get_rows(idxs);
	if (rows)
		return rows;

	// If we can't get rows from a matrix directly, we can extract them with
	// get_rows of dense_matrix.
	dense_matrix::ptr mat = dense_matrix::create(store);
	dense_matrix::ptr sub = mat->get_rows(idxs);
	assert(sub);
	return sub->get_raw_store();
}

matrix_store::const_ptr combined_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	// Get rows from a tall matrix.
	if (!is_wide()) {
		if (idxs.size() > get_num_cols()) {
			std::vector<matrix_store::const_ptr> sub_stores(mats.size());
			for (size_t i = 0; i < sub_stores.size(); i++)
				sub_stores[i] = mats[i]->get_rows(idxs);
			return combined_matrix_store::create(sub_stores, store_layout());
		}
		else {
			mem_matrix_store::ptr ret = mem_matrix_store::create(idxs.size(),
					get_num_cols(), store_layout(), get_type(), -1);
			size_t start_col = 0;
			for (size_t i = 0; i < mats.size(); i++) {
				matrix_store::const_ptr sub_store = mats[i]->get_rows(idxs);
				local_matrix_store::const_ptr portion = sub_store->get_portion(0);
				assert(sub_store->get_num_rows() == portion->get_num_rows());
				assert(sub_store->get_num_cols() == portion->get_num_cols());
				local_matrix_store::ptr ret_portion = ret->get_portion(0,
						start_col, idxs.size(), sub_store->get_num_cols());
				ret_portion->copy_from(*portion);
				start_col += sub_store->get_num_cols();
			}
			return ret;
		}

		// TODO I'll deal with this later.
		assert(0);
	}

	// If all the idxs are sorted.
	if (std::is_sorted(idxs.begin(), idxs.end())) {
		std::vector<off_t> mat_idxs(idxs.size());
		std::vector<off_t> local_row_idxs(idxs.size());
		std::vector<off_t> start_rows(mats.size());
		for (size_t i = 1; i < mats.size(); i++)
			start_rows[i] = start_rows[i - 1] + mats[i - 1]->get_num_rows();
		// Find the right matrix.
		auto loc = find_mat(mats.begin(), mats.end(), start_rows[0], idxs[0]);
		mat_idxs[0] = loc.first;
		local_row_idxs[0] = loc.second;
		for (size_t i = 1; i < idxs.size(); i++) {
			off_t mat_idx = mat_idxs[i - 1];
			loc = find_mat(mats.begin() + mat_idx, mats.end(),
					start_rows[mat_idx], idxs[i]);
			assert(loc.first >= 0);
			local_row_idxs[i] = loc.second;
			// find_mat returns a relative idx in the vector of matrices.
			mat_idxs[i] = loc.first + mat_idx;
		}

		std::vector<matrix_store::const_ptr> ret;
		size_t i = 0;
		while (i < idxs.size()) {
			size_t last = find_last(mat_idxs, i);
			size_t local_nrow = last - i + 1;
			size_t mat_idx = mat_idxs[i];

			// Test if we need all rows from a matrix.
			bool all_rows = false;
			if (local_nrow == mats[mat_idx]->get_num_rows()) {
				std::unordered_set<off_t> uniq(local_row_idxs.begin() + i,
						local_row_idxs.begin() + i + local_nrow);
				all_rows = uniq.size() == local_row_idxs.size();
			}

			if (all_rows)
				ret.push_back(mats[mat_idx]);
			else {
				std::vector<off_t> tmp_rows(local_row_idxs.begin() + i,
						local_row_idxs.begin() + i + local_nrow);
				auto rows = _get_rows(mats[mat_idx], tmp_rows);
				ret.push_back(rows);
			}
			i = last + 1;
		}
		return combined_matrix_store::create(ret, store_layout());
	}

	// Figure out the physical matrix and the local row index for each row.
	std::vector<off_t> all_mat_idxs(get_num_rows());
	std::vector<off_t> all_lrow_idxs(get_num_rows());
	off_t mat_idx = 0;
	off_t lrow_idx = 0;
	for (size_t i = 0; i < get_num_rows(); i++) {
		assert((size_t) mat_idx < mats.size());
		assert((size_t) lrow_idx < mats[mat_idx]->get_num_rows());
		all_mat_idxs[i] = mat_idx;
		all_lrow_idxs[i] = lrow_idx;
		lrow_idx++;
		if ((size_t) lrow_idx == mats[mat_idx]->get_num_rows()) {
			mat_idx++;
			lrow_idx = 0;
		}
	}

	// The requried rows aren't in the sorted order. We need to figure
	// out the physical matrix and the row index in the matrix for each
	// required row.
	std::vector<off_t> mat_idxs(idxs.size());
	std::vector<off_t> lrow_idxs(idxs.size());
	bool same_mat = true;
	for (size_t i = 0; i < idxs.size(); i++) {
		if ((size_t) idxs[i] >= get_num_rows()) {
			BOOST_LOG_TRIVIAL(error) << "row idxs are out of bounds";
			return matrix_store::const_ptr();
		}
		mat_idxs[i] = all_mat_idxs[idxs[i]];
		lrow_idxs[i] = all_lrow_idxs[idxs[i]];
		// If they aren't in the same matrix.
		if (mat_idxs[i] != mat_idxs[0])
			same_mat = false;
	}

	if (same_mat)
		return _get_rows(mats[mat_idxs[0]], lrow_idxs);

	// If rows are from different matrices, we get rows individually and
	// then combine them to form a new matrix.
	std::vector<matrix_store::const_ptr> rows(idxs.size());
	for (size_t i = 0; i < idxs.size(); i++) {
		std::vector<off_t> lrow_idx(1, lrow_idxs[i]);
		rows[i] = _get_rows(mats[mat_idxs[i]], lrow_idx);
	}
	return combined_matrix_store::create(rows, matrix_layout_t::L_ROW);
}

bool combined_matrix_store::is_sparse() const
{
	for (size_t i = 0; i < mats.size(); i++)
		if (mats[i]->is_sparse())
			return true;
	return false;
}

bool combined_matrix_store::set_persistent(const std::string &name) const
{
	std::vector<const EM_object *> objs;
	for (size_t i = 0; i < mats.size(); i++) {
		const EM_object *obj = dynamic_cast<const EM_object *>(mats[i].get());
		if (obj)
			objs.push_back(obj);
	}

	for (size_t i = 0; i < objs.size(); i++) {
		bool ret = objs[i]->set_persistent(name + "-" + ltoa(i));
		if (!ret) {
			// Unset all previous matrices.
			for (size_t j = 0; j < i; j++)
				objs[j]->unset_persistent();
			return false;
		}
	}
	return true;
}

void combined_matrix_store::unset_persistent() const
{
	for (size_t i = 0; i < mats.size(); i++) {
		const EM_object *obj = dynamic_cast<const EM_object *>(mats[i].get());
		if (obj)
			obj->unset_persistent();
	}
}

bool combined_matrix_store::share_data(const matrix_store &store) const
{
	const combined_matrix_store *combined_store
		= dynamic_cast<const combined_matrix_store *>(&store);
	if (combined_store == NULL)
		return false;
	if (mats.size() != combined_store->mats.size())
		return false;
	for (size_t i = 0; i < mats.size(); i++)
		if (!mats[i]->share_data(*combined_store->mats[i]))
			return false;
	return true;
}

}

}
