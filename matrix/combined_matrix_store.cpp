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
	const scalar_type &type = mats.front()->get_type();
	bool is_wide = mats.front()->is_wide();
	for (size_t i = 1; i < mats.size(); i++) {
		if (mats[i]->get_type() != type) {
			BOOST_LOG_TRIVIAL(error)
				<< "can't combine matrices of different element types";
			return combined_matrix_store::ptr();
		}
		if (mats[i]->is_wide() != is_wide) {
			BOOST_LOG_TRIVIAL(error)
				<< "can't combine matrices of different row/col length";
			return combined_matrix_store::ptr();
		}
	}
	if (is_wide) {
		size_t num_cols = mats.front()->get_num_cols();
		for (size_t i = 1; i < mats.size(); i++) {
			if (mats[i]->get_num_cols() != num_cols) {
				BOOST_LOG_TRIVIAL(error)
					<< "can't combine matrices with different row lengths";
				return combined_matrix_store::ptr();
			}
		}
		return ptr(new combined_matrix_store(mats, layout));
	}
	else {
		size_t num_rows = mats.front()->get_num_rows();
		for (size_t i = 1; i < mats.size(); i++) {
			if (mats[i]->get_num_rows() != num_rows) {
				BOOST_LOG_TRIVIAL(error)
					<< "can't combine matrices with different col lengths";
				return combined_matrix_store::ptr();
			}
		}
		return ptr(new combined_matrix_store(mats, layout));
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
	for (size_t i = 0; i < ins.size(); i++) {
		out.resize(0, col_idx, out.get_num_rows(), ins[i]->get_num_cols());
		out.copy_from(*ins[i]);
		col_idx += ins[i]->get_num_cols();
	}
	out.reset_size();
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
	for (size_t i = 0; i < ins.size(); i++) {
		out.resize(row_idx, 0, ins[i]->get_num_rows(), out.get_num_cols());
		// TODO I need to test the performance of memory copy
		// It might be slow.
		out.copy_from(*ins[i]);
		row_idx += ins[i]->get_num_rows();
	}
	out.reset_size();
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

std::string combined_matrix_store::get_name() const
{
	std::string name = std::string("combine(") + mats[0]->get_name();
	for (size_t i = 1; mats.size(); i++)
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
	return matrix_store::const_ptr(new combined_matrix_store(tmp, layout));
}

matrix_store::const_ptr combined_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	// TODO
	return matrix_store::const_ptr();
}

}

}
