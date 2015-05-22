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

#include <boost/format.hpp>

#include "log.h"

#include "mem_data_frame.h"
#include "mem_vector.h"
#include "mem_vector_vector.h"
#include "local_vec_store.h"

namespace fm
{

void expose_portion(const mem_data_frame &sorted_df, off_t loc, size_t length,
		sub_data_frame &sub_df)
{
	// TODO This is an very inefficient implementation.
	sub_df.resize(sorted_df.get_num_vecs());
	for (size_t i = 0; i < sub_df.size(); i++)
		sub_df[i] = sorted_df.get_vec_ref(i).get_portion(loc, length);
}

vector_vector::ptr mem_data_frame::groupby(const std::string &col_name,
		gr_apply_operate<sub_data_frame> &op) const
{
	const detail::vec_store &col = get_vec_ref(col_name);
	detail::mem_vec_store::ptr sorted_col;
	mem_data_frame::ptr sorted_df;
	if (!col.is_sorted()) {
		// If the column where we group by isn't sorted, we'll create a copy
		// of this data frame and sort on the replicated data frame.
		sorted_col = detail::mem_vec_store::cast(col.deep_copy());
		detail::mem_vec_store::ptr idxs = detail::mem_vec_store::cast(
				sorted_col->sort_with_index());

		sorted_df = mem_data_frame::create();
		sorted_df->add_vec(col_name, sorted_col);
		for (size_t i = 0; i < get_num_vecs(); i++) {
			if (get_vec_name(i) == col_name)
				continue;
			const detail::mem_vec_store &mem_vec
				= static_cast<const detail::mem_vec_store &>(get_vec_ref(i));
			sorted_df->add_vec(get_vec_name(i), mem_vec.get(*idxs));
		}
	}
	else {
		// We aren't going to change the data in the vector.
		// Discard const qualifer.
		sorted_col = detail::mem_vec_store::cast(
				const_cast<detail::vec_store &>(col).shallow_copy());
		// If the column has been sorted, we still need to create
		// a data frame, but it'll reference the vectors in the original
		// data frame.
		sorted_df = mem_data_frame::create();
		sorted_df->add_vec(col_name, sorted_col);
		for (size_t i = 0; i < get_num_vecs(); i++) {
			if (get_vec_name(i) == col_name)
				continue;
			detail::vec_store &vec = const_cast<detail::vec_store &>(
					get_vec_ref(i));
			sorted_df->add_vec(get_vec_name(i), vec.shallow_copy());
		}
	}

	sub_data_frame sub_df;
	size_t out_size;
	// If the user can predict the number of output elements, we can create
	// a buffer of the expected size.
	if (op.get_num_out_eles() > 0)
		out_size = op.get_num_out_eles();
	else
		// If the user can't, we create a small buffer.
		out_size = 16;
	local_buf_vec_store row(0, out_size, op.get_output_type(), -1);

	detail::mem_vv_store::ptr ret = detail::mem_vv_store::create(
			op.get_output_type());
	const agg_operate &find_next
		= sorted_col->get_type().get_agg_ops().get_find_next();
	size_t loc = 0;
	size_t col_len = sorted_col->get_length();
	const char *start = sorted_col->get_raw_arr();
	size_t entry_size = sorted_col->get_entry_size();
	while (loc < col_len) {
		size_t curr_length = col_len - loc;
		const char *curr_ptr = start + entry_size * loc;
		size_t rel_end;
		find_next.run(curr_length, curr_ptr, &rel_end);
		// This expose a portion of the data frame.
		expose_portion(*sorted_df, loc, rel_end, sub_df);
		// The first argument is the key and the second one is the value
		// (a data frame)
		op.run(curr_ptr, sub_df, row);
		if (row.get_length() > 0)
			ret->append(row);
		loc += rel_end;
	}

	return mem_vector_vector::create(ret);
}

bool mem_data_frame::sort(const std::string &col_name)
{
	detail::vec_store::ptr sorted_col = get_vec(col_name);
	if (sorted_col == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"The column %1% doesn't exist") % col_name;
		return false;
	}
	if (sorted_col->is_sorted())
		return true;

	detail::mem_vec_store::ptr idxs = detail::mem_vec_store::cast(
			sorted_col->sort_with_index());
	for (size_t i = 0; i < get_num_vecs(); i++) {
		detail::mem_vec_store::ptr mem_vec = detail::mem_vec_store::cast(get_vec(i));
		if (mem_vec == sorted_col)
			continue;
		detail::mem_vec_store::ptr tmp = mem_vec->get(*idxs);
		assert(!tmp->is_sorted());
		set_vec(i, tmp);
	}
	return true;
}

data_frame::ptr merge_data_frame(const std::vector<data_frame::const_ptr> dfs)
{
	assert(!dfs.empty());
	size_t num_vecs = dfs[0]->get_num_vecs();
	for (size_t i = 1; i < dfs.size(); i++) {
		if (num_vecs != dfs[i]->get_num_vecs()) {
			BOOST_LOG_TRIVIAL(error)
				<< "The data frames have different numbers of vectors";
			return data_frame::ptr();
		}
	}

	mem_data_frame::ptr df = mem_data_frame::create();
	for (size_t vec_idx = 0; vec_idx < num_vecs; vec_idx++) {
		std::string vec_name = dfs[0]->named_vecs[vec_idx].first;
		const scalar_type &vec_type
			= dfs[0]->named_vecs[vec_idx].second->get_type();
		detail::vec_store::ptr vec
			= dfs[0]->named_vecs[vec_idx].second->deep_copy();

		// Here we collect the same column from all the data frame
		// except the first one.
		std::vector<detail::vec_store::const_ptr> vecs;
		for (size_t df_idx = 1; df_idx < dfs.size(); df_idx++) {
			if (vec_name != dfs[df_idx]->named_vecs[vec_idx].first
					|| vec_type != dfs[df_idx]->named_vecs[vec_idx].second->get_type()) {
				BOOST_LOG_TRIVIAL(error)
					<< "The data frames have different vectors";
				return data_frame::ptr();
			}
			vecs.push_back(dfs[df_idx]->named_vecs[vec_idx].second);
		}
		vec->append(vecs.begin(), vecs.end());
		df->add_vec(vec_name, vec);
	}
	return df;
}

bool mem_data_frame::is_sorted(const std::string &col_name) const
{
	return get_vec_ref(col_name).is_sorted();
}

}
