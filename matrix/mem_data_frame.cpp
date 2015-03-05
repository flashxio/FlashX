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

namespace fm
{

vector_vector::ptr mem_data_frame::groupby(const std::string &col_name,
		gr_apply_operate<data_frame> &op) const
{
	vector::ptr col = get_vec(col_name);
	if (col == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"The column %1% doesn't exist") % col_name;
		return vector_vector::ptr();
	}
	mem_vector::ptr sorted_col = mem_vector::cast(col->deep_copy());
	type_mem_vector<off_t>::ptr idxs = type_mem_vector<off_t>::cast(
			sorted_col->sort_with_index());

	mem_data_frame::ptr sorted_df = mem_data_frame::create();
	sorted_df->add_vec(col_name, sorted_col);
	for (size_t i = 0; i < get_num_vecs(); i++) {
		mem_vector::ptr mem_vec = mem_vector::cast(get_vec(i));
		if (mem_vec == col)
			continue;
		sorted_df->add_vec(get_vec_name(i), mem_vec->get(*idxs));
	}

	vector_vector::ptr ret = std::static_pointer_cast<vector_vector>(
			op.get_output_type().create_mem_vec_vec());
	mem_vector::ptr row = op.get_output_type().create_mem_vec(0);
	const agg_operate &find_next
		= sorted_col->get_type().get_agg_ops().get_find_next();
	size_t loc = 0;
	const mem_vector *const_sorted_col = sorted_col.get();
	size_t col_len = sorted_col->get_length();
	const char *start = const_sorted_col->get_raw_arr();
	size_t entry_size = sorted_col->get_entry_size();
	while (loc < col_len) {
		size_t curr_length = col_len - loc;
		const char *curr_ptr = start + entry_size * loc;
		size_t rel_end;
		find_next.run(curr_length, curr_ptr, &rel_end);
		// This expose a portion of the data frame.
		sorted_df->expose_portion(loc, rel_end);
		// The first argument is the key and the second one is the value
		// (a data frame)
		op.run(curr_ptr, *sorted_df, *row);
		if (row->get_length() > 0)
			ret->append(*row);
		loc += rel_end;
	}

	return ret;
}

}
