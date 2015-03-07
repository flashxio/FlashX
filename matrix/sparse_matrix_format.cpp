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

#include "sparse_matrix_format.h"

namespace fm
{

void sparse_block_2d::verify(const block_2d_size &block_size) const
{
	assert(num_rows <= block_size.get_num_rows() && num_rows > 0);
	row_part_iterator it = get_iterator();
	size_t rel_row_id = it.get_curr().get_rel_row_idx();
	while (it.has_next()) {
		const sparse_row_part &part = it.next();
		assert(part.get_num_non_zeros() <= block_size.get_num_cols());
		assert(rel_row_id < part.get_rel_row_idx());
		rel_row_id = part.get_rel_row_idx();
	}
}

row_part_iterator sparse_block_2d::append(const row_part_iterator &it,
		const sparse_row_part &part)
{
	assert(it.get_row_idx() == num_rows);
	assert(it.get_row_idx() <= part.get_rel_row_idx());
	// Discard the const qualifier.
	sparse_row_part *end = const_cast<sparse_row_part *>(&it.get_curr());
	memcpy(end, &part, part.get_size());
	num_rows++;
	row_part_iterator end_it = it;
	end_it.next();
	return end_it;
}

}
