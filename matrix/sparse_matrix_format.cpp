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

#include "sparse_matrix_format.h"
#include "matrix_config.h"

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

sparse_matrix_index::ptr sparse_matrix_index::create(const matrix_header &header,
		const std::vector<off_t> &offs)
{
	block_2d_size block_size = header.get_2d_block_size();
	if (offs.size() != block_size.cal_num_row_blocks(header.get_num_rows()) + 1) {
		BOOST_LOG_TRIVIAL(error) << "There are an incorrect number of offsets";
		return sparse_matrix_index::ptr();
	}
	void *buf = NULL;
	int ret = posix_memalign(&buf, PAGE_SIZE, get_size(offs.size()));
	if (ret) {
		BOOST_LOG_TRIVIAL(error) << "Can't allocate memory";
		return sparse_matrix_index::ptr();
	}
	sparse_matrix_index *idx = new (buf) sparse_matrix_index(header);
	for (size_t i = 0; i < offs.size(); i++) 
		idx->offs[i] = offs[i];
	return sparse_matrix_index::ptr(idx, deleter());
}

size_t sparse_matrix_index::get_num_entries() const
{
	block_2d_size block_size = header.get_2d_block_size();
	return block_size.cal_num_row_blocks(header.get_num_rows()) + 1;
}

void sparse_matrix_index::dump(const std::string &file) const
{
	FILE *f = fopen(file.c_str(), "w");
	if (f == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format("can't open %1%: %2%")
			% file % strerror(errno);
		return;
	}

	size_t ret = fwrite(this, get_size(get_num_entries()), 1, f);
	if (ret == 0) {
		BOOST_LOG_TRIVIAL(error) << boost::format("can't write to %1%: %2%")
			% file % strerror(errno);
		return;
	}
	fclose(f);
}

}
