#ifndef __SPARSE_MATRIX_FORMAT_H__
#define __SPARSE_MATRIX_FORMAT_H__

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

#include "vertex.h"
#include "matrix_header.h"

namespace fm
{

/*
 * This stores part of a row. It doesn't contain any attributes.
 */
class sparse_row_part
{
	uint16_t rel_row_idx;
	uint16_t num_non_zeros;
	uint16_t rel_col_idxs[0];
public:
	static size_t get_size(size_t num_non_zeros) {
		return sizeof(sparse_row_part) + sizeof(uint16_t) * num_non_zeros;
	}

	static size_t get_col_entry_size() {
		return sizeof(uint16_t);
	}

	sparse_row_part(uint16_t rel_row_idx) {
		this->rel_row_idx = rel_row_idx;
		num_non_zeros = 0;
	}

	size_t get_rel_row_idx() const {
		return rel_row_idx;
	}

	size_t get_num_non_zeros() const {
		return num_non_zeros;
	}

	void add(const block_2d_size &block_size, size_t cidx) {
		rel_col_idxs[num_non_zeros] = cidx & block_size.get_ncol_mask();
		num_non_zeros++;
	}

	size_t get_size() const {
		return get_size(num_non_zeros);
	}
};

/*
 * This iterates row parts in a 2D-partitioned block. It only allows reads.
 */
class row_part_iterator
{
	size_t num_rows;
	size_t curr_idx;
	const char *curr_row_part;
public:
	row_part_iterator(const char *row_part_start, size_t num_rows) {
		this->num_rows = num_rows;
		curr_idx = 0;
		curr_row_part = row_part_start;
	}

	size_t get_row_idx() const {
		return curr_idx;
	}

	bool has_next() const {
		return curr_idx < num_rows - 1;
	}

	const sparse_row_part &get_curr() const {
		return *(sparse_row_part *) curr_row_part;
	}

	const sparse_row_part &next() {
		curr_row_part += get_curr().get_size();
		curr_idx++;
		return *(sparse_row_part *) curr_row_part;
	}
};

/*
 * This stores the block of a sparse matrix partitioned in 2D.
 * The block is stored in row major.
 */
class sparse_block_2d
{
	uint32_t block_row_idx;
	uint32_t block_col_idx;
	uint32_t num_rows;
	// This is where the row parts are serialized.
	char row_parts[0];
public:
	sparse_block_2d(uint32_t block_row_idx, uint32_t block_col_idx) {
		this->block_row_idx = block_row_idx;
		this->block_col_idx = block_col_idx;
		num_rows = 0;
	}

	row_part_iterator get_iterator() const {
		return row_part_iterator(row_parts, num_rows);
	}

	row_part_iterator append(const row_part_iterator &it,
			const sparse_row_part &part);

	/*
	 * Calculate the offset (in bytes) of the current row
	 * in this 2D partitioned block.
	 */
	size_t get_byte_off(const row_part_iterator &it) const {
		return ((char *) &it.get_curr()) - ((char *) this);
	}

	void verify(const block_2d_size &block_size) const;
};

}

#endif
