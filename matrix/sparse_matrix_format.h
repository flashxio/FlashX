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

	// The row part can only be allocated in the heap. The copy constructor
	// doesn't make any sense for this class.
	sparse_row_part(const sparse_row_part &part) = delete;
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

	size_t get_rel_col_idx(size_t idx) const {
		return rel_col_idxs[idx];
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
	size_t curr_idx;
	const char *curr_row_part;
	const char *rpart_end;
public:
	row_part_iterator(const char *row_part_start, size_t rparts_size) {
		curr_idx = 0;
		curr_row_part = row_part_start;
		rpart_end = row_part_start + rparts_size;
	}

	size_t get_row_idx() const {
		return curr_idx;
	}

	bool has_next() const {
		return curr_row_part < rpart_end;
	}

	const sparse_row_part &get_curr() const {
		return *(sparse_row_part *) curr_row_part;
	}

	const sparse_row_part &next() {
		const char *orig = curr_row_part;
		curr_row_part += get_curr().get_size();
		curr_idx++;
		return *(sparse_row_part *) orig;
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
	// TODO I need to make sure 32-bits are enough. Normally, they should be.
	// This is the total size of all row parts in the block.
	uint32_t rparts_size;
	// This is where the row parts are serialized.
	char row_parts[0];

	// The 2D-partitioned block has to be allowed in the heap. The copy
	// constructor doesn't make sense for it.
	sparse_block_2d(const sparse_block_2d &block) = delete;
public:
	sparse_block_2d(uint32_t block_row_idx, uint32_t block_col_idx) {
		this->block_row_idx = block_row_idx;
		this->block_col_idx = block_col_idx;
		rparts_size = 0;
	}

	size_t get_block_row_idx() const {
		return block_row_idx;
	}

	size_t get_block_col_idx() const {
		return block_col_idx;
	}

	size_t get_rparts_size() const {
		return rparts_size;
	}

	bool is_empty() const {
		return get_rparts_size() == 0;
	}

	size_t get_size() const {
		return sizeof(*this) + get_rparts_size();
	}

	row_part_iterator get_iterator() const {
		return row_part_iterator(row_parts, get_rparts_size());
	}

	void append(const sparse_row_part &part);

	/*
	 * Calculate the offset (in bytes) of the current row
	 * in this 2D partitioned block.
	 */
	size_t get_byte_off(const row_part_iterator &it) const {
		return ((char *) &it.get_curr()) - ((char *) this);
	}

	void verify(const block_2d_size &block_size) const;
};

/*
 * This allows users to iterate the blocks in a block row.
 */
class block_row_iterator
{
	const sparse_block_2d *block;
	const sparse_block_2d *end;
public:
	/*
	 * `first' is the first block in the block row.
	 * `end' is the end of the last block in the block row.
	 */
	block_row_iterator(const sparse_block_2d *first, const sparse_block_2d *end) {
		block = first;
		this->end = end;
	}

	bool has_next() const {
		return block < end;
	}

	const sparse_block_2d &get_curr() const {
		return *block;
	}

	const sparse_block_2d &next() {
		const sparse_block_2d *orig = block;
		block = (const sparse_block_2d *) ((const char *) block)
			+ block->get_size();
		return *orig;
	}
};

/*
 * This indexes a sparse matrix for easy access to the 2D-partitioned blocks.
 */
class SpM_2d_index
{
	struct deleter {
		void operator()(SpM_2d_index *p) const{
			free(p);
		}
	};

	matrix_header header;
	off_t offs[0];

	SpM_2d_index(const matrix_header &_header): header(_header) {
	}

	/*
	 * The number of offset entries in the index.
	 */
	size_t get_num_entries() const {
		return get_num_block_rows() + 1;
	}
public:
	typedef std::shared_ptr<SpM_2d_index> ptr;

	/*
	 * Calculate the storage size of the index.
	 */
	static size_t get_size(size_t num_entries) {
		return sizeof(SpM_2d_index)
			+ num_entries * sizeof(SpM_2d_index::offs[0]);
	}

	static ptr create(const matrix_header &header,
			const std::vector<off_t> &offs);
	static ptr load(const std::string &idx_file);

	size_t get_num_block_rows() const;
	void dump(const std::string &file) const;
	off_t get_block_row_off(size_t idx) const;
	const matrix_header &get_header() const {
		return header;
	}
};

class SpM_2d_storage
{
	std::unique_ptr<char[]> data;
	SpM_2d_index::ptr index;

	SpM_2d_storage(std::unique_ptr<char[]> data,
			SpM_2d_index::ptr index) {
		this->data = std::move(data);
		this->index = index;
	}
public:
	typedef std::shared_ptr<SpM_2d_storage> ptr;

	static ptr load(const std::string &mat_file,
			SpM_2d_index::ptr index);

	block_row_iterator get_block_row_it(size_t idx) const {
		char *start = data.get() + index->get_block_row_off(idx);
		char *end = data.get() + index->get_block_row_off(idx + 1);
		return block_row_iterator((const sparse_block_2d *) start,
				(const sparse_block_2d *) end);
	}

	size_t get_num_block_rows() const {
		return index->get_num_block_rows();
	}

	void verify() const;
};

}

#endif
