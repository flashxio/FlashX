#ifndef __MATRIX_IO_H__
#define __MATRIX_IO_H__

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

#include <stdlib.h>

#include <memory>

#include "vertex.h"
#include "sparse_matrix.h"

class matrix_loc
{
	std::pair<off_t, off_t> pair;
public:
	matrix_loc(off_t row_idx, off_t col_idx): pair(row_idx, col_idx) {
	}

	off_t get_row_idx() const {
		return pair.first;
	}

	off_t get_col_idx() const {
		return pair.second;
	}
};

class matrix_io
{
	matrix_loc top_left;
	size_t num_rows;
	size_t num_cols;

	data_loc_t loc;
	size_t size;
public:
	matrix_io(const matrix_loc &_top_left, size_t num_rows,
			size_t num_cols, data_loc_t loc, size_t size): top_left(_top_left) {
		this->num_rows = num_rows;
		this->num_cols = num_cols;
		this->loc = loc;
		this->size = size;
	}

	const data_loc_t &get_loc() const {
		return loc;
	}

	size_t get_size() const {
		return size;
	}

	const matrix_loc &get_top_left() const {
		return top_left;
	}

	size_t get_num_rows() const {
		return num_rows;
	}

	size_t get_num_cols() const {
		return num_cols;
	}
};

class sparse_matrix;

class matrix_io_generator
{
public:
	typedef std::shared_ptr<matrix_io_generator> ptr;

	virtual matrix_io get_next_io() = 0;
	virtual bool has_next_io() = 0;
};

// In the number of row blocks.
static const size_t MIN_ROW_IO_SIZE = 1024;

// TODO we need to test the code.
class row_io_generator: public matrix_io_generator
{
	const std::vector<row_block> &blocks;
	// The current offset in the row_block vector.
	off_t curr_block_off;
	// The number of I/O generators there exist for the matrix.
	int num_gens;
	int file_id;
	size_t tot_num_rows;
	size_t tot_num_cols;

	row_io_generator(const std::vector<row_block> &_blocks, size_t tot_num_rows,
			size_t tot_num_cols, int file_id, int gen_id,
			int num_gens): blocks(_blocks) {
		curr_block_off = gen_id * MIN_ROW_IO_SIZE;
		this->tot_num_rows = tot_num_rows;
		this->tot_num_cols = tot_num_cols;
		this->num_gens = num_gens;
		this->file_id = file_id;
	}

	// Get the current row in the matrix.
	off_t get_curr_row() const {
		return curr_block_off * ROW_BLOCK_SIZE;
	}

	// Get the offset of the current row I/O in the matrix file.
	off_t get_curr_offset() const {
		return blocks[curr_block_off].get_offset();
	}

	// Get the number of rows in the current row I/O.
	size_t get_curr_num_rows() const {
		// `next_off' identifies the row block behind the I/O request.
		off_t next_off = curr_block_off + MIN_ROW_IO_SIZE;
		// blocks[blocks.size() - 1] is an empty block. It only indicates
		// the end of the matrix file.
		// blocks[blocks.size() - 2] is the last row block. It's possible that
		// it's smaller than the full-size row block.

		// If the row block behind the I/O request is the last row block,
		// we get the maximal number of rows in an I/O request.
		if (next_off < blocks.size() - 1)
			return MIN_ROW_IO_SIZE * ROW_BLOCK_SIZE;
		else
			return tot_num_rows - curr_block_off * ROW_BLOCK_SIZE;
	}

	// Get the current row I/O size.
	size_t get_curr_size() const {
		off_t next_off = curr_block_off + MIN_ROW_IO_SIZE;
		next_off = next_off < blocks.size() ? next_off : blocks.size() - 1;
		return blocks[next_off].get_offset()
			- blocks[curr_block_off].get_offset();
	}
public:
	static matrix_io_generator::ptr create(
			const std::vector<row_block> &_blocks, size_t tot_num_rows,
			size_t tot_num_cols, int file_id, int gen_id, int num_gens) {
		return matrix_io_generator::ptr(new row_io_generator(_blocks,
					tot_num_rows, tot_num_cols, file_id, gen_id, num_gens));
	}

	virtual matrix_io get_next_io() {
		matrix_loc top_left(get_curr_row(), 0);
		matrix_io ret(top_left, get_curr_num_rows(), tot_num_cols,
					data_loc_t(file_id, get_curr_offset()), get_curr_size());
		// We now need to jump to the next row block that belongs to
		// the current I/O generator.
		curr_block_off += MIN_ROW_IO_SIZE * num_gens;
		return ret;
	}

	virtual bool has_next_io() {
		// The last entry in the row_block vector is the size of the matrix
		// file.
		return curr_block_off < blocks.size() - 1;
	}
};

#endif
