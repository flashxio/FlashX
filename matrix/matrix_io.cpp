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

#include "sparse_matrix.h"
#include "matrix_io.h"
#include "matrix_config.h"

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

	// Get the current row in the matrix.
	off_t get_curr_row() const {
		return curr_block_off * matrix_conf.get_row_block_size();
	}

	// Get the offset of the current row I/O in the matrix file.
	off_t get_curr_offset() const {
		return blocks[curr_block_off].get_offset();
	}

	// Get the number of rows in the current row I/O.
	size_t get_curr_num_rows() const {
		// `next_off' identifies the row block behind the I/O request.
		off_t next_off = curr_block_off + matrix_conf.get_rb_io_size();
		// blocks[blocks.size() - 1] is an empty block. It only indicates
		// the end of the matrix file.
		// blocks[blocks.size() - 2] is the last row block. It's possible that
		// it's smaller than the full-size row block.

		// If the row block behind the I/O request is the last row block,
		// we get the maximal number of rows in an I/O request.
		if (next_off < blocks.size() - 1)
			return matrix_conf.get_rb_io_size() * matrix_conf.get_row_block_size();
		else
			return tot_num_rows - curr_block_off * matrix_conf.get_row_block_size();
	}

	// Get the current row I/O size.
	size_t get_curr_size() const {
		off_t next_off = curr_block_off + matrix_conf.get_rb_io_size();
		next_off = next_off < blocks.size() ? next_off : blocks.size() - 1;
		return blocks[next_off].get_offset()
			- blocks[curr_block_off].get_offset();
	}
public:
	row_io_generator(const std::vector<row_block> &_blocks, size_t tot_num_rows,
			size_t tot_num_cols, int file_id, int gen_id,
			int num_gens): blocks(_blocks) {
		curr_block_off = gen_id * matrix_conf.get_rb_io_size();
		this->tot_num_rows = tot_num_rows;
		this->tot_num_cols = tot_num_cols;
		this->num_gens = num_gens;
		this->file_id = file_id;
	}

	virtual matrix_io get_next_io() {
		matrix_loc top_left(get_curr_row(), 0);
		matrix_io ret(top_left, get_curr_num_rows(), tot_num_cols,
					data_loc_t(file_id, get_curr_offset()), get_curr_size());
		// We now need to jump to the next row block that belongs to
		// the current I/O generator.
		curr_block_off += matrix_conf.get_rb_io_size() * num_gens;
		return ret;
	}

	virtual bool has_next_io() {
		// The last entry in the row_block vector is the size of the matrix
		// file.
		return curr_block_off < blocks.size() - 1;
	}
};

matrix_io_generator::ptr matrix_io_generator::create(
		const std::vector<row_block> &_blocks, size_t tot_num_rows,
		size_t tot_num_cols, int file_id, int gen_id, int num_gens)
{
	return matrix_io_generator::ptr(new row_io_generator(_blocks,
				tot_num_rows, tot_num_cols, file_id, gen_id, num_gens));
}
