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

namespace fm
{

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

/*
 * This represents a matrix block that we want to access. This data structure
 * maintains the mapping of the matrix block to a data region on disks.
 */
class matrix_io
{
	matrix_loc top_left;
	size_t num_rows;
	size_t num_cols;

	safs::data_loc_t loc;
	size_t size;
public:
	matrix_io(): top_left(-1, -1) {
		num_rows = 0;
		num_cols = 0;
		size = 0;
	}

	matrix_io(const matrix_loc &_top_left, size_t num_rows,
			size_t num_cols, safs::data_loc_t loc, size_t size): top_left(_top_left) {
		this->num_rows = num_rows;
		this->num_cols = num_cols;
		this->loc = loc;
		this->size = size;
	}

	const safs::data_loc_t &get_loc() const {
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

	bool is_valid() const {
		return loc.get_offset() >= 0;
	}
};

class sparse_matrix;

/*
 * This represents a minimal row block that we can access from disks.
 */
class row_block
{
	off_t off;
public:
	row_block(off_t off) {
		this->off = off;
	}

	off_t get_offset() const {
		return off;
	}
};

/*
 * This class generates I/O accesses for a worker thread.
 * It defines which matrix blocks are accessed by the current worker thread.
 * An I/O generator is referenced by multiple threads, so all of its methods
 * need to be thread-safe.
 */
class matrix_io_generator
{
public:
	typedef std::shared_ptr<matrix_io_generator> ptr;

	/*
	 * This generates an I/O generator that accesses a matrix by rows.
	 */
	static matrix_io_generator::ptr create(
			const std::vector<row_block> &_blocks, size_t tot_num_rows,
			size_t tot_num_cols, int file_id, int gen_id, int num_gens);

	// Get the next I/O access in the current worker thread.
	virtual matrix_io get_next_io() = 0;
	// Get an I/O access from a worker thread that doesn't own the I/O generator.
	// The I/O size returned from this method is typically much smaller than
	// the one returned by get_next_io().
	virtual matrix_io steal_io() = 0;
	virtual bool has_next_io() = 0;
};

}

#endif
