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

#include "io_request.h"

#include "vertex.h"
#include "sparse_matrix_format.h"

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
 * It's used for the matrix partitioned in one dimension.
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

class SpM_2d_index;

/*
 * This maps row blocks to different I/O generators. It also defines how
 * row blocks are distributed across threads.
 */
class row_block_mapper
{
public:
	/*
	 * A range of row blocks that are accessed together.
	 */
	struct rb_range
	{
		// The index of the row block in array.
		off_t idx;
		// The number of row blocks in the range.
		size_t num;
	};
private:
	std::vector<rb_range> ranges;

	void init(size_t num_rbs, size_t range_size);
public:
	row_block_mapper(const std::vector<row_block> &rblocks, size_t range_size) {
		// The last row block in the vector indicates the end of the row block
		// vector, so the true number of row blocks is one fewer.
		size_t num_rbs = rblocks.size() - 1;
		init(num_rbs, range_size);
	}
	row_block_mapper(const SpM_2d_index &index, size_t range_size);

	size_t get_num_ranges() const {
		return ranges.size();
	}

	rb_range get_range(off_t idx) const {
		return ranges[idx];
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
			size_t tot_num_cols, int file_id, const row_block_mapper &mapper);
	/*
	 * This generates an I/O generator that accesses a matrix partitioned
	 * into blocks.
	 */
	static matrix_io_generator::ptr create(SpM_2d_index::ptr idx, int file_id,
			size_t io_num_brows, size_t min_num_brows);

	// Get the next I/O access in the current worker thread.
	virtual matrix_io get_next_io() = 0;
	virtual bool has_next_io() = 0;
};

}

#endif
