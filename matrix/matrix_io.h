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

class matrix_io_generator
{
public:
	typedef std::shared_ptr<matrix_io_generator> ptr;

	static matrix_io_generator::ptr create(
			const std::vector<row_block> &_blocks, size_t tot_num_rows,
			size_t tot_num_cols, int file_id, int gen_id, int num_gens);

	virtual matrix_io get_next_io() = 0;
	virtual bool has_next_io() = 0;
};

#endif
