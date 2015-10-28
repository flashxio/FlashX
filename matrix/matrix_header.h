#ifndef __MATRIX_HEADER_H__
#define __MATRIX_HEADER_H__

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
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "generic_type.h"

namespace fm
{

enum matrix_type
{
	VECTOR,
	DENSE,
	SPARSE,
};

enum matrix_layout_t
{
	L_COL,
	L_ROW,
	L_ROW_2D,
	// It indicates that the layout isn't defined.
	L_NONE,
};

static const size_t block_max_num_rows = ((size_t) std::numeric_limits<int16_t>::max()) + 1;
static const size_t block_max_num_cols = ((size_t) std::numeric_limits<int16_t>::max()) + 1;

class block_2d_size
{
	uint16_t nrow_log;
	uint16_t ncol_log;
public:
	block_2d_size() {
		nrow_log = 0;
		ncol_log = 0;
	}

	block_2d_size(size_t num_rows, size_t num_cols) {
		assert(num_rows <= block_max_num_rows);
		assert(num_cols <= block_max_num_cols);
		nrow_log = log2(num_rows);
		ncol_log = log2(num_cols);
		assert((1U << nrow_log) == num_rows);
		assert((1U << ncol_log) == num_cols);
	}

	size_t get_nrow_log() const {
		return nrow_log;
	}

	size_t get_ncol_log() const {
		return ncol_log;
	}

	size_t get_nrow_mask() const {
		return get_num_rows() - 1;
	}

	size_t get_ncol_mask() const {
		return get_num_cols() - 1;
	}

	size_t get_num_rows() const {
		return 1 << nrow_log;
	}

	size_t get_num_cols() const {
		return 1 << ncol_log;
	}

	size_t cal_num_block_rows(size_t tot_num_rows) const {
		return ceil(((double) tot_num_rows) / get_num_rows());
	}
};

/*
 * The matrix header contains the basic information about the matrix
 * when the matrix is stored on disks.
 */
class matrix_header
{
	static const int64_t MAGIC_NUMBER = 0x123456789FFFFFFL;
	static const int CURR_VERSION = 2;

	union {
		struct info {
			int64_t magic_number;
			int version_number;
			matrix_type type;
			size_t entry_size;
			size_t nrows;
			size_t ncols;
			matrix_layout_t layout;
			prim_type data_type;
			// These two fields are valid when the matrix layout is
			// 2D partitioned.
			size_t block_2d_height;
			size_t block_2d_width;
			// It doesn't include the entry type, which should be determined by
			// users at runtime.
		} d;

		char page[4096];
	} u;

	void init(matrix_type mat_type, size_t entry_size, size_t nrows,
			size_t ncols, matrix_layout_t layout, prim_type data_type);
public:
	static size_t get_header_size() {
		return sizeof(matrix_header);
	}

	matrix_header() {
		memset(u.page, 0, sizeof(u.page));
	}

	matrix_header(matrix_type mat_type, size_t entry_size, size_t nrows,
			size_t ncols, matrix_layout_t layout, prim_type data_type) {
		init(mat_type, entry_size, nrows, ncols, layout, data_type);
	}

	matrix_header(matrix_type mat_type, size_t entry_size, size_t nrows,
			size_t ncols, matrix_layout_t layout, prim_type data_type,
			const block_2d_size &block_size);

	matrix_type get_matrix_type() const {
		return u.d.type;
	}

	bool is_sparse() const {
		return get_matrix_type() == matrix_type::SPARSE;
	}

	size_t get_entry_size() const {
		return u.d.entry_size;
	}

	size_t get_num_rows() const {
		return u.d.nrows;
	}

	size_t get_num_cols() const {
		return u.d.ncols;
	}

	matrix_layout_t get_layout() const {
		return u.d.layout;
	}

	const scalar_type &get_data_type() const {
		return get_scalar_type(u.d.data_type);
	}

	bool is_matrix_file() const {
		return u.d.magic_number == MAGIC_NUMBER;
	}

	bool is_right_version() const {
		return u.d.version_number == CURR_VERSION;
	}

	void verify() const;

	block_2d_size get_2d_block_size() const {
		return block_2d_size(u.d.block_2d_height, u.d.block_2d_width);
	}
};

}

#endif
