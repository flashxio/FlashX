#ifndef __CRS_HEADER_H__
#define __CRS_HEADER_H__

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

class crs_header
{
	union {
		char page[4096];
		struct {
			size_t nrows;
			size_t ncols;
			size_t nnz;
		} s;
	} data;
public:
	crs_header() {
		memset(this, 0, sizeof(*this));
	}

	crs_header(size_t nrows, size_t ncols, size_t nnz) {
		data.s.nrows = nrows;
		data.s.ncols = ncols;
		data.s.nnz = nnz;
	}

	size_t get_num_rows() const {
		return data.s.nrows;
	}

	size_t get_num_cols() const {
		return data.s.ncols;
	}

	size_t get_num_non_zeros() const {
		return data.s.nnz;
	}
};

typedef size_t crs_idx_t;

#endif
