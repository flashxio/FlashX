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

#include "matrix_header.h"
#include "matrix_exception.h"

namespace fm
{

void matrix_header::init(matrix_type mat_type, size_t entry_size, size_t nrows,
		size_t ncols, matrix_layout_t layout, prim_type data_type)
{
	memset(u.page, 0, sizeof(u.page));

	u.d.magic_number = MAGIC_NUMBER;
	u.d.version_number = CURR_VERSION;
	u.d.type = mat_type;
	u.d.entry_size = entry_size;
	u.d.nrows = nrows;
	u.d.ncols = ncols;
	u.d.layout = layout;
	u.d.data_type = data_type;
}

matrix_header::matrix_header(matrix_type mat_type, size_t entry_size,
		size_t nrows, size_t ncols, matrix_layout_t layout,
		prim_type data_type, const block_2d_size &block_size)
{
	init(mat_type, entry_size, nrows, ncols, layout, data_type);
	u.d.block_2d_height = block_size.get_num_rows();
	u.d.block_2d_width = block_size.get_num_cols();
}

void matrix_header::verify() const
{
	if (!is_matrix_file())
		throw wrong_format(boost::str(boost::format("wrong magic number: %1%")
					% u.d.magic_number));
	if (!is_right_version())
		throw wrong_format(boost::str(boost::format("wrong version number: %1%")
					% u.d.version_number));
}

}
