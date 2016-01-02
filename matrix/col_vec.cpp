/*
 * Copyright 2016 Open Connectome Project (http://openconnecto.me)
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

#include "col_vec.h"

namespace fm
{

col_vec::ptr col_vec::create(dense_matrix::ptr mat)
{
	if (mat->get_num_cols() > 1) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't convert a matrix with more than one col to a vector";
		return ptr();
	}
	if (mat->get_data().store_layout() == matrix_layout_t::L_ROW)
		mat = mat->conv2(matrix_layout_t::L_COL);
	return ptr(new col_vec(mat->get_raw_store()));
}

}
