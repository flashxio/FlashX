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

#include "factor.h"

namespace fm
{

factor_col_vector::ptr factor_col_vector::create(const factor &f,
		dense_matrix::ptr mat)
{
	if (mat->get_num_cols() > 1) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't convert a matrix with more than one col to a vector";
		return ptr();
	}
	if (mat->get_data().store_layout() == matrix_layout_t::L_ROW)
		mat = mat->conv2(matrix_layout_t::L_COL);
	return ptr(new factor_col_vector(f, mat->get_raw_store()));
}

static detail::matrix_store::ptr create_factor_store(size_t len,
		int num_nodes, bool in_mem, const set_operate &op)
{
	detail::matrix_store::ptr ret = detail::matrix_store::create(len, 1,
			matrix_layout_t::L_COL, get_scalar_type<factor_value_t>(),
			num_nodes, in_mem);
	ret->set_data(op);
	return ret;
}

factor_col_vector::factor_col_vector(const factor &_f, size_t len,
		int num_nodes, bool in_mem, const set_operate &op): col_vec(
			create_factor_store(len, num_nodes, in_mem, op)), f(_f) {
}

}
