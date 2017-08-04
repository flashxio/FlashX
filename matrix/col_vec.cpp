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
#include "vector.h"
#include "mem_worker_thread.h"
#include "generic_hashtable.h"

namespace fm
{

col_vec::ptr col_vec::create(detail::matrix_store::const_ptr store)
{
	if (store->get_num_cols() > 1 && store->get_num_rows() > 1) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't convert a matrix store with multiple cols&rows to a vector";
		assert(0);
		return ptr();
	}
	dense_matrix::ptr mat = dense_matrix::create(store);
	if (mat->get_num_cols() > 1)
		mat = mat->transpose();
	if (mat->get_data().store_layout() == matrix_layout_t::L_ROW)
		mat = mat->conv2(matrix_layout_t::L_COL);
	return ptr(new col_vec(mat->get_raw_store()));
}

col_vec::ptr col_vec::create(detail::vec_store::const_ptr store)
{
	auto mat = store->conv2mat(store->get_length(), 1, false);
	return col_vec::create(mat);
}

col_vec::ptr col_vec::create(dense_matrix::ptr mat)
{
	if (mat->get_num_cols() > 1 && mat->get_num_rows() > 1) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't convert a matrix with multiple cols&rows to a vector";
		return ptr();
	}
	if (mat->get_num_cols() > 1)
		mat = mat->transpose();

	if (mat->get_data().store_layout() == matrix_layout_t::L_ROW)
		mat = mat->conv2(matrix_layout_t::L_COL);
	return ptr(new col_vec(mat->get_raw_store()));
}

col_vec::ptr col_vec::create(vector::const_ptr vec)
{
	dense_matrix::ptr mat = vec->conv2mat(vec->get_length(), 1, false);
	assert(mat->store_layout() == matrix_layout_t::L_COL);
	if (mat == NULL)
		return col_vec::ptr();
	else
		return col_vec::create(mat);
}

bool col_vec::is_seq() const
{
	// TODO
	return true;
}

}
