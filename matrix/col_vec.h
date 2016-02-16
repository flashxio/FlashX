#ifndef __FM_COL_VEC_H__
#define __FM_COL_VEC_H__

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

#include "dense_matrix.h"

namespace fm
{

/*
 * This represents a vector with a one-col matrix.
 * As such, a vector can contain data that doesn't physically exist.
 */
class col_vec: public dense_matrix
{
protected:
	col_vec(detail::matrix_store::const_ptr mat): dense_matrix(mat) {
		assert(mat->get_num_cols() == 1);
	}
public:
	typedef std::shared_ptr<col_vec> ptr;
	typedef std::shared_ptr<const col_vec> const_ptr;

	template<class T>
	static ptr create_randn(size_t len) {
		dense_matrix::ptr mat = dense_matrix::create_randn<T>(0, 1, len, 1,
				matrix_layout_t::L_COL);
		return ptr(new col_vec(mat->get_raw_store()));
	}
	template<class T>
	static ptr create_randu(size_t len) {
		dense_matrix::ptr mat = dense_matrix::create_randu<T>(0, 1, len, 1,
				matrix_layout_t::L_COL);
		return ptr(new col_vec(mat->get_raw_store()));
	}

	static ptr create(dense_matrix::ptr mat);

	col_vec(): dense_matrix(NULL) {
	}

	col_vec(size_t len, const scalar_type &type): dense_matrix(len, 1,
			matrix_layout_t::L_COL, type) {
	}

	size_t get_length() const {
		return get_num_rows();
	}

	col_vec operator=(const dense_matrix &mat) {
		assert(mat.get_num_cols() == 1);
		assign(mat);
		return *this;
	}
};

}

#endif
