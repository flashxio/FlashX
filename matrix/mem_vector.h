#ifndef __MEM_VECTOR_H__
#define __MEM_VECTOR_H__

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

#include <memory>

#include "mem_dense_matrix.h"

namespace fm
{

/**
 * This vector implementation is a wrapper on a in-memory single-column
 * or single-row dense matrix.
 */
template<class T>
class mem_vector
{
	T *arr;
	size_t length;
	mem_dense_matrix::ptr data;

	mem_vector(mem_dense_matrix::ptr data) {
		this->data = data;
		// The length of the vector is the size of the dimension that isn't 1.
		if (data->get_num_rows() == 1)
			length = data->get_num_cols();
		else
			length = data->get_num_rows();
		// The data buffer is the first row or column of the matrix, so
		// the shape of the matrix doesn't matter.
		// TODO this may not work with submatrix.
		if (data->store_layout() == matrix_layout_t::L_ROW)
			arr = (T *) mem_row_dense_matrix::cast(data)->get_row(0);
		else if (data->store_layout() == matrix_layout_t::L_COL)
			arr = (T *) mem_col_dense_matrix::cast(data)->get_col(0);
		else
			BOOST_LOG_TRIVIAL(error) << "wrong matrix layout";
	}

	mem_vector(size_t length) {
		// Maybe the column form may be more useful.
		mem_col_dense_matrix::ptr tmp = mem_col_dense_matrix::create(length,
				1, sizeof(T));
		this->arr = (T *) tmp->get_col(0);
		this->data = std::static_pointer_cast<mem_dense_matrix>(tmp);
		this->length = length;
	}
public:
	typedef std::shared_ptr<mem_vector<T> > ptr;

	static ptr create(mem_dense_matrix::ptr data) {
		if (data->get_num_rows() > 1 && data->get_num_cols() > 1) {
			BOOST_LOG_TRIVIAL(error)
				<< "Can't convert a matrix with more than one row/column into a vector";
			return ptr();
		}
		else if (data->get_entry_size() != sizeof(T)) {
			BOOST_LOG_TRIVIAL(error) << "The matrix has a wrong data type";
			return ptr();
		}
		return ptr(new mem_vector(data));
	}

	static ptr create(size_t length) {
		return ptr(new mem_vector(length));
	}

	T get(off_t idx) const {
		return arr[idx];
	}

	void set(off_t idx, T v) {
		arr[idx] = v;
	}

	size_t get_length() const {
		return length;
	}

	mem_dense_matrix::ptr get_data() {
		return data;
	}
};

}

#endif
