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

#include "mem_vector.h"

namespace fm
{

mem_vector::mem_vector(mem_dense_matrix::ptr data): vector(
		// The length of the vector is the size of the dimension that isn't 1.
		data->get_num_rows() == 1 ? data->get_num_cols(): data->get_num_rows(),
		data->get_entry_size(), true)
{
	this->data = data;
	// The data buffer is the first row or column of the matrix, so
	// the shape of the matrix doesn't matter.
	// TODO this may not work with submatrix.
	if (data->store_layout() == matrix_layout_t::L_ROW)
		arr = mem_row_dense_matrix::cast(data)->get_row(0);
	else if (data->store_layout() == matrix_layout_t::L_COL)
		arr = mem_col_dense_matrix::cast(data)->get_col(0);
	else
		BOOST_LOG_TRIVIAL(error) << "wrong matrix layout";
}

mem_vector::mem_vector(size_t length, size_t entry_size): vector(length,
		entry_size, true)
{
	// Maybe the column form may be more useful.
	mem_col_dense_matrix::ptr tmp = mem_col_dense_matrix::create(length,
			1, entry_size);
	this->arr = tmp->get_col(0);
	this->data = std::static_pointer_cast<mem_dense_matrix>(tmp);
}

mem_vector::ptr mem_vector::cast(vector::ptr vec)
{
	if (!vec->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't cast a non-in-mem vector to in-mem vector";
		return mem_vector::ptr();
	}
	return std::static_pointer_cast<mem_vector>(vec);
}

mem_vector::const_ptr mem_vector::cast(vector::const_ptr vec)
{
	if (!vec->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't cast a non-in-mem vector to in-mem vector";
		return mem_vector::const_ptr();
	}
	return std::static_pointer_cast<const mem_vector>(vec);
}

bool mem_vector::append(std::vector<vector::ptr>::const_iterator vec_it,
		std::vector<vector::ptr>::const_iterator vec_end)
{
	// Get the total size of the result vector.
	size_t tot_res_size = this->get_length();
	for (auto it = vec_it; it != vec_end; it++) {
		tot_res_size += (*it)->get_length();
		if (!(*it)->is_in_mem()) {
			BOOST_LOG_TRIVIAL(error)
				<< "Not support appending an ext-mem vector to an in-mem vector";
			return false;
		}
	}

	// Merge all results to a single vector.
	off_t loc = this->get_length();
	this->resize(tot_res_size);
	for (auto it = vec_it; it != vec_end; it++) {
		assert(loc + (*it)->get_length() <= this->get_length());
		this->set_sub_vec(loc, **it);
		loc += (*it)->get_length();
	}
	return true;
}

bool mem_vector::resize(size_t new_length)
{
	// Keep the old information of the vector.
	mem_dense_matrix::ptr old_data = data;
	char *old_arr = arr;
	size_t old_length = get_length();

	// We realloate memory regardless of whether we increase or decrease
	// the length of the vector.
	// TODO we probably don't want to reallocate memory when shrinking
	// the vector size.
	mem_col_dense_matrix::ptr tmp = mem_col_dense_matrix::create(new_length,
			1, get_entry_size());
	if (tmp == NULL) {
		BOOST_LOG_TRIVIAL(error) << "can't allocate memory to resize the vector";
		return false;
	}

	this->arr = tmp->get_col(0);
	this->data = std::static_pointer_cast<mem_dense_matrix>(tmp);
	memcpy(arr, old_arr, std::min(old_length, new_length) * get_entry_size());
	return vector::resize(new_length);
}

bool mem_vector::set_sub_vec(off_t start, const vector &vec)
{
	if (!vec.is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "Not support setting a subvector from ext-mem vector";
		return false;
	}
	if (get_entry_size() != vec.get_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The two vectors don't have the same entry size";
		return false;
	}
	if (start + vec.get_length() > get_length()) {
		BOOST_LOG_TRIVIAL(error) << "set_sub_vec: out of range";
		return false;
	}

	const mem_vector &mem_vec = (const mem_vector &) vec;
	memcpy(get_raw_arr() + start * get_entry_size(), mem_vec.get_raw_arr(),
			mem_vec.get_length() * mem_vec.get_entry_size());
	return true;
}

vector::const_ptr mem_vector::get_sub_vec(off_t start, size_t length) const
{
	if (start + length > get_length()) {
		BOOST_LOG_TRIVIAL(error) << "get_sub_vec: out of range";
		return vector::ptr();
	}

	// We need to discard the const from the "this" pointer.
	mem_vector *mutable_this = (mem_vector *) this;

	mem_vector::ptr mem_vec = mem_vector::cast(mutable_this->shallow_copy());
	mem_vec->resize(length);
	mem_vec->arr = mutable_this->get_raw_arr() + start * get_entry_size();
	mem_vec->data = mutable_this->get_data();
	return std::static_pointer_cast<const vector>(mem_vec);
}

vector::ptr mem_vector::clone() const
{
	// We need to discard the const from the "this" pointer.
	mem_vector *mutable_this = (mem_vector *) this;
	mem_vector::ptr mem_vec = mem_vector::cast(mutable_this->shallow_copy());
	mem_vec->data = mem_dense_matrix::cast(data->clone());
	if (mem_vec->data->store_layout() == matrix_layout_t::L_ROW)
		mem_vec->arr = mem_row_dense_matrix::cast(mem_vec->data)->get_row(0);
	else if (mem_vec->data->store_layout() == matrix_layout_t::L_COL)
		mem_vec->arr = mem_col_dense_matrix::cast(mem_vec->data)->get_col(0);
	else
		BOOST_LOG_TRIVIAL(error) << "wrong matrix layout";
	return std::static_pointer_cast<vector>(mem_vec);
}

}
