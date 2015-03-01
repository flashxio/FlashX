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
#if defined(_OPENMP)
#include <parallel/algorithm>
#else
#include <algorithm>
#endif

#include "mem_dense_matrix.h"
#include "vector.h"

namespace fm
{

class data_frame;

class mem_vector: public vector
{
	char *arr;
	mem_dense_matrix::ptr data;

protected:
	mem_vector(mem_dense_matrix::ptr data);
	mem_vector(size_t length, size_t entry_size);

	const char *get_raw_arr() const {
		return arr;
	}

	char *get_raw_arr() {
		return arr;
	}

	mem_dense_matrix::ptr get_data() const {
		return data;
	}

	virtual std::shared_ptr<mem_vector> create_int(size_t length) const = 0;
public:
	typedef std::shared_ptr<mem_vector> ptr;
	typedef std::shared_ptr<const mem_vector> const_ptr;
	static ptr cast(vector::ptr vec);
	static const_ptr cast(vector::const_ptr vec);

	mem_dense_matrix::ptr get_data() {
		return data;
	}

	char *get(off_t idx) {
		return arr + idx * get_entry_size();
	}

	const char *get(off_t idx) const {
		return arr + idx * get_entry_size();
	}

	bool verify_groupby(const agg_operate &find_next,
		const agg_operate &agg_op, const vec_creator &create) const;
	std::shared_ptr<data_frame> serial_groupby(const agg_operate &find_next,
		const agg_operate &agg_op, const vec_creator &create) const;
	virtual std::shared_ptr<data_frame> groupby(const agg_operate &find_next,
		const agg_operate &agg_op, const vec_creator &create) const;

	virtual bool append(std::vector<vector::ptr>::const_iterator vec_it,
			std::vector<vector::ptr>::const_iterator vec_end);
	virtual bool resize(size_t new_length);
	virtual bool set_sub_vec(off_t start, const vector &vec);
	virtual vector::const_ptr get_sub_vec(off_t start, size_t length) const;
	virtual vector::ptr clone() const;
};

/**
 * This vector implementation is a wrapper on a in-memory single-column
 * or single-row dense matrix.
 */
template<class T>
class type_mem_vector: public mem_vector
{
	bool sorted;

	type_mem_vector(mem_dense_matrix::ptr data): mem_vector(data) {
		sorted = false;
	}

	type_mem_vector(size_t length): mem_vector(length, sizeof(T)) {
		sorted = false;
	}

	virtual mem_vector::ptr create_int(size_t length) const {
		return mem_vector::ptr(new type_mem_vector<T>(length));
	}
public:
	typedef std::shared_ptr<type_mem_vector<T> > ptr;

	static ptr create(mem_dense_matrix::ptr data) {
		if (data->get_num_rows() > 1 && data->get_num_cols() > 1) {
			BOOST_LOG_TRIVIAL(error)
				<< "Can't convert a matrix with more than one row/column into a vector";
			return ptr();
		}
		else if (!data->is_type<T>()) {
			BOOST_LOG_TRIVIAL(error) << "The matrix has a wrong data type";
			return ptr();
		}
		return ptr(new type_mem_vector<T>(data));
	}

	static ptr create(size_t length) {
		return ptr(new type_mem_vector<T>(length));
	}

	static type_mem_vector<T>::ptr cast(vector::ptr vec) {
		if (!vec->is_in_mem())
			return type_mem_vector<T>::ptr();
		else if (!vec->is_type<T>())
			return type_mem_vector<T>::ptr();
		else
			return std::static_pointer_cast<type_mem_vector<T> >(vec);
	}

	T get(off_t idx) const {
		return ((T *) get_raw_arr())[idx];
	}

	void set(off_t idx, T v) {
		((T *) get_raw_arr())[idx] = v;
	}

	virtual bool is_sorted() const {
		return sorted;
	}

	virtual void sort() {
		T *start = (T *) get_raw_arr();
		T *end = start + get_length();
#if defined(_OPENMP)
		__gnu_parallel::sort(start, end);
#else
		std::sort(start, end);
#endif
		sorted = true;
	}

	virtual void serial_sort() {
		T *start = (T *) get_raw_arr();
		T *end = start + get_length();
		std::sort(start, end);
		sorted = true;
	}

	virtual vector::ptr shallow_copy() {
		type_mem_vector<T>::ptr ret = type_mem_vector<T>::create(get_data());
		ret->sorted = this->sorted;
		return std::static_pointer_cast<vector>(ret);
	}

	virtual vector::const_ptr shallow_copy() const {
		type_mem_vector<T>::ptr ret = type_mem_vector<T>::create(get_data());
		ret->sorted = this->sorted;
		return std::static_pointer_cast<const vector>(ret);
	}

	virtual prim_type get_type() const {
		return fm::get_type<T>();
	}
};

}

#endif
