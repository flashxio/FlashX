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

#include <boost/format.hpp>

#include "mem_dense_matrix.h"
#include "vector.h"

namespace fm
{

class data_frame;

template<class T> class type_mem_vector;

class mem_vector: public vector
{
	char *arr;
	mem_dense_matrix::ptr data;

protected:
	mem_vector(mem_dense_matrix::ptr data);
	mem_vector(size_t length, size_t entry_size);
	mem_vector(std::shared_ptr<char> data, size_t len, size_t entry_size);

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

	char *get_raw_arr() {
		return arr;
	}

	const char *get_raw_arr() const {
		return arr;
	}

	char *get(off_t idx) {
		return arr + idx * get_entry_size();
	}

	const char *get(off_t idx) const {
		return arr + idx * get_entry_size();
	}
	virtual mem_vector::ptr get(type_mem_vector<off_t> &idxs) const;

	virtual bool equals(const mem_vector &vec) const;

	bool verify_groupby(const gr_apply_operate<mem_vector> &op) const;
	std::shared_ptr<data_frame> serial_groupby(
			const gr_apply_operate<mem_vector> &op, bool with_val) const;
	using vector::groupby;
	virtual std::shared_ptr<data_frame> groupby(
			const gr_apply_operate<mem_vector> &op, bool with_val) const;

	virtual bool append(std::vector<vector::ptr>::const_iterator vec_it,
			std::vector<vector::ptr>::const_iterator vec_end);
	virtual bool append(const vector &vec);
	virtual bool resize(size_t new_length);
	virtual bool set_sub_vec(off_t start, const vector &vec);
	virtual vector::const_ptr get_sub_vec(off_t start, size_t length) const;
	virtual vector::ptr get_sub_vec(off_t start, size_t length);
	size_t get_sub_start() const;
	virtual bool expose_sub_vec(off_t start, size_t length);
	virtual vector::ptr deep_copy() const;

	bool export2(FILE *f) const;
};

/**
 * This vector implementation is a wrapper on a in-memory single-column
 * or single-row dense matrix.
 */
template<class T>
class type_mem_vector: public mem_vector
{
	bool sorted;

	virtual mem_vector::ptr create_int(size_t length) const {
		return mem_vector::ptr(new type_mem_vector<T>(length));
	}

protected:
	type_mem_vector(std::shared_ptr<char> data, size_t len): mem_vector(
			data, len, sizeof(T)) {
		sorted = false;
	}

	type_mem_vector(mem_dense_matrix::ptr data): mem_vector(data) {
		sorted = false;
	}

	type_mem_vector(size_t length): mem_vector(length, sizeof(T)) {
		sorted = false;
	}
public:
	typedef std::shared_ptr<type_mem_vector<T> > ptr;

	static ptr create(std::shared_ptr<char> data, size_t num_bytes) {
		if (num_bytes % sizeof(T) != 0) {
			BOOST_LOG_TRIVIAL(error)
				<< "The data array has a wrong number of bytes";
			return type_mem_vector<T>::ptr();
		}
		size_t len = num_bytes / sizeof(T);
		return ptr(new type_mem_vector(data, len));
	}

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

	static ptr cast(vector::ptr vec) {
		if (!vec->is_in_mem()) {
			BOOST_LOG_TRIVIAL(error) << "the vector isn't in memory";
			return ptr();
		}
		if (!vec->is_type<T>()) {
			BOOST_LOG_TRIVIAL(error) << "the vector isn't of the right type";
			return ptr();
		}
		return std::static_pointer_cast<type_mem_vector<T> >(vec);
	}

	T get(off_t idx) const {
		return ((T *) get_raw_arr())[idx];
	}

	bool equals(const mem_vector &vec) const {
		const type_mem_vector<T> &t_vec = (const type_mem_vector<T> &) vec;
		for (size_t i = 0; i < t_vec.get_length(); i++)
			if (t_vec.get(i) != this->get(i))
				return false;
		return true;
	}

	void set(off_t idx, T v) {
		((T *) get_raw_arr())[idx] = v;
	}

	virtual void set_data(const type_set_operate<T> &op) {
		get_data()->set_data(op);
		T *start = (T *) get_raw_arr();
		sorted = std::is_sorted(start, start + get_length());
	}

	virtual bool is_sorted() const {
		return sorted;
	}

	virtual vector::ptr sort_with_index() {
		struct indexed_entry {
			T val;
			off_t idx;
			bool operator<(const indexed_entry &e) const {
				return val < e.val;
			}
		};
		std::unique_ptr<indexed_entry[]> entries
			= std::unique_ptr<indexed_entry[]>(new indexed_entry[get_length()]);
		for (size_t i = 0; i < get_length(); i++) {
			entries[i].val = get(i);
			entries[i].idx = i;
		}

		indexed_entry *start = entries.get();
		indexed_entry *end = start + get_length();
#if defined(_OPENMP)
		__gnu_parallel::sort(start, end);
#else
		std::sort(start, end);
#endif
		type_mem_vector<off_t>::ptr indexes
			= type_mem_vector<off_t>::create(get_length());
		for (size_t i = 0; i < get_length(); i++) {
			set(i, start[i].val);
			indexes->set(i, start[i].idx);
		}
		sorted = true;
		return indexes;
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
		// The vector might be a sub vector.
		off_t start = get_sub_start();
		if (start == 0)
			return std::static_pointer_cast<vector>(ret);
		else
			return ret->get_sub_vec(start, get_length());
	}

	virtual vector::const_ptr shallow_copy() const {
		type_mem_vector<T>::ptr ret = type_mem_vector<T>::create(get_data());
		ret->sorted = this->sorted;
		// The vector might be a sub vector.
		off_t start = get_sub_start();
		if (start == 0)
			return std::static_pointer_cast<const vector>(ret);
		else
			return std::static_pointer_cast<const vector>(
					ret->get_sub_vec(start, get_length()));
	}

	virtual const scalar_type &get_type() const {
		return get_scalar_type<T>();
	}

	T max() const {
		static basic_ops_impl<T, T, T> ops;
		scalar_variable_impl<T> res;
		get_data()->aggregate(*ops.get_op(basic_ops::op_idx::MAX), res);
		return res.get();
	}
};

template<class T>
class seq_set_operate: public set_operate
{
	long n;
	T from;
	T by;
public:
	seq_set_operate(long n, T from, T by) {
		this->n = n;
		this->from = from;
		this->by = by;
	}

	virtual void set(void *raw_arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		T *arr = (T *) raw_arr;
		// We are initializing a single-column matrix.
		T v = from + row_idx * by;
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = v;
			v += by;
		}
	}

	virtual size_t entry_size() const {
		return sizeof(T);
	}
};

/*
 * Create a sequence of values in [start, end]. `end' is inclusive.
 */
template<class EntryType>
vector::ptr create_vector(EntryType start, EntryType end, EntryType stride)
{
	if ((end < start && stride > 0) || (end > stride && stride < 0)) {
		BOOST_LOG_TRIVIAL(error)
			<< "There are a negative number of elements in the sequence";
		return vector::ptr();
	}
	// TODO let's just use in-memory dense matrix first.
	long n = (end - start) / stride;
	// We need to count the start element.
	n++;
	typename type_mem_vector<EntryType>::ptr v
		= type_mem_vector<EntryType>::create(n);
	v->get_data()->set_data(seq_set_operate<EntryType>(n, start, stride));
	return std::static_pointer_cast<vector>(v);
}

template<>
vector::ptr create_vector<double>(double start, double end,
		double stride);

template<class EntryType>
class set_const_operate: public set_operate
{
	EntryType v;
public:
	set_const_operate(EntryType v) {
		this->v = v;
	}

	virtual void set(void *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		EntryType *ele_p = (EntryType *) arr;
		for (size_t i = 0; i < num_eles; i++)
			ele_p[i] = v;
	}

	virtual size_t entry_size() const {
		return sizeof(EntryType);
	}
};

template<class EntryType>
vector::ptr create_vector(size_t length, EntryType initv)
{
	// TODO let's just use in-memory dense matrix first.
	typename type_mem_vector<EntryType>::ptr v
		= type_mem_vector<EntryType>::create(length);
	v->get_data()->set_data(set_const_operate<EntryType>(initv));
	return std::static_pointer_cast<vector>(v);
}

}

#endif
