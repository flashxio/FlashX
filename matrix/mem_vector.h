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

#include <boost/format.hpp>

#include "mem_dense_matrix.h"
#include "vector.h"

namespace fm
{

class data_frame;
class scalar_type;

class mem_vector: public vector
{
	char *arr;
	mem_dense_matrix::ptr data;
	bool sorted;

protected:
	mem_vector(mem_dense_matrix::ptr data);
	mem_vector(size_t length, const scalar_type &type);
	mem_vector(std::shared_ptr<char> data, size_t len, const scalar_type &type);
	mem_vector(const mem_vector &vec): vector(vec.get_length(), vec.get_entry_size(),
			true) {
		this->arr = vec.arr;
		this->data = vec.data;
		this->sorted = vec.sorted;
	}

	mem_dense_matrix::ptr get_data() const {
		return data;
	}
public:
	typedef std::shared_ptr<mem_vector> ptr;
	typedef std::shared_ptr<const mem_vector> const_ptr;
	static ptr cast(vector::ptr vec);
	static const_ptr cast(vector::const_ptr vec);

	static ptr create(std::shared_ptr<char> data, size_t num_bytes,
			const scalar_type &type) {
		if (num_bytes % type.get_size() != 0) {
			BOOST_LOG_TRIVIAL(error)
				<< "The data array has a wrong number of bytes";
			return mem_vector::ptr();
		}
		size_t len = num_bytes / type.get_size();
		return ptr(new mem_vector(data, len, type));
	}

	static ptr create(mem_dense_matrix::ptr data) {
		if (data->get_num_rows() > 1 && data->get_num_cols() > 1) {
			BOOST_LOG_TRIVIAL(error)
				<< "Can't convert a matrix with more than one row/column into a vector";
			return ptr();
		}
		return ptr(new mem_vector(data));
	}

	static ptr create(size_t length, const scalar_type &type) {
		return ptr(new mem_vector(length, type));
	}

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
	virtual mem_vector::ptr get(const mem_vector &idxs) const;

	virtual bool equals(const mem_vector &vec) const;
	virtual const scalar_type &get_type() const {
		return data->get_type();
	}

	virtual bool is_sorted() const {
		return sorted;
	}

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

	virtual vector::ptr shallow_copy() {
		return vector::ptr(new mem_vector(*this));
	}

	virtual vector::const_ptr shallow_copy() const {
		return vector::const_ptr(new mem_vector(*this));
	}

	bool export2(FILE *f) const;

	void set_data(const set_operate &op) {
		get_data()->set_data(op);
		sorted = get_type().get_sorter().is_sorted(get_raw_arr(),
				get_length(), false);
	}

	virtual vector::ptr sort_with_index();

	virtual void sort() {
		get_type().get_sorter().sort(get_raw_arr(), get_length(), false);
		sorted = true;
	}

	virtual void serial_sort() {
		get_type().get_sorter().serial_sort(get_raw_arr(), get_length(), false);
		sorted = true;
	}

	template<class T>
	T get(off_t idx) const {
		return *(const T*) get(idx);
	}

	template<class T>
	void set(off_t idx, T v) {
		*(T*) get(idx) = v;
	}

	template<class T>
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
	mem_vector::ptr v = mem_vector::create(n, get_scalar_type<EntryType>());
	v->get_data()->set_data(seq_set_operate<EntryType>(n, start, stride));
	return v;
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
	mem_vector::ptr v = mem_vector::create(length, get_scalar_type<EntryType>());
	v->get_data()->set_data(set_const_operate<EntryType>(initv));
	return v;
}

}

#endif
