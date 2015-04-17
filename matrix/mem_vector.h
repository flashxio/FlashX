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

#include "log.h"
#include "common.h"

#include "vector.h"
#include "raw_data_array.h"

namespace fm
{

class data_frame;
class scalar_type;

class mem_vector: public vector
{
	char *arr;
	detail::raw_data_array data;
	bool sorted;

	const char *get(off_t idx) const {
		return arr + idx * get_entry_size();
	}
protected:
	mem_vector(size_t length, const scalar_type &type);
	mem_vector(const detail::raw_data_array &data, size_t len, const scalar_type &type);
	mem_vector(const mem_vector &vec): vector(vec.get_length(), vec.get_type(),
			true) {
		this->arr = vec.arr;
		this->data = vec.data;
		this->sorted = vec.sorted;
	}
public:
	typedef std::shared_ptr<mem_vector> ptr;
	static ptr cast(vector::ptr vec);

	static ptr create(const detail::raw_data_array &data, size_t num_bytes,
			const scalar_type &type) {
		if (num_bytes % type.get_size() != 0) {
			BOOST_LOG_TRIVIAL(error)
				<< "The data array has a wrong number of bytes";
			return mem_vector::ptr();
		}
		size_t len = num_bytes / type.get_size();
		return ptr(new mem_vector(data, len, type));
	}

	static ptr create(size_t length, const scalar_type &type) {
		return ptr(new mem_vector(length, type));
	}

	char *get_raw_arr() {
		return arr;
	}

	const char *get_raw_arr() const {
		return arr;
	}

	virtual mem_vector::ptr get(const mem_vector &idxs) const;

	virtual bool equals(const mem_vector &vec) const;

	virtual bool is_sorted() const {
		return sorted;
	}

	bool verify_groupby(const gr_apply_operate<mem_vector> &op) const;
	std::shared_ptr<data_frame> serial_groupby(
			const gr_apply_operate<mem_vector> &op, bool with_val) const;
	using vector::groupby;
	virtual std::shared_ptr<data_frame> groupby(
			const gr_apply_operate<mem_vector> &op, bool with_val) const;
	virtual scalar_variable::ptr aggregate(const bulk_operate &op) const;

	virtual bool append(std::vector<vector::ptr>::const_iterator vec_it,
			std::vector<vector::ptr>::const_iterator vec_end);
	virtual bool append(const vector &vec);
	virtual bool resize(size_t new_length);
	virtual vector::ptr get_sub_vec(off_t start, size_t length) const;
	size_t get_sub_start() const;
	virtual bool expose_sub_vec(off_t start, size_t length);
	virtual vector::ptr deep_copy() const;

	virtual void reset_data() {
		data.reset_data();
	}

	virtual vector::ptr shallow_copy() const {
		return vector::ptr(new mem_vector(*this));
	}

	bool export2(FILE *f) const;

	void set_data(const set_operate &op);

	virtual vector::ptr sort_with_index();

	virtual void sort() {
		get_type().get_sorter().sort(arr, get_length(), false);
		sorted = true;
	}

	virtual void serial_sort() {
		get_type().get_sorter().serial_sort(arr, get_length(), false);
		sorted = true;
	}

	void set(const std::vector<const char *> &locs) {
		assert(locs.size() <= get_length());
		get_type().get_sg().gather(locs, arr);
	}

	template<class T>
	T get(off_t idx) const {
		const T *eles = (const T *) arr;
		return eles[idx];
	}

	template<class T>
	void set(off_t idx, T v) {
		T *eles = (T *) arr;
		eles[idx] = v;
	}

	template<class T>
	T max() const {
		const bulk_operate &max_op = *get_type().get_basic_ops().get_op(
				basic_ops::op_idx::MAX);
		scalar_variable::ptr res = aggregate(max_op);
		return *(T *) res->get_raw();
	}
};

template<class T>
class seq_set_operate: public type_set_operate<T>
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

	virtual void set(T *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		// We are initializing a single-column matrix.
		T v = from + row_idx * by;
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = v;
			v += by;
		}
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
	v->set_data(seq_set_operate<EntryType>(n, start, stride));
	return v;
}

template<>
vector::ptr create_vector<double>(double start, double end,
		double stride);

template<class EntryType>
class set_const_operate: public type_set_operate<EntryType>
{
	EntryType v;
public:
	set_const_operate(EntryType v) {
		this->v = v;
	}

	virtual void set(EntryType *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = v;
	}
};

/*
 * Create a vector filled with a constant value.
 */
template<class EntryType>
vector::ptr create_vector(size_t length, EntryType initv)
{
	// TODO let's just use in-memory dense matrix first.
	mem_vector::ptr v = mem_vector::create(length, get_scalar_type<EntryType>());
	v->set_data(set_const_operate<EntryType>(initv));
	return v;
}

}

#endif
