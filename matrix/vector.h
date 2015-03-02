#ifndef __FM_VECTOR_H__
#define __FM_VECTOR_H__

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

#include "generic_type.h"

namespace fm
{

class bulk_operate;
class vec_creator;
class data_frame;
class agg_operate;

class vector
{
	size_t length;
	size_t entry_size;
	bool in_mem;

protected:
	vector(size_t length, size_t entry_size, bool in_mem) {
		this->length = length;
		this->entry_size = entry_size;
		this->in_mem = in_mem;
	}
public:
	typedef std::shared_ptr<vector> ptr;
	typedef std::shared_ptr<const vector> const_ptr;

	virtual ~vector() {
	}

	bool is_in_mem() const {
		return in_mem;
	}

	size_t get_length() const {
		return length;
	}

	template<class T>
	bool is_type() const {
		return this->get_type() == fm::get_type<T>();
	}

	size_t get_entry_size() const {
		return entry_size;
	}

	virtual bool resize(size_t new_length) {
		this->length = new_length;
		return true;
	}

	virtual prim_type get_type() const = 0;
	virtual bool set_sub_vec(off_t start, const vector &vec) = 0;
	virtual vector::const_ptr get_sub_vec(off_t start, size_t length) const = 0;
	virtual bool append(std::vector<vector::ptr>::const_iterator vec_it,
			std::vector<vector::ptr>::const_iterator vec_end) = 0;

	virtual void sort() = 0;
	virtual vector::ptr sort_with_index() = 0;
	virtual bool is_sorted() const = 0;
	// It should return data frame instead of vector.
	virtual std::shared_ptr<data_frame> groupby(const agg_operate &find_next,
		const agg_operate &agg_op, const vec_creator &create) const = 0;

	/**
	 * This method copies all members of the vector object as well as
	 * the data in the vector;
	 */
	virtual vector::ptr deep_copy() const = 0;
	/**
	 * This only copies all members of the vector object.
	 */
	virtual vector::ptr shallow_copy() = 0;
	virtual vector::const_ptr shallow_copy() const = 0;
};

/*
 * This interface allows us to create a vector with type.
 */
class vec_creator
{
public:
	virtual size_t get_entry_size() const = 0;
	virtual vector::ptr create(size_t length) const = 0;
};

}

#endif
