/**
 * Copyright 2014 Da Zheng
 *
 * This file is part of SA-GraphLib.
 *
 * SA-GraphLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SA-GraphLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SA-GraphLib.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __MY_BITMAP_H__
#define __MY_BITMAP_H__

#include <stdlib.h>

#include <vector>

#include "common.h"

/**
 * The functionality of this bitmap is very similar to std::vector<bool>
 * in STL. But this one is optimized for extra operations such as merge
 * and iterating on true values or on false values.
 */
class bitmap
{
	static const int NUM_BITS_LONG = sizeof(long) * 8;
	size_t max_num_bits;
	long *ptr;

	size_t get_num_longs() const {
		return ROUNDUP(max_num_bits, NUM_BITS_LONG) / NUM_BITS_LONG;
	}

	template<class T>
	static void get_set_bits_long(long value, size_t idx, std::vector<T> &v) {
		for (int i = 0; i < NUM_BITS_LONG; i++) {
			if (value & (1L << i))
				v.push_back(i + idx * NUM_BITS_LONG);
		}
	}
public:
#if 0
	class const_iterator
	{
		size_t idx;
		long *ptr;
	public:
		const_iterator(long *ptr) {
			this->ptr = ptr;
			idx = 0;
		}

		bool operator*() const {
		}

		const_iterator &operator++() {
		}

		bool operator==(const const_iterator &it) const {
			return this->idx == it.idx;
		}

		bool operator!=(const const_iterator &it) const {
			return this->idx != it.idx;
		}

		const_iterator &operator+=(int num) {
			idx += num;
			return *this;
		}
	};
#endif

	bitmap(size_t max_num_bits, int node_id) {
		this->max_num_bits = max_num_bits;
		size_t num_longs = get_num_longs();
		ptr = (long *) numa_alloc_onnode(num_longs * sizeof(ptr[0]), node_id);
		memset(ptr, 0, sizeof(ptr[0]) * num_longs);
	}

	~bitmap() {
		numa_free(ptr, get_num_longs() * sizeof(ptr[0]));
	}

	size_t get_num_bits() const {
		return max_num_bits;
	}

	void set(size_t idx) {
		assert(idx < max_num_bits);
		size_t arr_off = idx / NUM_BITS_LONG;
		size_t inside_off = idx % NUM_BITS_LONG;
		ptr[arr_off] |= (1L << inside_off);
	}

	bool get(size_t idx) const {
		assert(idx < max_num_bits);
		size_t arr_off = idx / NUM_BITS_LONG;
		size_t inside_off = idx % NUM_BITS_LONG;
		return ptr[arr_off] & (1L << inside_off);
	}

	void merge(const bitmap &map) {
		assert(this->get_num_bits() >= map.get_num_bits());
		size_t size = map.get_num_longs();
		for (size_t i = 0; i < size; i++) {
			ptr[i] |= map.ptr[i];
		}
	}

	void clear() {
		memset(ptr, 0, sizeof(ptr[0]) * get_num_longs());
	}

	/**
	 * This method collects all bits that have been set to 1.
	 */
	template<class T>
	size_t get_set_bits(std::vector<T> &v) const {
		size_t size = get_num_longs();
		for (size_t i = 0; i < size; i++) {
			if (ptr[i])
				get_set_bits_long(ptr[i], i, v);
		}
		return v.size();
	}
};

#endif
