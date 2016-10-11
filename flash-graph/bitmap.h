/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashGraph.
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

#ifndef __MY_BITMAP_H__
#define __MY_BITMAP_H__

#include <stdlib.h>

#ifdef USE_NUMA
#include <numa.h>
#endif

#include <vector>
#include <atomic>

#include "common.h"

static const int NUM_BITS_LONG = sizeof(long) * 8;

/*
 * The functionality of this bitmap is very similar to std::vector<bool>
 * in STL. But this one is optimized for extra operations such as merge
 * and iterating on true values or on false values.
 */
class bitmap
{
	long num_set_bits;
	size_t max_num_bits;
	long *ptr;

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
		this->num_set_bits = 0;
		if (max_num_bits > 0) {
			size_t num_longs = get_num_longs();
			ptr = (long *) malloc_large(num_longs * sizeof(ptr[0]));
			memset(ptr, 0, sizeof(ptr[0]) * num_longs);
		}
		else
			ptr = NULL;
	}

	~bitmap() {
		if (ptr)
			free_large(ptr, get_num_longs() * sizeof(ptr[0]));
	}

	size_t get_num_longs() const {
		return ROUNDUP(max_num_bits, NUM_BITS_LONG) / NUM_BITS_LONG;
	}

	size_t get_num_bits() const {
		return max_num_bits;
	}

	size_t get_num_set_bits() const {
		return num_set_bits;
	}

	void set_all() {
		if (max_num_bits > 0) {
			memset(ptr, 0xff, sizeof(ptr[0]) * (get_num_longs() - 1));
			num_set_bits = NUM_BITS_LONG * (get_num_longs() - 1);
			for (size_t i = num_set_bits; i < max_num_bits; i++)
				set(i);
		}
	}

	void reset(size_t idx) {
		assert(idx < max_num_bits);
		size_t arr_off = idx / NUM_BITS_LONG;
		size_t inside_off = idx % NUM_BITS_LONG;
		// If the bit has beeen set, we need to decrease the count.
		if ((ptr[arr_off] & (1L << inside_off))) {
			num_set_bits--;
			ptr[arr_off] &= ~(1L << inside_off);
		}
	}

	void set(size_t idx) {
		assert(idx < max_num_bits);
		size_t arr_off = idx / NUM_BITS_LONG;
		size_t inside_off = idx % NUM_BITS_LONG;
		// If the bit hasn't been set, we now need to increase the count.
		if (!(ptr[arr_off] & (1L << inside_off))) {
			num_set_bits++;
			ptr[arr_off] |= (1L << inside_off);
		}
	}

	bool get(size_t idx) const {
		assert(idx < max_num_bits);
		size_t arr_off = idx / NUM_BITS_LONG;
		size_t inside_off = idx % NUM_BITS_LONG;
		return ptr[arr_off] & (1L << inside_off);
	}

	void clear() {
		memset(ptr, 0, sizeof(ptr[0]) * get_num_longs());
		num_set_bits = 0;
	}

	/*
	 * This method collects all bits that have been set to 1.
	 */
	template<class T>
	size_t get_set_bits(std::vector<T> &v) const {
		size_t size = get_num_longs();
		for (size_t i = 0; i < size; i++) {
			if (ptr[i])
				get_set_bits_long(ptr[i], i, v);
		}
		assert(v.size() == (size_t) num_set_bits);
		return v.size();
	}

	template<class T>
	size_t get_reset_set_bits(std::vector<T> &v) {
		size_t size = get_num_longs();
		for (size_t i = 0; i < size; i++) {
			if (ptr[i]) {
				get_set_bits_long(ptr[i], i, v);
				ptr[i] = 0;
			}
		}
		assert(v.size() == (size_t) num_set_bits);
		num_set_bits = 0;
		return v.size();
	}

	/*
	 * This method collects a specified number of bits that have been
	 * set to 1.
	 */
	template<class T>
	size_t get_set_bits(size_t begin_idx, size_t end_idx, std::vector<T> &v) const {
		// For simplicity, begin_index has to referent to the beginning
		// of a long.
		assert(begin_idx % NUM_BITS_LONG == 0);
		if (end_idx == get_num_bits())
			end_idx = ROUNDUP(end_idx, NUM_BITS_LONG);
		assert(end_idx % NUM_BITS_LONG == 0);
		// Find the last long we should visit (excluded).
		size_t long_end = end_idx / NUM_BITS_LONG;
		if (long_end > get_num_longs())
			long_end = get_num_longs();
		size_t orig_size = v.size();
		for (size_t i = begin_idx / NUM_BITS_LONG; i < long_end; i++) {
			if (ptr[i])
				get_set_bits_long(ptr[i], i, v);
		}
		return v.size() - orig_size;
	}

	template<class T>
	size_t get_reset_set_bits(size_t begin_idx, size_t end_idx, std::vector<T> &v) {
		// For simplicity, begin_index has to referent to the beginning
		// of a long.
		assert(begin_idx % NUM_BITS_LONG == 0);
		if (end_idx == get_num_bits())
			end_idx = ROUNDUP(end_idx, NUM_BITS_LONG);
		assert(end_idx % NUM_BITS_LONG == 0);
		// Find the last long we should visit (excluded).
		size_t long_end = end_idx / NUM_BITS_LONG;
		if (long_end > get_num_longs())
			long_end = get_num_longs();
		size_t orig_size = v.size();
		for (size_t i = begin_idx / NUM_BITS_LONG; i < long_end; i++) {
			if (ptr[i]) {
				get_set_bits_long(ptr[i], i, v);
				ptr[i] = 0;
			}
		}
		num_set_bits -= (v.size() - orig_size);
		assert(num_set_bits >= 0);
		return v.size() - orig_size;
	}

	void copy_to(bitmap &map) const {
		assert(max_num_bits == map.max_num_bits);
		map.num_set_bits = num_set_bits;
		memcpy(map.ptr, ptr, map.get_num_longs() * sizeof(long));
	}
};

/*
 * This is a thread-safe bitmap.
 * All set/clear operations on the bitmap is atomic. However, users of
 * the bitmap need to insert memory barrier themselves.
 */
class thread_safe_bitmap
{
	size_t max_num_bits;
	std::atomic_ulong *ptr;
public:
	thread_safe_bitmap(size_t max_num_bits, int node_id) {
		this->max_num_bits = max_num_bits;
		size_t num_longs = ROUNDUP(max_num_bits, NUM_BITS_LONG) / NUM_BITS_LONG;
#ifdef USE_NUMA
		ptr = (std::atomic_ulong *) numa_alloc_onnode(
				num_longs * sizeof(*ptr), node_id);
#else
		ptr = (std::atomic_ulong *) malloc_aligned(num_longs * sizeof(*ptr), 4096);
#endif
		clear();
	}

	~thread_safe_bitmap() {
		size_t num_longs = ROUNDUP(max_num_bits, NUM_BITS_LONG) / NUM_BITS_LONG;
#ifdef USE_NUMA
		numa_free(ptr, num_longs * sizeof(*ptr));
#else
		free(ptr);
#endif
	}

	size_t get_num_bits() const {
		return max_num_bits;
	}

	void set(size_t idx) {
		assert(idx < max_num_bits);
		size_t arr_off = idx / NUM_BITS_LONG;
		size_t inside_off = idx % NUM_BITS_LONG;
		// We only want atomicity here.
		ptr[arr_off].fetch_or(1L << inside_off, std::memory_order_relaxed);
	}

	bool get(size_t idx) const {
		assert(idx < max_num_bits);
		size_t arr_off = idx / NUM_BITS_LONG;
		size_t inside_off = idx % NUM_BITS_LONG;
		return ptr[arr_off].load(std::memory_order_relaxed) & (1L << inside_off);
	}

	void clear(size_t idx) {
		assert(idx < max_num_bits);
		size_t arr_off = idx / NUM_BITS_LONG;
		size_t inside_off = idx % NUM_BITS_LONG;
		ptr[arr_off].fetch_and(~(1L << inside_off), std::memory_order_relaxed);
	}

	void clear() {
		size_t num_longs = ROUNDUP(max_num_bits, NUM_BITS_LONG) / NUM_BITS_LONG;
		for (size_t i = 0; i < num_longs; i++)
			new (ptr + i) std::atomic_ulong();
	}
};

#endif
