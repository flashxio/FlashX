#ifndef __MY_SORTER_H__
#define __MY_SORTER_H__

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

#include <assert.h>
#include <memory>
#if defined(_OPENMP)
#include <parallel/algorithm>
#else
#include <algorithm>
#endif

namespace fm
{

class sorter
{
public:
	virtual bool is_sorted(const char *data, size_t num, bool decreasing) const = 0;
	virtual void sort_with_index(char *data, off_t *offs, size_t num,
			bool decreasing) const = 0;
	virtual void sort(char *data, size_t num, bool decreasing) const = 0;
	virtual void serial_sort(char *data, size_t num, bool decreasing) const = 0;
	virtual void merge(
			const std::vector<std::pair<const char *, const char *> > &arrs,
			char *output, size_t out_num) const = 0;
	virtual void merge(
			const std::vector<std::pair<const char *, const char *> > &arrs,
			const std::vector<std::pair<int, off_t> > &merge_index,
			char *output, size_t out_num) const = 0;
	virtual void merge_with_index(
			const std::vector<std::pair<const char *, const char *> > &arrs,
			char *output, size_t out_num,
			std::vector<std::pair<int, off_t> > &merge_index) const = 0;
};

template<class T>
class type_sorter: public sorter
{
	struct {
		bool operator()(const T &e1, const T &e2) const {
			return e1 < e2;
		}
	} entry_less;
	struct {
		bool operator()(const T &e1, const T &e2) const {
			return e1 > e2;
		}
	} entry_greater;

public:
	virtual bool is_sorted(const char *data, size_t num, bool decreasing) const;
	virtual void sort_with_index(char *data, off_t *offs, size_t num,
			bool decreasing) const;
	virtual void sort(char *data, size_t num, bool decreasing) const;
	virtual void serial_sort(char *data, size_t num, bool decreasing) const;
	virtual void merge(
			const std::vector<std::pair<const char *, const char *> > &arrs,
			char *output, size_t out_num) const;
	virtual void merge(
			const std::vector<std::pair<const char *, const char *> > &arrs,
			const std::vector<std::pair<int, off_t> > &merge_index,
			char *output, size_t out_num) const;
	virtual void merge_with_index(
			const std::vector<std::pair<const char *, const char *> > &arrs,
			char *output, size_t out_num,
			std::vector<std::pair<int, off_t> > &merge_index) const;
};

template<class T>
bool type_sorter<T>::is_sorted(const char *data1, size_t num, bool decreasing) const
{
	T *data = (T *) data1;
	if (decreasing)
		return std::is_sorted(data, data + num, entry_greater);
	else
		return std::is_sorted(data, data + num, entry_less);
}

template<class T>
void type_sorter<T>::sort_with_index(char *data1, off_t *offs, size_t num,
		bool decreasing) const
{
	T *data = (T *) data1;
	struct indexed_entry {
		T val;
		off_t idx;
	};

	std::unique_ptr<indexed_entry[]> entries
		= std::unique_ptr<indexed_entry[]>(new indexed_entry[num]);
#pragma omp parallel for
	for (size_t i = 0; i < num; i++) {
		entries[i].val = data[i];
		entries[i].idx = i;
	}

	struct {
		bool operator()(const indexed_entry &e1, const indexed_entry &e2) const {
			return e1.val < e2.val;
		}
	} entry_less;
	struct {
		bool operator()(const indexed_entry &e1, const indexed_entry &e2) const {
			return e1.val > e2.val;
		}
	} entry_greater;

	indexed_entry *start = entries.get();
	indexed_entry *end = start + num;
#if defined(_OPENMP)
	if (decreasing)
		__gnu_parallel::sort(start, end, entry_greater);
	else
		__gnu_parallel::sort(start, end, entry_less);
#else
	if (decreasing)
		std::sort(start, end, entry_greater);
	else
		std::sort(start, end, entry_less);
#endif
#pragma omp parallel for
	for (size_t i = 0; i < num; i++) {
		data[i] = start[i].val;
		offs[i] = start[i].idx;
	}
}

template<class T>
void type_sorter<T>::sort(char *data1, size_t num, bool decreasing) const
{
	T *data = (T *) data1;
	T *start = (T *) data;
	T *end = start + num;
#if defined(_OPENMP)
	if (decreasing)
		__gnu_parallel::sort(start, end, entry_greater);
	else
		__gnu_parallel::sort(start, end, entry_less);
#else
	if (decreasing)
		std::sort(start, end, entry_greater);
	else
		std::sort(start, end, entry_less);
#endif
}

template<class T>
void type_sorter<T>::serial_sort(char *data1, size_t num, bool decreasing) const
{
	T *data = (T *) data1;
	T *start = (T *) data;
	T *end = start + num;
	if (decreasing)
		std::sort(start, end, entry_greater);
	else
		std::sort(start, end, entry_less);
}

template<class T>
void type_sorter<T>::merge(
		const std::vector<std::pair<const char *, const char *> > &raw_arrs,
		char *output, size_t out_num) const
{
	std::vector<std::pair<T *, T *> > arrs(raw_arrs.size());
	for (size_t i = 0; i < arrs.size(); i++)
		arrs[i] = std::pair<T *, T *>((T *) raw_arrs[i].first,
				(T *) raw_arrs[i].second);
	__gnu_parallel::multiway_merge(arrs.begin(), arrs.end(), (T *) output,
			out_num, entry_less);
}

/*
 * Merge multiple arrays according to the specified locations.
 * Here I assume there are a few number of arrays to merge.
 */
template<class T>
void type_sorter<T>::merge(
		const std::vector<std::pair<const char *, const char *> > &arrs,
		const std::vector<std::pair<int, off_t> > &merge_index,
		char *output, size_t out_num) const
{
	T *t_output = (T *) output;
#pragma omp parallel for
	for (size_t i = 0; i < out_num; i++) {
		int arr_idx = merge_index[i].first;
		off_t off_in_arr = merge_index[i].second;
		const T *t_arr = (const T *) arrs[arr_idx].first;
		assert(&t_arr[off_in_arr] <= (const T *) arrs[arr_idx].second);
		t_output[i] = t_arr[off_in_arr];
	}
}

/*
 * Get the length of an array indicated by the pair (`first' is the beginning
 * of the array and `second' is the end of the array.
 */
template<class T>
size_t get_length(const std::pair<const char *, const char *> &arr)
{
	return (arr.second - arr.first) / sizeof(T);
}

/*
 * Merge multiple arrays and return the merged result as well as how
 * the arrays are merged.
 * Here I assume there are a few number of arrays to merge.
 */
template<class T>
void type_sorter<T>::merge_with_index(
		const std::vector<std::pair<const char *, const char *> > &arrs,
		char *output, size_t out_num,
		std::vector<std::pair<int, off_t> > &merge_index) const
{
	struct indexed_entry {
		T val;
		int arr_idx;
		off_t off_in_arr;
	};
	std::unique_ptr<indexed_entry[]> buf(new indexed_entry[out_num]);
	// Move data from `arrs' to `buf' in parallel.
#pragma omp parallel
	{
		size_t avg_part_len = ceil(((double) out_num) / omp_get_num_threads());
		size_t thread_id = omp_get_thread_num();
		size_t start = thread_id * avg_part_len;
		size_t part_len = std::min(out_num - start, avg_part_len);

		// Find the first array for the current thread.
		size_t curr_arr_idx = 0;
		size_t i = 0;
		while (true) {
			// If the array is empty, it works fine.
			size_t num_eles = get_length<T>(arrs[curr_arr_idx]);
			if (i + num_eles > start)
				break;
			i += num_eles;
			curr_arr_idx++;
		}
		assert(start >= i);
		off_t curr_off_in_arr = start - i;
		for (size_t i = 0; i < part_len; i++) {
			const T *curr_ptr
				= ((const T *) arrs[curr_arr_idx].first) + curr_off_in_arr;
			assert(curr_ptr < (const T *) arrs[curr_arr_idx].second);
			buf[start + i].val = *curr_ptr;
			buf[start + i].arr_idx = curr_arr_idx;
			buf[start + i].off_in_arr = curr_off_in_arr;
			// If the current pointer points to the last element in the array,
			// switch to the next array.
			if (curr_ptr == ((const T *) arrs[curr_arr_idx].second) - 1) {
				curr_arr_idx++;
				// We need to skip the empty arrays.
				while (curr_arr_idx < arrs.size()
						&& get_length<T>(arrs[curr_arr_idx]) == 0)
					curr_arr_idx++;
				curr_off_in_arr = 0;
				if (i + 1 < part_len)
					assert(curr_arr_idx < arrs.size());
			}
			else
				curr_off_in_arr++;
		}
	}

	std::vector<std::pair<indexed_entry *, indexed_entry *> > indexed_arrs(
			arrs.size());
	size_t off = 0;
	for (size_t i = 0; i < arrs.size(); i++) {
		size_t len = get_length<T>(arrs[i]);
		indexed_arrs[i] = std::pair<indexed_entry *, indexed_entry *>(
				&buf[off], &buf[off + len]);
		off += len;
	}
	assert(off == out_num);

	struct {
		bool operator()(const indexed_entry &e1, const indexed_entry &e2) const {
			return e1.val < e2.val;
		}
	} entry_less;
	std::unique_ptr<indexed_entry[]> merge_res(new indexed_entry[out_num]);
	__gnu_parallel::multiway_merge(indexed_arrs.begin(), indexed_arrs.end(),
			merge_res.get(), out_num, entry_less);
	T *t_output = (T *) output;
	for (size_t i = 0; i < out_num; i++) {
		t_output[i] = merge_res[i].val;
		merge_index[i].first = merge_res[i].arr_idx;
		merge_index[i].second = merge_res[i].off_in_arr;
	}
}

}

#endif
