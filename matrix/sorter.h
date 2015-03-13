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
	virtual bool is_sorted(char *data, size_t num, bool decreasing) const = 0;
	virtual void sort_with_index(char *data, off_t *offs, size_t num,
			bool decreasing) const = 0;
	virtual void sort(char *data, size_t num, bool decreasing) const = 0;
	virtual void serial_sort(char *data, size_t num, bool decreasing) const = 0;
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
	virtual bool is_sorted(char *data, size_t num, bool decreasing) const;
	virtual void sort_with_index(char *data, off_t *offs, size_t num,
			bool decreasing) const;
	virtual void sort(char *data, size_t num, bool decreasing) const;
	virtual void serial_sort(char *data, size_t num, bool decreasing) const;
};

template<class T>
bool type_sorter<T>::is_sorted(char *data1, size_t num, bool decreasing) const
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

}

#endif
