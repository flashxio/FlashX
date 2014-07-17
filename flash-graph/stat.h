#ifndef __MY_STAT_H__
#define __MY_STAT_H__

/**
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

#include <limits.h>
#include <math.h>
#include <assert.h>

#include <map>
#include <vector>
#include <boost/foreach.hpp>

#include "FG_basic_types.h"

class log_hist_bucket
{
	size_t lower_bound;
	size_t upper_bound;
	size_t count;
public:
	log_hist_bucket() {
		count = 0;
		lower_bound = 0;
		upper_bound = INT_MAX;
	}

	log_hist_bucket(int idx) {
		count = 0;
		lower_bound = pow(10, idx);
		if (idx == 0)
			lower_bound = 0;
		upper_bound = pow(10, idx + 1);
	}

	size_t get_lower_bound() const {
		return lower_bound;
	}

	size_t get_upper_bound() const {
		return upper_bound;
	}

	size_t get_count() const {
		return count;
	}

	void inc_count(size_t num) {
		count += num;
	}
};

class log_histogram
{
	std::vector<log_hist_bucket> buckets;
public:
	log_histogram(int num_buckets): buckets(num_buckets) {
		for (int i = 0; i < num_buckets; i++)
			buckets[i] = log_hist_bucket(i);
	}

	log_hist_bucket &find_bucket(size_t v) {
		for (size_t i = 0; i < buckets.size(); i++)
			if (buckets[i].get_lower_bound() <= v
					&& v < buckets[i].get_upper_bound())
				return buckets[i];
		printf("can't find a bucket for %ld. All buckets cover [%ld, %ld).\n",
				v, buckets.front().get_lower_bound(),
				buckets.back().get_upper_bound());
		assert(0);
	}

	void add_value(size_t v) {
		find_bucket(v).inc_count(1);
	}

	int get_num_buckets() const {
		return buckets.size();
	}

	const log_hist_bucket &get_bucket(int idx) const {
		return buckets[idx];
	}

	void print(FILE *f) const {
		int non_empty = get_num_buckets() - 1;
		for (; non_empty > 0 && get_bucket(non_empty).get_count() == 0; non_empty--);
		for (int i = 0; i <= non_empty; i++) {
			fprintf(f, "[%ld, %ld): %ld\n", get_bucket(i).get_lower_bound(),
					get_bucket(i).get_upper_bound(),
					get_bucket(i).get_count());
		}
	}
};

/**
 * This count the apperances of each value of type T.
 */
template<class T>
class count_map
{
	typedef std::map<T, size_t> map_t;
	map_t map;
public:
	void add(const T &value) {
		typename map_t::iterator it = map.find(value);
		if (it == map.end())
			map.insert(std::pair<T, vsize_t>(value, 1));
		else
			it->second++;
	}

	void merge(const count_map &map) {
		BOOST_FOREACH(typename map_t::value_type v, map.map) {
			typename map_t::iterator it = this->map.find(v.first);
			if (it == this->map.begin())
				this->map.insert(std::pair<T, size_t>(v.first, v.second));
			else
				it->second += v.second;
		}
	}

	std::pair<T, size_t> get_max_count() const {
		size_t max_counts = 0;
		vertex_id_t max_v = T();
		BOOST_FOREACH(typename map_t::value_type v, map) {
			if (max_counts < v.second) {
				max_counts = v.second;
				max_v = v.first;
			}
		}
		return std::pair<T, size_t>(max_v, max_counts);
	}

	template<class ApplyFunc>
	void apply(ApplyFunc func) const {
		BOOST_FOREACH(const typename map_t::value_type v, map) {
			func(v);
		}
	}

	size_t get_size() const {
		return map.size();
	}

	bool exists(const T &key) const {
		typename map_t::const_iterator it = map.find(key);
		return it != map.end();
	}

	size_t get(const T &key) const {
		typename map_t::const_iterator it = map.find(key);
		assert(it != map.end());
		return it->second;
	}

	void remove(const T &key) {
		map.erase(key);
	}
};

#endif
