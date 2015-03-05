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

#ifndef __DATA_FRAME_H__
#define __DATA_FRAME_H__

#include <assert.h>

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

#include "vector.h"

namespace fm
{

typedef std::pair<std::string, vector::ptr> named_vec_t;

template<class T> class gr_apply_operate;
class vector_vector;

/**
 * This implements the data frame in R.
 * It also resembles a SQL table. Data is stored column-wise and each column
 * has its own data type
 */
class data_frame
{
	std::vector<named_vec_t> named_vecs;
	std::unordered_map<std::string, vector::ptr> vec_map;

protected:
	bool expose_portion(off_t loc, size_t length);

	void set_vec(size_t off, vector::ptr vec) {
		std::string name = named_vecs[off].first;
		named_vecs[off].second = vec;
		vec_map[name] = vec;
	}

	data_frame() {
	}

	data_frame(const std::vector<named_vec_t> &named_vecs) {
		this->named_vecs = named_vecs;
		for (auto it = named_vecs.begin(); it != named_vecs.end(); it++)
			vec_map.insert(*it);
	}
public:
	typedef std::shared_ptr<data_frame> ptr;

	bool add_vec(const std::string &name, vector::ptr vec) {
		if (get_num_vecs() > 0 && vec->get_length() != get_num_entries())
			return false;
		named_vecs.push_back(named_vec_t(name, vec));
		vec_map.insert(named_vec_t(name, vec));
		return true;
	}

	/**
	 * This method appends multiple data frames to this data frame.
	 * All data frames should have the same number of columns and the columns
	 * should have the same names.
	 */
	bool append(std::vector<data_frame::ptr>::const_iterator begin,
			std::vector<data_frame::ptr>::const_iterator end);
	bool append(data_frame::ptr df);

	vector::ptr get_vec(size_t off) {
		return named_vecs[off].second;
	}

	const std::string &get_vec_name(size_t off) const {
		return named_vecs[off].first;
	}

	vector::ptr get_vec(const std::string &name) {
		auto it = vec_map.find(name);
		if (it == vec_map.end())
			return vector::ptr();
		else
			return it->second;
	}

	const vector &get_vec_ref(size_t off) const {
		return *named_vecs[off].second;
	}

	const vector &get_vec_ref(const std::string &name) const {
		auto it = vec_map.find(name);
		assert(it != vec_map.end());
		return *it->second;
	}

	size_t get_num_vecs() const {
		return named_vecs.size();
	}

	size_t get_num_entries() const {
		return named_vecs[0].second->get_length();
	}

	/**
	 * We group the rows of the data frame by the values in the specified column.
	 */
	virtual std::shared_ptr<vector_vector> groupby(const std::string &col_name,
			gr_apply_operate<data_frame> &op) const = 0;
	virtual bool sort(const std::string &col_name) = 0;
	virtual bool is_sorted(const std::string &col_name) const = 0;
};

}

#endif
