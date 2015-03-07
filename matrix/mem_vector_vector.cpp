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

#include <string.h>
#include <assert.h>

#include "mem_vector_vector.h"
#include "mem_vector.h"
#include "data_frame.h"
#include "factor.h"

namespace fm
{

void mem_vector_vector::expand(size_t min)
{
	for (; capacity < min; capacity *= 2);
	std::shared_ptr<char> new_data = std::shared_ptr<char>(
			(char *) malloc(capacity), deleter());
	memcpy(new_data.get(), data.get(), get_num_bytes());
	data = new_data;
}

bool mem_vector_vector::append(std::vector<vector::ptr>::const_iterator vec_it,
			std::vector<vector::ptr>::const_iterator vec_end)
{
	size_t tot_bytes = 0;
	for (auto it = vec_it; it != vec_end; it++) {
		if (!(*it)->is_in_mem()) {
			BOOST_LOG_TRIVIAL(error)
				<< "Not support appending an ext-mem vector";
			return false;
		}
		if (get_type() != (*it)->get_type()) {
			BOOST_LOG_TRIVIAL(error) << "The two vectors don't have the same type";
			return false;
		}
		tot_bytes += (*it)->get_length() * (*it)->get_type().get_size();
	}

	vector::resize(get_num_vecs() + (vec_end - vec_it));
	if (get_num_bytes() + tot_bytes > capacity)
		expand(get_num_bytes() + tot_bytes);
	for (auto it = vec_it; it != vec_end; it++) {
		const mem_vector &mem_vec = (const mem_vector &) **it;
		size_t vec_num_bytes = mem_vec.get_length() * mem_vec.get_type().get_size();
		assert(get_num_bytes() + vec_num_bytes <= capacity);
		memcpy(get_end(), mem_vec.get_raw_arr(), vec_num_bytes);
		off_t new_off = get_num_bytes() + vec_num_bytes;
		vec_offs.push_back(new_off);
	}
	return true;
}

bool mem_vector_vector::append(const vector &vec)
{
	if (!vec.is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "Not support appending an ext-mem vector";
		return false;
	}
	if (get_type() != vec.get_type()) {
		BOOST_LOG_TRIVIAL(error) << "The two vectors don't have the same type";
		return false;
	}

	vector::resize(get_num_vecs() + 1);
	const mem_vector &mem_vec = (const mem_vector &) vec;
	size_t vec_num_bytes = mem_vec.get_length() * mem_vec.get_entry_size();
	if (get_num_bytes() + vec_num_bytes > capacity)
		expand(get_num_bytes() + vec_num_bytes);
	assert(get_num_bytes() + vec_num_bytes <= capacity);
	memcpy(get_end(), mem_vec.get_raw_arr(), vec_num_bytes);
	off_t new_off = get_num_bytes() + vec_num_bytes;
	vec_offs.push_back(new_off);
	return true;
}

vector::ptr mem_vector_vector::cat() const
{
	mem_vector::ptr ret = get_type().create_mem_vec(data, get_num_bytes());
	return std::static_pointer_cast<vector>(ret);
}

class sorted_gr_label_operate: public gr_apply_operate<mem_vector>
{
	const factor_vector &labels;
	const gr_apply_operate<sub_vector_vector> &op;
	const mem_vector_vector &vv;
public:
	sorted_gr_label_operate(const factor_vector &_labels,
			const gr_apply_operate<sub_vector_vector> &_op,
			const mem_vector_vector &_vv): labels(_labels), op(_op), vv(_vv) {
	}

	virtual void run(const void *key, const mem_vector &val, mem_vector &vec) const;

	virtual const scalar_type &get_key_type() const {
		return labels.get_type();
	}

	virtual const scalar_type &get_output_type() const {
		return op.get_output_type();
	}

	virtual size_t get_num_out_eles() const {
		return 0;
	}
};

void sorted_gr_label_operate::run(const void *key, const mem_vector &val,
		mem_vector &vec) const
{
	off_t off_bytes = ((const char *) key) - labels.get_raw_arr();
	off_t off = off_bytes / labels.get_type().get_size();
	size_t len = val.get_length();
	std::vector<off_t> vec_idxs(len);
	for (size_t i = 0; i < len; i++)
		vec_idxs[i] = i + off;
	op.run(key, sub_vector_vector(vv, vec_idxs), vec);
}

vector_vector::ptr mem_vector_vector::groupby(const factor_vector &labels,
			const gr_apply_operate<sub_vector_vector> &op) const
{
	if (labels.get_length() != get_num_vecs()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The label vector has a different length from the vector vector";
		return vector_vector::ptr();
	}
	assert(labels.is_sorted());

	data_frame::ptr df = labels.groupby(
			sorted_gr_label_operate(labels, op, *this), false);
	return vector_vector::cast(df->get_vec("agg"));
}

}
