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

#include "common.h"

#include "mem_vector_vector.h"
#include "mem_vector.h"
#include "data_frame.h"
#include "factor.h"

namespace fm
{

void mem_vector_vector::append_mem_vector(const mem_vector &mem_vec)
{
	vector::resize(get_num_vecs() + 1);
	size_t vec_num_bytes = mem_vec.get_length() * mem_vec.get_entry_size();
	if (get_num_bytes() + vec_num_bytes > data.get_num_bytes())
		data.expand(get_num_bytes() + vec_num_bytes);
	assert(get_num_bytes() + vec_num_bytes <= data.get_num_bytes());
	memcpy(get_end(), mem_vec.get_raw_arr(), vec_num_bytes);
	off_t new_off = get_num_bytes() + vec_num_bytes;
	vec_offs.push_back(new_off);
}

void mem_vector_vector::append_mem_vectors(
		std::vector<vector::ptr>::const_iterator vec_it,
		std::vector<vector::ptr>::const_iterator vec_end)
{
	size_t tot_bytes = 0;
	for (auto it = vec_it; it != vec_end; it++)
		tot_bytes += (*it)->get_length() * (*it)->get_type().get_size();
	vector::resize(get_num_vecs() + (vec_end - vec_it));
	if (get_num_bytes() + tot_bytes > data.get_num_bytes())
		data.expand(get_num_bytes() + tot_bytes);

	for (auto it = vec_it; it != vec_end; it++) {
		const mem_vector &mem_vec = (const mem_vector &) **it;
		size_t vec_num_bytes = mem_vec.get_length() * mem_vec.get_type().get_size();
		assert(get_num_bytes() + vec_num_bytes <= data.get_num_bytes());
		memcpy(get_end(), mem_vec.get_raw_arr(), vec_num_bytes);
		off_t new_off = get_num_bytes() + vec_num_bytes;
		vec_offs.push_back(new_off);
	}
}

void mem_vector_vector::append_mem_vv(const mem_vector_vector &vv)
{
	// Copy all data in `vv' to this vector vector.
	size_t vv_num_bytes = vv.get_num_bytes();
	if (get_num_bytes() + vv_num_bytes > data.get_num_bytes())
		data.expand(get_num_bytes() + vv_num_bytes);
	assert(get_num_bytes() + vv_num_bytes <= data.get_num_bytes());
	memcpy(get_end(), vv.get_raw_data(), vv_num_bytes);

	// Construct the offset metadata of this vector vector.
	size_t num_vecs = vv.get_num_vecs();
	vector::resize(get_num_vecs() + num_vecs);
	off_t off_end = get_num_bytes();
	for (size_t i = 0; i < num_vecs; i++) {
		off_end += vv.get_num_bytes(i);
		vec_offs.push_back(off_end);
	}
}

void mem_vector_vector::append_mem_vvs(
		std::vector<vector_vector::ptr>::const_iterator vec_it,
		std::vector<vector_vector::ptr>::const_iterator vec_end)
{
	size_t tot_bytes = 0;
	size_t tot_num_vecs = 0;
	for (auto it = vec_it; it != vec_end; it++) {
		mem_vector_vector &vv = (mem_vector_vector &) **it;
		tot_bytes += vv.get_num_bytes();
		tot_num_vecs += vv.get_num_vecs();
	}
	vector::resize(get_num_vecs() + tot_num_vecs);
	if (get_num_bytes() + tot_bytes > data.get_num_bytes())
		data.expand(get_num_bytes() + tot_bytes);

	for (auto it = vec_it; it != vec_end; it++) {
		const mem_vector_vector &vv = (const mem_vector_vector &) **it;
		size_t vv_num_bytes = vv.get_num_bytes();
		assert(get_num_bytes() + vv_num_bytes <= data.get_num_bytes());
		memcpy(get_end(), vv.get_raw_data(), vv_num_bytes);

		off_t off_end = get_num_bytes();
		size_t num_vecs = vv.get_num_vecs();
		for (size_t i = 0; i < num_vecs; i++) {
			off_end += vv.get_num_bytes(i);
			vec_offs.push_back(off_end);
		}
	}
}

bool mem_vector_vector::append(std::vector<vector::ptr>::const_iterator vec_it,
			std::vector<vector::ptr>::const_iterator vec_end)
{
	if (vec_it == vec_end)
		return true;

	bool is_vv = is_vector_vector(**vec_it);
	std::vector<vector_vector::ptr> vvs;
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
		if (is_vv != is_vector_vector(**it)) {
			BOOST_LOG_TRIVIAL(error) << "Not all vectors contain the same";
			return false;
		}
		if (is_vector_vector(**it))
			vvs.push_back(vector_vector::cast(*it));
	}
	if (is_vv)
		append_mem_vvs(vvs.begin(), vvs.end());
	else
		append_mem_vectors(vec_it, vec_end);

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

	if (is_vector_vector(vec))
		append_mem_vv((const mem_vector_vector &) vec);
	else
		append_mem_vector((const mem_vector &) vec);

	return true;
}

vector::ptr mem_vector_vector::cat() const
{
	mem_vector::ptr ret = mem_vector::create(data, get_num_bytes(), get_type());
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

vector_vector::ptr mem_vector_vector::apply(const arr_apply_operate &op) const
{
	int num_parts = std::min((size_t) get_num_omp_threads(), get_num_vecs());
	std::vector<vector_vector::ptr> vvs(num_parts);
	size_t part_num_vecs = ceil(((double) get_num_vecs()) / num_parts);
#pragma omp parallel for
	for (int i = 0; i < num_parts; i++) {
		// TODO Here I divide the vector_vector evenly in terms of the number
		// of vectors in each partition. It potentially has load balancing
		// problems.
		off_t start = part_num_vecs * i;
		off_t end = std::min(part_num_vecs * (i + 1), get_num_vecs());
		// When there aren't many vectors, this may happen.
		if (end <= start)
			continue;
		mem_vector_vector::ptr sub_vv = this->get_sub_vec_vec(start,
				end - start);
		vvs[i] = mem_vector_vector::cast(sub_vv)->serial_apply(op);
	}

	mem_vector_vector::ptr ret = mem_vector_vector::cast(vvs[0]);
	// It's possible that some threads didn't generate results.
	// It's guaranteed that the non-empty results are in the front.
	int num_non_empty = 0;
	for (int i = 0; i < num_parts; i++) {
		if (vvs[i]) {
			assert(i == num_non_empty);
			num_non_empty++;
		}
	}
	if (num_parts > 1)
		ret->append_mem_vvs(vvs.begin() + 1, vvs.begin() + num_non_empty);
	return ret;
}

vector_vector::ptr mem_vector_vector::serial_apply(
		const arr_apply_operate &op) const
{
	const scalar_type &output_type = op.get_output_type();
	mem_vector_vector::ptr ret = mem_vector_vector::create(output_type);
	mem_vector::ptr vec = mem_vector::cast(this->flatten());
	mem_vector::ptr buf = output_type.create_mem_vec(1);
	// It's possible that this is a sub-vector.
	off_t off = vec->get_sub_start();
	for (size_t i = 0; i < get_num_vecs(); i++) {
		// TODO the off is the absolute offset in the vector.
		vec->expose_sub_vec(off, this->get_length(i));
		off += this->get_length(i);
		op.run(*vec, *buf);
		ret->append(*buf);
	}
	return ret;
}

mem_vector_vector::ptr mem_vector_vector::get_sub_vec_vec(off_t start,
		size_t len) const
{
	std::vector<off_t> offs(vec_offs.begin() + start,
			// The last entry shows the end of the last vector.
			vec_offs.begin() + start + len + 1);
	return mem_vector_vector::ptr(new mem_vector_vector(data, offs,
				get_type()));
}

vector::ptr mem_vector_vector::flatten() const
{
	assert(vec_offs.size() > 1);
	size_t tot_len = vec_offs.back() / get_type().get_size();
	size_t off = vec_offs.front() / get_type().get_size();
	mem_vector::ptr ret = mem_vector::create(data, vec_offs.back(), get_type());
	// This might be a sub vector vector.
	if (vec_offs.front() > 0)
		ret->expose_sub_vec(off, tot_len - off);
	return ret;
}

}
