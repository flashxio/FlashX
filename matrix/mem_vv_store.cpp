/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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

#include "log.h"

#include "vector_vector.h"
#include "mem_vv_store.h"

namespace fm
{

namespace detail
{

mem_vv_store::mem_vv_store(const scalar_type &type): mem_vec_store(0, 0, type)
{
	store = smp_vec_store::create(0, type);
	// The offset of the first vector in the vector vector is 0.
	vec_offs.push_back(0);
}

mem_vv_store::mem_vv_store(const detail::raw_data_array &data,
		const std::vector<off_t> &offs, const scalar_type &type): mem_vec_store(
			vec_offs.size() - 1, 0, type)
{
	// There must be at least two locations to identify the vectors
	// in the data array.
	assert(offs.size() > 1);
	assert((size_t) offs.back() <= data.get_num_bytes());
	this->vec_offs = offs;
	store = smp_vec_store::create(data, type);
}

mem_vv_store::ptr mem_vv_store::cast(vec_store::ptr store)
{
	if (!store->is_in_mem() || store->get_entry_size() != 0) {
		BOOST_LOG_TRIVIAL(error) << "This isn't mem vv store\n";
		return mem_vv_store::ptr();
	}
	return std::static_pointer_cast<mem_vv_store>(store);
}

bool mem_vv_store::append(
		std::vector<vec_store::const_ptr>::const_iterator vec_it,
		std::vector<vec_store::const_ptr>::const_iterator vec_end)
{
	if (vec_it == vec_end)
		return true;

	bool is_vv = is_vector_vector(**vec_it);
	// This contains the real vector store (not vv store).
	std::vector<vec_store::const_ptr> vecs(vec_end - vec_it);
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
		if (is_vector_vector(**it)) {
			const mem_vv_store &vv = static_cast<const mem_vv_store &>(**it);
			vecs[it - vec_it] = vv.store;
			off_t off_end = vec_offs.back();
			size_t num_vecs = vv.get_num_vecs();
			for (size_t i = 0; i < num_vecs; i++) {
				off_end += vv.get_num_bytes(i);
				vec_offs.push_back(off_end);
			}
		}
		else {
			off_t off_end = vec_offs.back()
				+ (*it)->get_length() * (*it)->get_type().get_size();
			vec_offs.push_back(off_end);
			vecs[it - vec_it] = *it;
		}
	}

	// We have to append real vector store.
	store->append(vecs.begin(), vecs.end());
	vec_store::resize(vec_offs.size() - 1);

	return true;
}

bool mem_vv_store::append(const vec_store &vec)
{
	bool ret;

	if (is_vector_vector(vec)) {
		const mem_vv_store &vv = static_cast<const mem_vv_store &>(vec);
		// Construct the offset metadata of this vector vector.
		size_t num_vecs = vv.get_num_vecs();
		off_t off_end = get_num_bytes();
		for (size_t i = 0; i < num_vecs; i++) {
			off_end += vv.get_num_bytes(i);
			vec_offs.push_back(off_end);
		}

		ret = store->append(*vv.store);
		assert(ret);
		vec_store::resize(vec_offs.size() - 1);
	}
	else {
		size_t vec_num_bytes = vec.get_length() * vec.get_entry_size();
		vec_offs.push_back(get_num_bytes() + vec_num_bytes);

		ret = store->append(vec);
		assert(ret);
		vec_store::resize(vec_offs.size() - 1);
	}

	return ret;
}

mem_vv_store::const_ptr mem_vv_store::get_sub_vec_vec(off_t start,
		size_t len) const
{
	std::vector<off_t> offs(vec_offs.begin() + start,
			// The last entry shows the end of the last vector.
			vec_offs.begin() + start + len + 1);
	mem_vv_store::ptr ret = mem_vv_store::ptr(new mem_vv_store(get_type()));
	ret->vec_offs = offs;
	ret->store = store;
	ret->vec_store::resize(ret->vec_offs.size() - 1);
	return ret;
}

vec_store::const_ptr mem_vv_store::cat() const
{
	smp_vec_store::ptr ret = smp_vec_store::cast(store->shallow_copy());
	off_t start = vec_offs.front() / get_type().get_size();
	size_t len = (vec_offs.back() - vec_offs.front()) / get_type().get_size();
	ret->expose_sub_vec(start, len);
	return ret;
}

}

}
