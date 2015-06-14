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

#include "vv_store.h"

namespace fm
{

namespace detail
{

vv_store::vv_store(const scalar_type &type, bool in_mem): vec_store(
		0, 0, type, in_mem)
{
	store = vec_store::create(0, type, in_mem);
	// The offset of the first vector in the vector vector is 0.
	vec_offs.push_back(0);
}

vv_store::vv_store(const std::vector<off_t> &offs,
		vec_store::ptr store): vec_store(offs.size() - 1, 0,
			store->get_type(), store->is_in_mem())
{
	this->vec_offs = offs;
	this->store = store;
}

vec_store::const_ptr vv_store::cat() const
{
	return store->shallow_copy();
}

vec_store::ptr vv_store::deep_copy() const
{
	vv_store::ptr ret
		= vv_store::cast(const_cast<vv_store *>(this)->shallow_copy());
	ret->vec_offs = this->vec_offs;
	ret->store = this->store->deep_copy();
	return ret;
}

bool vv_store::append(
		std::vector<vec_store::const_ptr>::const_iterator vec_it,
		std::vector<vec_store::const_ptr>::const_iterator vec_end)
{
	if (vec_it == vec_end)
		return true;

	size_t orig_num_offs = vec_offs.size();
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
			const vv_store &vv = static_cast<const vv_store &>(**it);
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
	bool ret = store->append(vecs.begin(), vecs.end());
	if (ret)
		vec_store::resize(vec_offs.size() - 1);
	else
		// recover from the failure and go back to the original state.
		vec_offs.resize(orig_num_offs);

	return ret;
}

bool vv_store::append(const vec_store &vec)
{
	if (!vec.is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "Not support appending an ext-mem vector";
		return false;
	}

	bool ret;
	size_t orig_num_offs = vec_offs.size();
	if (is_vector_vector(vec)) {
		const vv_store &vv = static_cast<const vv_store &>(vec);
		// Construct the offset metadata of this vector vector.
		size_t num_vecs = vv.get_num_vecs();
		off_t off_end = get_num_bytes();
		for (size_t i = 0; i < num_vecs; i++) {
			off_end += vv.get_num_bytes(i);
			vec_offs.push_back(off_end);
		}

		ret = store->append(*vv.store);
	}
	else {
		size_t vec_num_bytes = vec.get_length() * vec.get_entry_size();
		vec_offs.push_back(get_num_bytes() + vec_num_bytes);

		ret = store->append(vec);
	}
	if (ret)
		vec_store::resize(vec_offs.size() - 1);
	else
		// recover from the failure and go back to the original state.
		vec_offs.resize(orig_num_offs);

	return ret;
}

}

}
