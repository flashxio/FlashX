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

bool mem_vector_vector::append(const mem_vector &vec)
{
	vector::resize(get_num_vecs() + 1);
	size_t vec_num_bytes = vec.get_length() * vec.get_entry_size();
	if (get_num_bytes() + vec_num_bytes > capacity)
		expand(get_num_bytes() + vec_num_bytes);
	assert(get_num_bytes() + vec_num_bytes <= capacity);
	memcpy(get_end(), vec.get_raw_arr(), vec_num_bytes);
	off_t new_off = get_num_bytes() + vec_num_bytes;
	vec_offs.push_back(new_off);
	return true;
}

vector::ptr mem_vector_vector::cat() const
{
	mem_vector::ptr ret = get_type().create_mem_vec(data, get_num_bytes());
	return std::static_pointer_cast<vector>(ret);
}

}
