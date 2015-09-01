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

#include "vec_store.h"
#include "mem_vec_store.h"
#include "EM_vector.h"

namespace fm
{

namespace detail
{

vec_store::ptr vec_store::create(size_t length, const scalar_type &type,
		bool in_mem)
{
	if (in_mem)
		return smp_vec_store::create(length, type);
	else
		return EM_vec_store::create(length, type);
}

size_t vec_store::copy_to(char *data, size_t num_eles) const
{
	size_t num_copy_eles = std::min(get_length(), num_eles);
	size_t portion_size = get_portion_size();
	size_t entry_size = get_type().get_size();
	for (size_t idx = 0; idx < num_copy_eles; idx += portion_size) {
		size_t len = std::min(portion_size, num_copy_eles - idx);
		local_vec_store::const_ptr store = get_portion(idx, len);
		assert(store);
		memcpy(data + idx * entry_size, store->get_raw_arr(), len * entry_size);
	}
	return num_copy_eles;
}

}
}
