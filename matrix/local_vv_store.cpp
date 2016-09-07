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

#include "local_vv_store.h"
#include "local_vec_store.h"
#include "mem_vv_store.h"

namespace fm
{

local_vv_store::local_vv_store(off_t global_start,
		std::vector<off_t>::const_iterator start,
		std::vector<off_t>::const_iterator end,
		local_vec_store::ptr vec): local_vec_store(
			static_cast<const local_vec_store &>(*vec).get_raw_arr(),
			vec->get_raw_arr(), global_start, end - start - 1, 0,
			vec->get_type(), -1)
{
	assert(!detail::vv_store::is_vector_vector(*vec));
	this->off_start = start;
	this->off_end = end;
	this->vec = vec;
}

local_vec_store::ptr local_vv_store::get_portion(off_t loc, size_t size)
{
	assert(get_raw_arr());
	assert(loc + size <= get_length());
	off_t global_start = get_global_start() + loc;
	loc += get_local_start();
	off_t rel_vec_start = get_vec_off(loc) / get_type().get_size();
	return local_vec_store::ptr(new local_vv_store(global_start, get_off_it(loc),
				get_off_it(loc + size + 1),
				vec->get_portion(rel_vec_start, get_num_eles(loc, size))));
}

local_vec_store::const_ptr local_vv_store::get_portion(off_t loc,
		size_t size) const
{
	assert(get_raw_arr());
	assert(loc + size <= get_length());
	off_t global_start = get_global_start() + loc;
	loc += get_local_start();
	off_t rel_vec_start = get_vec_off(loc) / get_type().get_size();
	return local_vec_store::ptr(new local_vv_store(global_start, get_off_it(loc),
				get_off_it(loc + size + 1),
				vec->get_portion(rel_vec_start, get_num_eles(loc, size))));
}

namespace detail
{

detail::mem_vv_store::ptr apply(const local_vv_store &store,
		const arr_apply_operate &op)
{
	const scalar_type &output_type = op.get_output_type();
	size_t out_size = op.get_num_out_eles(store.get_length(0));
	// If the user can't predict the number of output elements, we create
	// a small buffer.
	if (out_size <= 0)
		out_size = 16;
	local_buf_vec_store buf(0, out_size, output_type, -1);

	detail::mem_vv_store::ptr ret = detail::mem_vv_store::create(output_type);
	for (size_t i = 0; i < store.get_num_vecs(); i++) {
		local_cref_vec_store lvec(store.get_raw_arr(i), 0,
				store.get_length(i), store.get_type(), -1);
		op.run(lvec, buf);
		ret->append(buf);
	}
	return ret;
}

}

}
