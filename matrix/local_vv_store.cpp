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

std::vector<off_t> local_vv_store::get_rel_offs(off_t start, size_t len) const
{
	// The last entry shows the end of the last vector.
	std::vector<off_t> new_offs(len + 1);
	off_t start_off = offs[start];
	for (size_t i = 0; i < new_offs.size(); i++)
		new_offs[i] = offs[i + start] - start_off;
	return new_offs;
}

local_vec_store::ptr local_vv_store::get_portion(off_t loc, size_t size)
{
	assert(get_raw_arr());
	assert(loc + size <= get_length());
	off_t start = get_global_start() + loc;
	std::vector<off_t> new_offs = get_rel_offs(loc, size);
	size_t entry_size = get_type().get_size();
	off_t rel_vec_start = offs[loc] / entry_size;
	off_t rel_vec_len = (new_offs.back() - new_offs.front()) / entry_size;
	return local_vec_store::ptr(new local_vv_store(start, new_offs,
				vec->get_portion(rel_vec_start, rel_vec_len)));
}

local_vec_store::const_ptr local_vv_store::get_portion(off_t loc,
		size_t size) const
{
	assert(get_raw_arr());
	assert(loc + size <= get_length());
	off_t start = get_global_start() + loc;
	std::vector<off_t> new_offs = get_rel_offs(loc, size);
	size_t entry_size = get_type().get_size();
	off_t rel_vec_start = offs[loc] / entry_size;
	off_t rel_vec_len = (new_offs.back() - new_offs.front()) / entry_size;
	return local_vec_store::ptr(new local_vv_store(start, new_offs,
				vec->get_portion(rel_vec_start, rel_vec_len)));
}

namespace detail
{

detail::mem_vv_store::ptr apply(const local_vv_store &store,
		const arr_apply_operate &op)
{
	const scalar_type &output_type = op.get_output_type();
	size_t out_size;
	// If the user can predict the number of output elements, we can create
	// a buffer of the expected size.
	if (op.get_num_out_eles() > 0)
		out_size = op.get_num_out_eles();
	else
		// If the user can't, we create a small buffer.
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
