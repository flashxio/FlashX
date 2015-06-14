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
