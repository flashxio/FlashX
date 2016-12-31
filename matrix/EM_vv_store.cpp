/*
 * Copyright 2016 Open Connectome Project (http://openconnecto.me)
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

#include <sys/time.h>

#include <boost/format.hpp>

#include "EM_vv_store.h"
#include "local_vv_store.h"

namespace fm
{

namespace detail
{

local_vec_store::ptr EM_vv_store::get_portion_async(off_t start,
		size_t len, portion_compute::ptr compute) const
{
	if (start + len > get_num_vecs()) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"can't get the portion [%1%, %2%)") % start % (start + len);
		return local_vec_store::ptr();
	}

	off_t start_ele = get_vec_off(start) / get_type().get_size();
	local_vec_store::ptr data = get_EM_data().get_portion_async(start_ele,
			get_num_eles(start, len), compute);
	return local_vv_store::ptr(new local_vv_store(start, get_off_it(start),
				get_off_it(start + len + 1), data));
}

}

}
