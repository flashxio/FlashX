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

#include <boost/format.hpp>

#include "log.h"

#include "vector_vector.h"
#include "mem_vv_store.h"
#include "local_vv_store.h"
#include "local_vec_store.h"

namespace fm
{

namespace detail
{

mem_vv_store::ptr mem_vv_store::cast(vec_store::ptr store)
{
	if (!store->is_in_mem() || store->get_entry_size() != 0) {
		BOOST_LOG_TRIVIAL(error) << "This isn't mem vv store\n";
		return mem_vv_store::ptr();
	}
	return std::static_pointer_cast<mem_vv_store>(store);
}

std::vector<off_t> mem_vv_store::get_rel_offs(off_t start, size_t len) const
{
	// The last entry shows the end of the last vector.
	std::vector<off_t> offs(len + 1);
	off_t start_off = get_vec_off(start);
	for (size_t i = 0; i < offs.size(); i++)
		offs[i] = get_vec_off(i + start) - start_off;
	return offs;
}

local_vec_store::const_ptr mem_vv_store::get_portion(off_t start,
		size_t len) const
{
	if (start + len > get_num_vecs()) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"can't get the portion [%1%, %2%)") % start % (start + len);
		return local_vec_store::const_ptr();
	}

	std::vector<off_t> offs = get_rel_offs(start, len);
	size_t num_eles = offs.back() / get_type().get_size();
	off_t start_ele = get_vec_off(start) / get_type().get_size();
	local_vec_store::ptr data(new local_cref_vec_store(get_mem_data().get(start_ele),
				start_ele, num_eles, get_type(), -1));

	return local_vv_store::ptr(new local_vv_store(start, offs, data));
}

local_vec_store::ptr mem_vv_store::get_portion(off_t start, size_t len)
{
	if (start + len > get_num_vecs()) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"can't get the portion [%1%, %2%)") % start % (start + len);
		return local_vec_store::ptr();
	}

	std::vector<off_t> offs = get_rel_offs(start, len);
	size_t num_eles = offs.back() / get_type().get_size();
	off_t start_ele = get_vec_off(start) / get_type().get_size();
	local_vec_store::ptr data(new local_ref_vec_store(get_mem_data().get(start_ele),
				start_ele, num_eles, get_type(), -1));

	return local_vv_store::ptr(new local_vv_store(start, offs, data));
}

}

}
