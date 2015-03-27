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

#include "NUMA_mapper.h"

namespace fm
{

namespace detail
{

std::vector<size_t> NUMA_mapper::cal_local_lengths(size_t len) const
{
	std::vector<size_t> ret(get_num_nodes());
	size_t last_range_id = len >> range_size_log;
	size_t last_range_size = len & range_mask;
	auto phy_loc = map2physical(last_range_id << range_size_log);
	ret[phy_loc.first] = phy_loc.second + last_range_size;
	size_t range_off = last_range_id << range_size_log;
	for (size_t i = 0; i < get_num_nodes() - 1 && range_off >= range_size; i++) {
		range_off -= range_size;
		phy_loc = map2physical(range_off);
		assert(ret[phy_loc.first] == 0);
		ret[phy_loc.first] = phy_loc.second + range_size;
	}
	return ret;
}

}

}
