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

#include <numa.h>

#include "raw_data_array.h"

namespace fm
{

namespace detail
{

namespace
{

class deleter
{
	size_t size;
public:
	deleter(size_t size) {
		this->size = size;
	}

	void operator()(char *addr) {
		numa_free(addr, size);
	}
};

}

raw_data_array::raw_data_array(size_t num_bytes, int node_id)
{
	this->node_id = node_id;
	this->num_bytes = num_bytes;
	void *addr = numa_alloc_onnode(num_bytes, node_id);
	data = std::shared_ptr<char>((char *) addr, deleter(num_bytes));
	start = data.get();
	num_used_bytes = num_bytes;
}

}

}
