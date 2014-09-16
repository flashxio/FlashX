#ifndef __VERTEX_POINTER_H__
#define __VERTEX_POINTER_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashGraph.
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
#include <stdlib.h>

class compute_vertex;

class compute_vertex_pointer
{
	static const int PART_BIT_LOC = 48;
	// TODO This is Linux specific.
	static const uint64_t PART_MASK = 0xFFFF000000000000UL;
	static const uint64_t ADDR_MASK = 0x0000FFFFFFFFFFFFUL;
	uint64_t addr;
public:
	static compute_vertex **conv(compute_vertex_pointer *arr) {
		assert(sizeof(compute_vertex_pointer) == sizeof(arr[0]));
		return (compute_vertex **) arr;
	}

	compute_vertex_pointer() {
		addr = 0;
	}

	explicit compute_vertex_pointer(compute_vertex *v) {
		assert(v);
		addr = (uint64_t) v;
		assert(addr > 100000);
	}

	compute_vertex_pointer(compute_vertex *v, bool part) {
		assert(v);
		assert(((uint64_t) v) > 100000);
		addr = ((uint64_t) v) + (((uint64_t) part) << PART_BIT_LOC);
	}

	bool is_part() const {
		return addr & PART_MASK;
	}

	bool is_valid() const {
		return (addr & ADDR_MASK) != 0;
	}

	compute_vertex *get() const {
		assert((addr & ADDR_MASK) != 0);
		return (compute_vertex *) (uintptr_t) (addr & ADDR_MASK);
	}

	compute_vertex *operator->() const {
		return get();
	}

	compute_vertex &operator*() const {
		return *get();
	}
};

#endif
