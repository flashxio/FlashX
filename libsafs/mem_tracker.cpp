/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
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

#include <stdio.h>

#include "mem_tracker.h"
#include "concurrency.h"

static atomic_number<size_t> alloc_objs;
static atomic_number<size_t> alloc_bytes;
static bool mem_trace;

namespace {
class global_max
{
	volatile size_t value;
public:
	global_max() {
		value = 0;
	}

	bool update(size_t new_v) {
		if (new_v <= value)
			return false;

		bool ret = false;
		// There is a race condition, but we don't care.
		// We don't need to get the exact max value, but roughly the max
		// value. We can't use lock here, because the init function of
		// a program may allocate memory before invoking the constructor
		// of this class.
		if (new_v > value) {
			value = new_v;
			ret = true;
		}
		return ret;
	}

	size_t get() const {
		return value;
	}
};
}

static global_max max_objs;
static global_max max_bytes;
static global_max max_alloc;

#ifdef ENABLE_MEM_TRACE

void *operator new(size_t n) throw (std::bad_alloc)
{
	if (mem_trace) {
		size_t ret = alloc_objs.inc(1);
		max_objs.update(ret);
		ret = alloc_bytes.inc(n);
		max_bytes.update(ret);
		max_alloc.update(n);
		if (n >= 10 * 1024 * 1024)
			printf("allocate %ld\n", n);
	}
	void *p = malloc(n + sizeof(size_t));
	size_t *size = (size_t *) p;
	*size = n;

	return (void *) (size + 1);
}

void operator delete(void *p) throw ()
{
	if (p == NULL)
		return;

	size_t *size = (size_t *) (((char *) p) - sizeof(size_t));
	if (mem_trace) {
		alloc_objs.dec(1);
		alloc_bytes.dec(*size);
	}
	free(size);
}

void *operator new[](size_t n) throw (std::bad_alloc)
{
	if (mem_trace) {
		size_t ret = alloc_objs.inc(1);
		max_objs.update(ret);
		ret = alloc_bytes.inc(n);
		max_bytes.update(ret);
		max_alloc.update(n);
	}
	void *p = malloc(n + sizeof(size_t));
	size_t *size = (size_t *) p;
	*size = n;

	return (void *) (size + 1);
}

void operator delete[](void *p) throw ()
{
	if (p == NULL)
		return;

	size_t *size = (size_t *) (((char *) p) - sizeof(size_t));
	if (mem_trace) {
		alloc_objs.dec(1);
		alloc_bytes.dec(*size);
	}
	free(size);
}

#endif

size_t get_alloc_objs()
{
	return alloc_objs.get();
}

size_t get_alloc_bytes()
{
	return alloc_bytes.get();
}

size_t get_max_alloc_objs()
{
	return max_objs.get();
}

size_t get_max_alloc_bytes()
{
	return max_bytes.get();
}

size_t get_max_alloc()
{
	return max_alloc.get();
}

void mem_trace_start()
{
	mem_trace = true;
}

void mem_trace_stop()
{
	mem_trace = false;
}
