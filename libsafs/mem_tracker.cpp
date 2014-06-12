/**
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
#include "debugger.h"

atomic_number<size_t> alloc_objs;
atomic_number<size_t> alloc_bytes;

#ifdef ENABLE_MEM_TRACE

class print_mem_tracker_task: public debug_task
{
public:
	void run() {
		printf("alloc %ld objs and %ld bytes\n", alloc_objs.get(),
				alloc_bytes.get());
	}
};

void init_mem_tracker()
{
	printf("register mem_tracker to debug\n");
	debug.register_task(new print_mem_tracker_task());
}

void *operator new(size_t n) throw (std::bad_alloc)
{
	alloc_objs.inc(1);
	alloc_bytes.inc(n);
	void *p = malloc(n + sizeof(size_t));
	size_t *size = (size_t *) p;
	*size = n;

	return (void *) (size + 1);
}

void operator delete(void *p) throw ()
{
	if (p == NULL)
		return;

	alloc_objs.dec(1);
	size_t *size = (size_t *) (((char *) p) - sizeof(size_t));
	alloc_bytes.dec(*size);
	free(size);
}

void *operator new[](size_t n) throw (std::bad_alloc)
{
	alloc_objs.inc(1);
	alloc_bytes.inc(n);
	void *p = malloc(n + sizeof(size_t));
	size_t *size = (size_t *) p;
	*size = n;

	return (void *) (size + 1);
}

void operator delete[](void *p) throw ()
{
	if (p == NULL)
		return;

	alloc_objs.dec(1);
	size_t *size = (size_t *) (((char *) p) - sizeof(size_t));
	alloc_bytes.dec(*size);
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
