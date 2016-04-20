#ifndef __ALIGNED_ALLOCATOR_H__
#define __ALIGNED_ALLOCATOR_H__

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

#include <malloc.h>

#include "common.h"

/*
 * The class can allocate aligned memory if the alignment provided
 * by the user is a power of two.
 */
class aligned_allocator
{
	size_t alignment;
public:
	aligned_allocator(size_t alignment) {
		this->alignment = alignment;
		bool aligned = align_check(alignment);
		printf("aligned allocator: alignment: %ld, aligned: %d\n",
				this->alignment, aligned);
		// If we find it's not a power of 2, indicate it's not aligned.
		if (!aligned)
			this->alignment = 0;
	}

	void *alloc(size_t size) {
		if (alignment) {
			void *buf = NULL;
			int ret = posix_memalign(&buf, alignment, size);
			if (ret == 0)
				return buf;
			else
				return NULL;
		}
		else
			return malloc(size);
	}

	void dealloc(void *p) {
		free(p);
	}
};

#endif
