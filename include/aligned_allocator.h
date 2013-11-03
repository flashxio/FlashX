#ifndef __ALIGNED_ALLOCATOR_H__
#define __ALIGNED_ALLOCATOR_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <malloc.h>

#include "common.h"

/**
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
			return memalign(alignment, size);
		}
		else
			return malloc(size);
	}

	void dealloc(void *p) {
		free(p);
	}
};

#endif
