#ifndef __ALIGNED_ALLOCATOR_H__
#define __ALIGNED_ALLOCATOR_H__

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
