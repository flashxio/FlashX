#ifndef __RAND_BUF_H__
#define __RAND_BUF_H__

#include "container.h"
#include "aligned_allocator.h"
#include "slab_allocator.h"

/**
 * The class maintains a set of buffers of a mixed size for a single thread.
 * It works the same as slab_allocator.
 * To reduce overhead, the class isn't thread-safe.
 */
class rand_buf
{
	int entry_size;
#ifdef MEMCHECK
	aligned_allocator allocator;
#else
	/* where the data read from the disk is stored */
	slab_allocator allocator;
#endif
public:
	rand_buf(int buf_size, int entry_size, int nodeid);

	void free_entry(char *buf);

	char *next_entry(int size);

	int get_entry_size() {
		return entry_size;
	}
};

#endif
