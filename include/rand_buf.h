#ifndef __RAND_BUF_H__
#define __RAND_BUF_H__

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
