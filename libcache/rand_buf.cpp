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

#include "rand_buf.h"

rand_buf::rand_buf(int buf_size, int entry_size, int nodeid)
#ifdef MEMCHECK
	: allocator(entry_size)
#else
	// We don't initialize pages but pin them.
	: allocator(std::string("rand_buf-") + itoa(nodeid), entry_size,
			buf_size, buf_size, nodeid, false, true, 0)
#endif
{
	this->entry_size = entry_size;
}

char *rand_buf::next_entry(int size) {
#ifdef MEMCHECK
	return (char *) allocator.alloc(size);
#else
	if (size > entry_size)
		return (char *) valloc(size);
	return allocator.alloc();
#if 0
	pthread_spin_lock(&lock);
	assert(!free_refs.is_empty());
	int off = free_refs.pop_front();
	assert(marks[off / entry_size] == 0);
	marks[off / entry_size] = 1;
	char *ret = &buf[off];
	pthread_spin_unlock(&lock);
	return ret;
#endif
#endif
}

void rand_buf::free_entry(char *buf) {
#ifdef MEMCHECK
	allocator.dealloc(buf);
#else
	if (allocator.contains(buf))
		allocator.free(buf);
	else
		free(buf);
#if 0
	pthread_spin_lock(&lock);
	off_t off = buf - this->buf;
	free_refs.push_back(off);
	off /= entry_size;
	if (marks[off] == 0)
		printf("%ld: free %p error\n", pthread_self(), buf);
	assert(marks[off]);
	marks[off] = 0;
	pthread_spin_unlock(&lock);
#endif
#endif
}
