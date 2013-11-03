#ifndef __MEM_TRACKER_H__
#define __MEM_TRACKER_H__

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

#include <stdlib.h>

#ifdef ENABLE_MEM_TRACE

#include <new>

void init_mem_tracker();

void *operator new(size_t n) throw (std::bad_alloc);
void operator delete(void *p) throw ();
void *operator new[](size_t n) throw (std::bad_alloc);
void operator delete[](void *p) throw ();

#endif

size_t get_alloc_objs();
size_t get_alloc_bytes();

#endif
