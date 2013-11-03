#ifndef __MEM_TRACKER_H__
#define __MEM_TRACKER_H__

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
