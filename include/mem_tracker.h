#ifndef __MEM_TRACKER_H__
#define __MEM_TRACKER_H__

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
