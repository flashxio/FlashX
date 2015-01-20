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

#include "memory_manager.h"

namespace safs
{

const long SHRINK_NPAGES = 1024;
const long INCREASE_SIZE = 1024 * 1024 * 128;

memory_manager::memory_manager(
		long max_size, int node_id): slab_allocator(
			std::string("mem_manager-") + itoa(node_id), PAGE_SIZE,
			INCREASE_SIZE <= max_size ? INCREASE_SIZE : max_size,
			// We don't initialize pages but we pin pages.
			max_size, node_id, false, true) {
}

/**
 * get `npages' pages for `request_cache'.
 * In the case of shrinking caches, it makes no sense
 * to shrink the cache that is requesting free pages.
 */
bool memory_manager::get_free_pages(int npages,
		char **pages, page_cache *request_cache) {
	int ret = slab_allocator::alloc(pages, npages);
	/* 
	 * slab_allocator allocates either all required number of 
	 * pages or none of them.
	 */
	if (ret == 0) {
		long size = 0;
		page_cache *cache = NULL;
		/* shrink the cache of the largest size */
		for (unsigned int i = 0; i < caches.size(); i++) {
			if (size < caches[i]->size()) {
				size = caches[i]->size();
				cache = caches[i];
			}
		}
		if (cache == NULL)
			return false;
		/* 
		 * if we are going to shrink the cache that requests
		 * free pages, it just fails.
		 */
		if (request_cache == cache) {
			return false;
		}
		int num_shrink = SHRINK_NPAGES;
		if (num_shrink < npages)
			num_shrink = npages;
		char *buf[num_shrink];
		if (!cache->shrink(num_shrink, buf)) {
			return false;
		}
		slab_allocator::free(buf, num_shrink);
		/* now it's guaranteed that we have enough free pages. */
		ret = slab_allocator::alloc(pages, npages);
	}
	return true;
}

void memory_manager::free_pages(int npages, char **pages) {
	slab_allocator::free(pages, npages);
}

}
