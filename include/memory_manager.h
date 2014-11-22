#ifndef __MEMORY_MANAGER__
#define __MEMORY_MANAGER__

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

#include <vector>

#include "cache.h"
#include "slab_allocator.h"

namespace safs
{

/**
 * manage free pages in the cache.
 * It also allocates pages from the operating system.
 */
class memory_manager: public slab_allocator
{
	std::vector<page_cache *> caches;

	memory_manager(long max_size, int node_id);

	~memory_manager() {
		// TODO
	}
public:
	static memory_manager *create(long max_size, int node_id) {
		assert(node_id >= 0);
		return new memory_manager(max_size, node_id);
	}

	static void destroy(memory_manager *m) {
		delete m;
	}

	void register_cache(page_cache *cache) {
		caches.push_back(cache);
	}

	void unregister_cache(page_cache *cache) {
		// TODO
	}

	bool get_free_pages(int npages, char **pages, page_cache *cache);
	void free_pages(int npages, char **pages);

	long average_cache_size() {
		return get_max_size() / caches.size();
	}
};

}

#endif
