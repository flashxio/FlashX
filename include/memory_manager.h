#ifndef __MEMORY_MANAGER__
#define __MEMORY_MANAGER__

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

#include <vector>

#include "cache.h"
#include "slab_allocator.h"

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

#endif
