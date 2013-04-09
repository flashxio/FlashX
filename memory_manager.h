#ifndef __MEMORY_MANAGER__
#define __MEMORY_MANAGER__

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
public:
	memory_manager(long max_size);

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
