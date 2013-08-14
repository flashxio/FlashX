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

	memory_manager(long max_size, int node_id);

	~memory_manager() {
		// TODO
	}
public:
	static memory_manager *create(long max_size, int node_id) {
		assert(node_id >= 0);
		void *addr = numa_alloc_onnode(sizeof(memory_manager), node_id);
		return new(addr) memory_manager(max_size, node_id);
	}

	static void destroy(memory_manager *m) {
		m->~memory_manager();
		numa_free(m, sizeof(*m));
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
