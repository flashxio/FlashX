#ifndef __NUMA_CACHE_H__
#define __NUMA_CACHE_H__

#include <tr1/unordered_map>
#include <vector>

#include "cache.h"
#include "associative_cache.h"
#include "hash_index_cache.h"
#include "LRU2Q.h"

/**
 * This cache divides a global cache to pieces of the same size,
 * and place them on the specified NUMA nodes.
 */
class NUMA_cache: public page_cache
{
	long cache_size;
	std::vector<page_cache *> caches;
public:
	NUMA_cache(int cache_type, long cache_size,
			const std::vector<int> &node_ids) {
		this->cache_size = cache_size;
		for (size_t i = 0; i < node_ids.size(); i++) {
			page_cache *cache = create_cache(cache_type,
					cache_size / node_ids.size(), node_ids[i], 1);
			caches.push_back(cache);
		}
	}

	static page_cache *create_cache(int cache_type, long cache_size,
			int node_id, int offset_factor) {
		page_cache *cache;
		struct bitmask *orig_bmp = numa_get_membind();
		bind_mem2node_id(node_id);
		switch (cache_type) {
			case LRU2Q_CACHE:
				cache = new LRU2Q_cache(cache_size);
				break;
			case ASSOCIATIVE_CACHE:
				cache = new associative_cache(cache_size, MAX_CACHE_SIZE,
						node_id, offset_factor);
				break;
			case HASH_INDEX_CACHE:
				cache = new hash_index_cache(cache_size, node_id);
				break;
			default:
				fprintf(stderr, "wrong cache type\n");
				exit(1);
		}
		numa_set_membind(orig_bmp);
		return cache;
	}

	virtual page *search(off_t offset, off_t &old_off) {
		int idx = offset / PAGE_SIZE % caches.size();
		return caches[idx]->search(offset, old_off);
	}

	virtual page *search(off_t offset) {
		int idx = offset / PAGE_SIZE % caches.size();
		return caches[idx]->search(offset);
	}

	virtual long size() {
		return cache_size;
	}

	virtual void init(io_interface *underlying) {
		for (size_t i = 0; i < caches.size(); i++)
			caches[i]->init(underlying);
	}

	virtual void print_stat() const {
		for (size_t i = 0; i < caches.size(); i++)
			caches[i]->print_stat();
	}
};

#endif
