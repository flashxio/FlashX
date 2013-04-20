#ifndef __NUMA_CACHE_H__
#define __NUMA_CACHE_H__

#include <tr1/unordered_map>
#include <vector>

#include "cache.h"
#include "cache_config.h"

/**
 * This cache divides a global cache to pieces of the same size,
 * and place them on the specified NUMA nodes.
 */
class NUMA_cache: public page_cache
{
	const cache_config *cache_conf;
	std::vector<page_cache *> caches;
public:
	NUMA_cache(const cache_config *config): caches(config->get_num_caches()) {
		cache_conf = config;
		std::vector<int> node_ids;
		cache_conf->get_node_ids(node_ids);
		for (size_t i = 0; i < node_ids.size(); i++) {
			page_cache *cache = cache_conf->create_cache_on_node(node_ids[i]);
			caches[i] = cache;
		}
	}

	virtual page *search(off_t offset, off_t &old_off) {
		int idx = cache_conf->page2cache(offset);
		return caches[idx]->search(offset, old_off);
	}

	virtual page *search(off_t offset) {
		int idx = cache_conf->page2cache(offset);
		return caches[idx]->search(offset);
	}

	virtual long size() {
		return cache_conf->get_size();
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
