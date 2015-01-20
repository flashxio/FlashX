#ifndef __NUMA_CACHE_H__
#define __NUMA_CACHE_H__

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

#include <tr1/unordered_map>
#include <vector>

#include "cache.h"
#include "cache_config.h"

namespace safs
{

/**
 * This cache divides a global cache to pieces of the same size,
 * and place them on the specified NUMA nodes.
 */
class NUMA_cache: public page_cache
{
	const cache_config *cache_conf;
	std::vector<page_cache::ptr> caches;

	NUMA_cache(const cache_config *config, int max_num_pending_flush): caches(
			config->get_num_cache_parts()) {
		cache_conf = config;
		std::vector<int> node_ids;
		cache_conf->get_node_ids(node_ids);
		cache_conf->create_cache_on_nodes(node_ids,
				max_num_pending_flush / node_ids.size(), caches);
	}
public:
	static page_cache::ptr create(const cache_config *config,
			int max_num_pending_flush) {
		return page_cache::ptr(new NUMA_cache(config, max_num_pending_flush));
	}

	page_cache &get_cache_on_node(int node_id) const {
		return *caches[node_id];
	}

	virtual page *search(const page_id_t &pg_id, page_id_t &old_id) {
		int idx = cache_conf->page2cache(pg_id);
		return caches[idx]->search(pg_id, old_id);
	}

	virtual page *search(const page_id_t &pg_id) {
		int idx = cache_conf->page2cache(pg_id);
		return caches[idx]->search(pg_id);
	}

	virtual long size() {
		return cache_conf->get_size();
	}

	// TODO shouldn't I use a different underlying IO for cache
	// on the different nodes.
	virtual void init(std::shared_ptr<io_interface> underlying) {
		for (size_t i = 0; i < caches.size(); i++) {
			// TODO this is a ugly hack. It makes sure that each flush
			// thread of a SA-cache has a reference to this NUMA cache.
			caches[i]->create_flusher(underlying, this);
			caches[i]->init(underlying);
		}
	}

	virtual void print_stat() const {
		for (size_t i = 0; i < caches.size(); i++)
			caches[i]->print_stat();
	}

	virtual void mark_dirty_pages(thread_safe_page *pages[], int num,
			io_interface &io) {
		for (int i = 0; i < num; i++) {
			page_id_t pg_id(pages[i]->get_file_id(), pages[i]->get_offset());
			int idx = cache_conf->page2cache(pg_id);
			caches[idx]->mark_dirty_pages(&pages[i], 1, io);
		}
	}

	virtual void flush_callback(io_request &req) {
		assert(req.within_1page());
		page_id_t pg_id(req.get_file_id(), req.get_offset());
		int idx = cache_conf->page2cache(pg_id);
		caches[idx]->flush_callback(req);
	}

	virtual int flush_dirty_pages(page_filter *filter, int max_num) {
		int tot = 0;
		for (size_t i = 0; i < caches.size(); i++) {
			tot += caches[i]->flush_dirty_pages(filter,
					max_num / caches.size());
		}
		return tot;
	}

	virtual void sanity_check() const {
		for (size_t i = 0; i < caches.size(); i++) {
			caches[i]->sanity_check();
		}
	}
};

}

#endif
