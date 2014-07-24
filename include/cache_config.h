#ifndef __CACHE_CONFIG_H__
#define __CACHE_CONFIG_H__

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
#include <tr1/unordered_map>

#include "common.h"
#include "cache.h"
#include "thread.h"

class page_cache;

enum {
	TREE_CACHE,
	ASSOCIATIVE_CACHE,
	HASH_INDEX_CACHE,
	CUCKOO_CACHE,
	LRU2Q_CACHE,
	GCLOCK_CACHE,
};

/**
 * This class defines the information about the cache.
 * It defines
 *	size,
 *	type,
 *	cache distribution: how cache is split and distributed in NUMA nodes;
 *	page mapping: how to map a page to a partition of the cache.
 */
class cache_config
{
	long size;
	int type;
	// node id <-> the size of each partition
	std::tr1::unordered_map<int, long> part_sizes;

	class create_cache_thread: public thread
	{
		const cache_config &config;
		int max_num_pending_flush;
		page_cache *cache;
	public:
		create_cache_thread(const cache_config &_config, int max_num_pending_flush,
				std::string name, int node_id): thread(name,
					node_id), config(_config) {
			this->max_num_pending_flush = max_num_pending_flush;
		}

		void run() {
			cache = config.__create_cache_on_node(get_node_id(),
					max_num_pending_flush);
			stop();
		}

		page_cache *get_cache() {
			return cache;
		}
	};

	page_cache *__create_cache_on_node(int node_id,
			int max_num_pending_flush) const;
protected:
	void init(const std::tr1::unordered_map<int, long> &part_sizes) {
		this->part_sizes = part_sizes;
	}

public:
	cache_config(long size, int type) {
		this->size = size;
		this->type = type;
	}

	virtual ~cache_config() {
	}

	long get_size() const {
		return size;
	}

	int get_type() const {
		return type;
	}

	int get_num_cache_parts() const {
		return (int) part_sizes.size();
	}

	long get_part_size(int node_id) const {
		std::tr1::unordered_map<int, long>::const_iterator it
			= part_sizes.find(node_id);
		assert(it != part_sizes.end());
		return it->second;
	}

	virtual int page2cache(const page_id_t &pg_id) const = 0;

	void get_node_ids(std::vector<int> &node_ids) const {
		for (std::tr1::unordered_map<int, long>::const_iterator it
				= part_sizes.begin(); it != part_sizes.end(); it++)
			node_ids.push_back(it->first);
	}

	page_cache *create_cache_on_node(int node_id, int max_num_pending_flush) const;
	int create_cache_on_nodes(const std::vector<int> &node_ids,
			int max_num_pending_flush, std::vector<page_cache *> &caches) const;
	void destroy_cache_on_node(page_cache *cache) const;
	page_cache *create_cache(int max_num_pending_flush) const;
	void destroy_cache(page_cache *cache) const;

};

/**
 * This cache config split a global cache and distributes them evenly into nodes.
 */
class even_cache_config: public cache_config
{
public:
	even_cache_config(long size, int type,
			const std::vector<int> &node_ids): cache_config(size, type) {
		long part_size = size / node_ids.size();
		std::tr1::unordered_map<int, long> part_sizes;
		for (size_t i = 0; i < node_ids.size(); i++)
			part_sizes.insert(std::pair<int, long>(node_ids[i], part_size));
		init(part_sizes);
	}

	virtual int page2cache(const page_id_t &pg_id) const {
		return (int) (pg_id.get_offset() / PAGE_SIZE) % get_num_cache_parts();
	}
};

/**
 * This cache config is for testing.
 * It actually can yield the best performance for parted global cache in the NUMA machine.
 */
class test_cache_config: public cache_config
{
	std::vector<int> map;
	int num_page_sets;
	int part_size;		// also in the number of page sets.
public:
	test_cache_config(long size, int type,
			const std::vector<int> &node_ids): cache_config(size, type), map(3) {
		num_page_sets = size / PAGE_SIZE / CELL_MIN_NUM_PAGES;
		part_size = num_page_sets / 3;
		map[0] = 0;
		map[1] = 2;
		map[2] = 3;
		assert(node_ids.size() == 4);
		std::tr1::unordered_map<int, long> part_sizes;
		part_sizes.insert(std::pair<int, long>(0, size));
		part_sizes.insert(std::pair<int, long>(1, 0));
		part_sizes.insert(std::pair<int, long>(2, 0));
		part_sizes.insert(std::pair<int, long>(3, 0));
		init(part_sizes);
	}

	virtual int page2cache(const page_id_t &pg_id) const {
		return 0;
//		return map[(off / PAGE_SIZE) % 3];
//		return (int) (off / PAGE_SIZE) % 2;
//		int idx = (int) (off / PAGE_SIZE);
//		if (idx % 7 == 0)
//			return 0;
//		else
//			return idx % 7 % 2 + 2;
	}
};

class file_mapper;

/**
 * This cache config split and distribute a global cache according to
 * the distribution of storage devices on the machine.
 */
class file_map_cache_config: public cache_config
{
	file_mapper *mapper;
	int shift;
public:
	file_map_cache_config(long size, int type, const std::vector<int> &node_ids,
			file_mapper *mapper, int shift = 0);

	virtual int page2cache(const page_id_t &pg_id) const;
};

#endif
