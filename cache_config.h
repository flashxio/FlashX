#ifndef __CACHE_CONFIG_H__
#define __CACHE_CONFIG_H__

#include <vector>

#include "common.h"
#include "file_mapper.h"

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
protected:
	void init(const std::tr1::unordered_map<int, long> &part_sizes) {
		this->part_sizes = part_sizes;
	}

public:
	cache_config(long size, int type) {
		this->size = size;
		this->type = type;
	}

	long get_size() const {
		return size;
	}

	int get_type() const {
		return type;
	}

	int get_num_caches() const {
		return (int) part_sizes.size();
	}

	long get_part_size(int node_id) const {
		std::tr1::unordered_map<int, long>::const_iterator it
			= part_sizes.find(node_id);
		fprintf(stderr, "node_id: %d\n", node_id);
		assert(it != part_sizes.end());
		return it->second;
	}

	virtual int page2cache(off_t off) const = 0;

	void get_node_ids(std::vector<int> &node_ids) const {
		for (std::tr1::unordered_map<int, long>::const_iterator it
				= part_sizes.begin(); it != part_sizes.end(); it++)
			node_ids.push_back(it->first);
	}

	page_cache *create_cache_on_node(int node_id) const;

	page_cache *create_cache() const;

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

	virtual int page2cache(off_t off) const {
		return (int) (off / PAGE_SIZE) % get_num_caches();
	}
};

/**
 * This cache config is for testing.
 * It actually can yield the best performance for parted global cache in the NUMA machine.
 */
class test_cache_config: public cache_config
{
public:
	test_cache_config(long size, int type,
			const std::vector<int> &node_ids): cache_config(size, type) {
		assert(node_ids.size() == 4);
		std::tr1::unordered_map<int, long> part_sizes;
		part_sizes.insert(std::pair<int, long>(0, 0));
		part_sizes.insert(std::pair<int, long>(1, 0));
		part_sizes.insert(std::pair<int, long>(2, size / 2));
		part_sizes.insert(std::pair<int, long>(3, size / 2));
		init(part_sizes);
	}

	virtual int page2cache(off_t off) const {
		return (int) (off / PAGE_SIZE) % 2 + 2;
	}
};

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

	virtual int page2cache(off_t off) const {
		int idx = mapper->map2file(off / PAGE_SIZE);
		return mapper->get_file_node_id(idx) + shift;
	}
};

#endif
