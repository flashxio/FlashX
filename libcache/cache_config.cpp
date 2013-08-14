#include <tr1/unordered_map>

#include "cache_config.h"
#include "associative_cache.h"
#include "hash_index_cache.h"
#include "LRU2Q.h"
#include "NUMA_cache.h"

page_cache *cache_config::create_cache_on_node(int node_id) const
{
	page_cache *cache;
	switch (get_type()) {
		case LRU2Q_CACHE:
			cache = new LRU2Q_cache(get_part_size(node_id));
			break;
		case ASSOCIATIVE_CACHE:
			cache = associative_cache::create(get_part_size(node_id),
					MAX_CACHE_SIZE, node_id, 1);
			break;
		case HASH_INDEX_CACHE:
			cache = new hash_index_cache(get_part_size(node_id), node_id);
			break;
		default:
			fprintf(stderr, "wrong cache type\n");
			exit(1);
	}
	return cache;
}

page_cache *cache_config::create_cache() const
{
	std::vector<int> node_ids;
	get_node_ids(node_ids);
	if (node_ids.size() == 1)
		return create_cache_on_node(node_ids[0]);
	else
		return new NUMA_cache(this);
}

static bool node_exist(const std::vector<int> &node_ids, int node_id)
{
	for (size_t i = 0; i < node_ids.size(); i++)
		if (node_ids[i] == node_id)
			return true;
	return false;
}

file_map_cache_config::file_map_cache_config(long size, int type,
		const std::vector<int> &node_ids,
		file_mapper *mapper, int shift): cache_config(size, type)
{
	this->mapper = mapper;
	this->shift = shift;
	// This counts the number of files connected to each node.
	std::map<int, int> node_files;
	for (int i = 0; i < mapper->get_num_files(); i++) {
		int node_id = mapper->get_file_node_id(i);
		std::map<int, int>::iterator it = node_files.find(node_id);
		if (it == node_files.end())
			node_files.insert(std::pair<int, int>(node_id, 1));
		else
			it->second++;
	}

	std::tr1::unordered_map<int, long> part_sizes;
	int tot_files = 0;
	for (size_t i = 0; i < node_ids.size(); i++) {
		int node_id = node_ids[i];
		std::map<int, int>::const_iterator it = node_files.find(node_id);
		int new_node_id = (node_id + shift) % node_ids.size();
		if (it == node_files.end()) {
			part_sizes.insert(std::pair<int, long>(new_node_id, 0));
			printf("file mapping: cache part %d: size: %d\n", new_node_id, 0);
		}
		else {
			int num_files = it->second;
			tot_files += num_files;
			long part_size = size * (((float) num_files)
					/ mapper->get_num_files());
			assert(node_exist(node_ids, new_node_id));
			part_sizes.insert(std::pair<int, long>(new_node_id, part_size));
			printf("file mapping: cache part %d: size: %ld\n",
					new_node_id, part_size);
		}
	}
	assert(tot_files == mapper->get_num_files());
	init(part_sizes);
}
