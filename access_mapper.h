#ifndef __ACCESS_MAPPER_H__
#define __ACCESS_MAPPER_H__

#include <vector>

#include "common.h"

class access_mapper
{
public:
	virtual int map(off_t off) const = 0;
	virtual access_mapper *clone() const = 0;
};

/**
 * This is to define a mapping function to map an access to
 * a particular NUMA node.
 */
class NUMA_access_mapper: public access_mapper
{
	std::vector<int> mapping_array;

	NUMA_access_mapper() {
	}
public:
	NUMA_access_mapper(const std::vector<file_info> &files): mapping_array(
			files.size()) {
		int size = (int) files.size();
		for (int i = 0; i < size; i++)
			mapping_array[i] = files[i].node_id;
	}

	int map(off_t off) const {
		int size = (int) mapping_array.size();
		int idx = file_hash(off / PAGE_SIZE, size);
		return mapping_array[idx];
	}

	access_mapper *clone() const {
		NUMA_access_mapper *mapper = new NUMA_access_mapper();
		mapper->mapping_array = this->mapping_array;
		return mapper;
	}
};

#endif
