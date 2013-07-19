#ifndef __RAID_CONFIG_H__
#define __RAID_CONFIG_H__

#include <vector>
#include <set>

#include "file_mapper.h"

enum {
	RAID0,
	RAID5,
	HASH,
};

class RAID_config
{
	int RAID_mapping_option;
	int RAID_block_size;
	std::vector<file_info> files;
public:
	RAID_config(std::vector<file_info> &files, int mapping_option, int block_size) {
		RAID_mapping_option = mapping_option;
		RAID_block_size = block_size;
		this->files = files;
	}

	file_mapper *create_file_mapper() const;

	/**
	 * This returns the nodes where the RAID attaches to.
	 */
	std::set<int> get_node_ids() const;

	const file_info &get_file(int idx) const {
		return files[idx];
	}
};

#endif
