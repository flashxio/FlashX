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
	std::string conf_file;
	int RAID_mapping_option;
	int RAID_block_size;
	std::vector<file_info> root_paths;
public:
	RAID_config() {
	}

	RAID_config(const std::string &conf_file, int mapping_option, int block_size) {
		this->conf_file = conf_file;
		retrieve_data_files(conf_file, root_paths);
		printf("there are %ld disks\n", root_paths.size());

		RAID_mapping_option = mapping_option;
		RAID_block_size = block_size;
	}

	file_mapper *create_file_mapper(const std::string &file_name) const;

	/**
	 * This returns the nodes where the RAID attaches to.
	 */
	std::set<int> get_node_ids() const;

	const file_info &get_disk(int idx) const {
		return root_paths[idx];
	}

	int get_num_disks() const {
		return root_paths.size();
	}
};

#endif
