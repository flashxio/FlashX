#ifndef __RAID_CONFIG_H__
#define __RAID_CONFIG_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

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

	/**
	 * Create a file mapper for the RAID directories.
	 */
	file_mapper *create_file_mapper() const;
	/**
	 * Create a file mapper for a file in the RAID.
	 */
	file_mapper *create_file_mapper(const std::string &file_name) const;

	/**
	 * This returns the nodes where the RAID attaches to.
	 */
	std::set<int> get_node_ids() const;

	const file_info &get_disk(int idx) const {
		return root_paths[idx];
	}

	int get_disks(std::vector<file_info> &disks) const {
		disks = root_paths;
		return disks.size();
	}

	int get_num_disks() const {
		return root_paths.size();
	}
};

#endif
