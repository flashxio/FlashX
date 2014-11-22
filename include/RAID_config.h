#ifndef __RAID_CONFIG_H__
#define __RAID_CONFIG_H__

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
#include <set>

#include "safs_file.h"

namespace safs
{

class file_mapper;

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
	std::vector<part_file_info> root_paths;
public:
	typedef std::shared_ptr<RAID_config> ptr;

	static ptr create(const std::string &conf_file, int mapping_option,
			int block_size);

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

	const part_file_info &get_disk(int idx) const {
		return root_paths[idx];
	}

	int get_disks(std::vector<part_file_info> &disks) const {
		disks = root_paths;
		return disks.size();
	}

	int get_num_disks() const {
		return root_paths.size();
	}
};

}

#endif
