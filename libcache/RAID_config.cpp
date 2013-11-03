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

#include "RAID_config.h"
#include "file_mapper.h"

file_mapper *RAID_config::create_file_mapper(const std::string &file_name) const
{
	std::vector<file_info> files = root_paths;
	for (unsigned i = 0; i < files.size(); i++) {
		files[i].name += std::string("/") + file_name;
	}
	switch (RAID_mapping_option) {
		case RAID0:
			return new RAID0_mapper(files, RAID_block_size);
		case RAID5:
			return new RAID5_mapper(files, RAID_block_size);
		case HASH:
			return new hash_mapper(files, RAID_block_size);
		default:
			fprintf(stderr, "wrong RAID mapping option\n");
			exit(1);
	}
}

std::set<int> RAID_config::get_node_ids() const
{
	std::set<int> node_ids;
	int num_paths = root_paths.size();
	for (int k = 0; k < num_paths; k++) {
		node_ids.insert(root_paths[k].node_id);
	}
	return node_ids;
}
