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

#include "RAID_config.h"
#include "file_mapper.h"
#include "native_file.h"

file_mapper *RAID_config::create_file_mapper(const std::string &file_name) const
{
	/*
	 * The individual files on the native file system are partitions of
	 * a logical SAFS file. They are organized as follows:
	 * in each SSD, there is a directory named after the SAFS file name;
	 * inside the directory, there is exactly one file that stores the data
	 * of a partition, and the file name is the partition ID.
	 */
	std::map<int, file_info> file_map;
	for (unsigned i = 0; i < root_paths.size(); i++) {
		std::string dir_name = root_paths[i].name + std::string("/") + file_name;
		native_dir dir(dir_name);
		if (!dir.exist()) {
			fprintf(stderr, "%s for the SAFS file %s doesn't exist\n",
					dir_name.c_str(), file_name.c_str());
			return NULL;
		}
		if (!dir.is_dir()) {
			fprintf(stderr, "%s for the SAFS file %s isn't a directory\n",
					dir_name.c_str(), file_name.c_str());
			return NULL;
		}
		std::vector<std::string> part_ids;
		dir.read_all_files(part_ids);
		printf("there are %ld files in %s\n", part_ids.size(), dir_name.c_str());
		if (part_ids.size() != 1) {
			fprintf(stderr,
					"wrong format of the SAFS file %s, check the directory %s\n",
					file_name.c_str(), dir_name.c_str());
			return NULL;
		}
		int part_id = atoi(part_ids[0].c_str());
		file_info info = root_paths[i];
		info.name = dir_name + std::string("/") + part_ids[0];
		file_map.insert(std::pair<int, file_info>(part_id, info));
	}
	if (file_map.size() < root_paths.size()) {
		fprintf(stderr, "duplicated partition id of the SAFS file %s\n",
				file_name.c_str());
		return NULL;
	}

	std::vector<file_info> files;
	for (std::map<int, file_info>::const_iterator it = file_map.begin();
			it != file_map.end(); it++) {
		files.push_back(it->second);
	}

	switch (RAID_mapping_option) {
		case RAID0:
			return new RAID0_mapper(file_name, files, RAID_block_size);
		case RAID5:
			return new RAID5_mapper(file_name, files, RAID_block_size);
		case HASH:
			return new hash_mapper(file_name, files, RAID_block_size);
		default:
			fprintf(stderr, "wrong RAID mapping option\n");
			exit(1);
	}
}

file_mapper *RAID_config::create_file_mapper() const
{
	switch (RAID_mapping_option) {
		case RAID0:
			return new RAID0_mapper("root", root_paths, RAID_block_size);
		case RAID5:
			return new RAID5_mapper("root", root_paths, RAID_block_size);
		case HASH:
			return new hash_mapper("root", root_paths, RAID_block_size);
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
