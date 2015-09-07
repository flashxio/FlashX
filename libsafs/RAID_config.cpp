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

#include <boost/format.hpp>

#include "log.h"
#include "RAID_config.h"
#include "file_mapper.h"
#include "native_file.h"

namespace safs
{

static safs_header get_safs_header(const RAID_config &conf,
		const std::string &file_name)
{
	safs_file f(conf, file_name);
	return f.get_header();
}

file_mapper *RAID_config::create_file_mapper(const std::string &file_name) const
{
	/*
	 * The individual files on the native file system are partitions of
	 * a logical SAFS file. They are organized as follows:
	 * in each SSD, there is a directory named after the SAFS file name;
	 * inside the directory, there is exactly one file that stores the data
	 * of a partition, and the file name is the partition ID.
	 */
	std::map<int, part_file_info> file_map;
	for (unsigned i = 0; i < root_paths.size(); i++) {
		std::string dir_name = root_paths[i].get_file_name()
			+ std::string("/") + file_name;
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
		if (part_ids.size() > 1)
			part_ids = safs_file::erase_header_file(part_ids);
		if (part_ids.size() != 1) {
			fprintf(stderr,
					"wrong format of the SAFS file %s, check the directory %s\n",
					file_name.c_str(), dir_name.c_str());
			return NULL;
		}
		int part_id = atoi(part_ids[0].c_str());
		part_file_info info(dir_name + std::string("/") + part_ids[0],
				root_paths[i].get_disk_id(), root_paths[i].get_node_id());
		file_map.insert(std::pair<int, part_file_info>(part_id, info));
	}
	if (file_map.size() < root_paths.size()) {
		fprintf(stderr, "duplicated partition id of the SAFS file %s\n",
				file_name.c_str());
		return NULL;
	}

	std::vector<part_file_info> files;
	for (std::map<int, part_file_info>::const_iterator it = file_map.begin();
			it != file_map.end(); it++) {
		files.push_back(it->second);
	}

	safs_header header = get_safs_header(*this, file_name);
	int block_size = RAID_block_size;
	int mapping_option = RAID_mapping_option;
	// The per-file config can overwrite the default config.
	if (header.is_valid()) {
		block_size = header.get_block_size();
		mapping_option = header.get_mapping_option();
	}
	switch (mapping_option) {
		case RAID0:
			return new RAID0_mapper(file_name, files, block_size);
		case RAID5:
			return new RAID5_mapper(file_name, files, block_size);
		case HASH:
			return new hash_mapper(file_name, files, block_size);
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
		node_ids.insert(root_paths[k].get_node_id());
	}
	return node_ids;
}

static int retrieve_data_files(std::string file_file,
		std::vector<part_file_info> &data_files)
{
	char *line = NULL;
	size_t size = 0;
	int line_length;
	FILE *fd = fopen(file_file.c_str(), "r");
	if (fd == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format("open RAID conf file %1%: %2%")
			% file_file % strerror(errno);
		return 0;
	}

	int disk_id = 0;
	while ((line_length = getline(&line, &size, fd)) > 0) {
		line[line_length - 1] = 0;
		// skip comment lines.
		if (*line == '#')
			continue;

		char *colon = strstr(line, ":");
		char *name = line;
		int node_id = 0;
		if (colon) {
			*colon = 0;
			std::string node_id_str = line;
			node_id_str.erase(std::remove_if(node_id_str.begin(),
						node_id_str.end(), isspace), node_id_str.end());
			node_id = atoi(node_id_str.c_str());
			colon++;
			name = colon;
		}
		std::string path_name = name;
		path_name.erase(std::remove_if(path_name.begin(), path_name.end(),
					isspace), path_name.end());
		data_files.emplace_back(path_name, disk_id, node_id);
		free(line);
		line = NULL;
		size = 0;
		disk_id++;
	}
	fclose(fd);
	return data_files.size();
}

RAID_config::ptr RAID_config::create(const std::string &conf_file,
		int mapping_option, int block_size)
{
	RAID_config::ptr conf = RAID_config::ptr(new RAID_config());
	conf->conf_file = conf_file;
	int ret = retrieve_data_files(conf->conf_file, conf->root_paths);
	if (ret == 0)
		return RAID_config::ptr();

	conf->RAID_mapping_option = mapping_option;
	conf->RAID_block_size = block_size;
	return conf;
}

}
