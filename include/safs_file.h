#ifndef __SAFS_FILE_H__
#define __SAFS_FILE_H__

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

#include <stdlib.h>

#include <string>
#include <vector>
#include <set>

#include "common.h"
#include "native_file.h"
#include "safs_header.h"
#include "parameters.h"

namespace safs
{

/*
 * This stores the local filesystem file information that stores
 * each partition of SAFS files.
 */
class part_file_info
{
	// The file name in the local filesystem.
	std::string name;
	// The disk ID where the block is stored.
	int disk_id;
	// The NUMA node id where the disk is connected to.
	int node_id;
public:
	part_file_info() {
		disk_id = 0;
		node_id = 0;
	}

	part_file_info(const std::string &name, int disk_id, int node_id) {
		this->name = name;
		this->disk_id = disk_id;
		this->node_id = node_id;
	}

	std::string get_file_name() const {
		return name;
	}

	int get_disk_id() const {
		return disk_id;
	}

	int get_node_id() const {
		return node_id;
	}
};

class RAID_config;
class safs_file_group;

class safs_file
{
	// The collection of native files.
	std::vector<part_file_info> native_dirs;
	// The path of the header file in the local Linux filesystem.
	std::string header_file;
	std::string name;

	std::string get_header_file() const;
public:
	static std::vector<std::string> erase_header_file(
			const std::vector<std::string> &files);

	safs_file(const RAID_config &conf, const std::string &file_name);

	safs_header get_header() const;

	/*
	 * An SAFS file allows a user to store user-defined metadata along with
	 * the data in the file.
	 */
	bool set_user_metadata(const std::vector<char> &data);
	std::vector<char> get_user_metadata() const;

	const std::string &get_name() const {
		return name;
	}

	bool exist() const;
	ssize_t get_size() const;
	bool create_file(size_t file_size,
			int block_size = params.get_RAID_block_size(),
			int mapping_option = params.get_RAID_mapping_option(),
			std::shared_ptr<safs_file_group> group = NULL);
	bool delete_file();
	bool rename(const std::string &new_name);
};

class safs_file_group
{
public:
	enum group_t {
		NAIVE,
		ROTATE,
		RAND_ROTATE,
	};
	typedef std::shared_ptr<safs_file_group> ptr;

	static ptr create(const RAID_config &conf, group_t type);
	virtual std::vector<int> add_file(safs_file &file) = 0;
	virtual std::string get_name() const = 0;
};

size_t get_all_safs_files(std::set<std::string> &files);

}

#endif
