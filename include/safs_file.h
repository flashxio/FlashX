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

#include "common.h"
#include "RAID_config.h"

class safs_file
{
	// The collection of native files.
	std::vector<file_info> native_dirs;
	std::string name;
public:
	safs_file(const RAID_config &conf, const std::string &file_name) {
		conf.get_disks(native_dirs);
		for (unsigned i = 0; i < native_dirs.size(); i++)
			native_dirs[i].name += "/" + file_name;
		this->name = file_name;
	}

	const std::string &get_name() const {
		return name;
	}

	bool exist() const;
	size_t get_file_size() const;
	bool create_file(size_t file_size);
	bool delete_file();
};

#endif
