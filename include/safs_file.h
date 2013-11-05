#ifndef __SAFS_FILE_H__
#define __SAFS_FILE_H__

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
