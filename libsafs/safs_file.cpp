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

#include <limits.h>

#include "log.h"
#include "native_file.h"
#include "safs_file.h"
#include "RAID_config.h"
#include "io_interface.h"

namespace safs
{

safs_file::safs_file(const RAID_config &conf, const std::string &file_name)
{
	conf.get_disks(native_dirs);
	for (unsigned i = 0; i < native_dirs.size(); i++)
		native_dirs[i].name += "/" + file_name;
	this->name = file_name;
}

bool safs_file::exist() const
{
	std::set<int> part_ids;
	for (unsigned i = 0; i < native_dirs.size(); i++) {
		native_dir dir(native_dirs[i].name);
		if (!dir.exist())
			return false;
		std::vector<std::string> files;
		dir.read_all_files(files);
		if (files.size() != 1) {
			fprintf(stderr, "%s doesn't have exactly one file\n",
					dir.get_name().c_str());
			return false;
		}
		part_ids.insert(atoi(files[0].c_str()));
	}
	if (part_ids.size() < native_dirs.size()) {
		fprintf(stderr, "there are duplicated partition ids in %s.\n",
				name.c_str());
		return false;
	}
	return true;
}

size_t safs_file::get_file_size() const
{
	if (!exist())
		return -1;
	size_t ret = 0;
	for (unsigned i = 0; i < native_dirs.size(); i++) {
		native_dir dir(native_dirs[i].name);
		std::vector<std::string> local_files;
		dir.read_all_files(local_files);
		assert(local_files.size() == 1);
		native_file f(dir.get_name() + "/" + local_files[0]);
		ret += f.get_size();
	}
	return ret;
}

bool safs_file::create_file(size_t file_size)
{
	size_t size_per_disk = file_size / native_dirs.size();
	if (file_size % native_dirs.size() > 0)
		size_per_disk++;
	size_per_disk = ROUNDUP(size_per_disk, 512);

	for (unsigned i = 0; i < native_dirs.size(); i++) {
		native_dir dir(native_dirs[i].name);
		bool ret = dir.create_dir(true);
		if (!ret)
			return false;
		native_file f(dir.get_name() + "/" + itoa(i));
		ret = f.create_file(size_per_disk);
		if (!ret)
			return false;
	}
	return true;
}

bool safs_file::delete_file()
{
	for (unsigned i = 0; i < native_dirs.size(); i++) {
		native_dir dir(native_dirs[i].name);
		bool ret = dir.delete_dir(true);
		if (!ret)
			return false;
	}
	return true;
}

size_t get_all_safs_files(std::set<std::string> &files)
{
	std::set<std::string> all_files;
	const RAID_config &conf = get_sys_RAID_conf();

	// First find all individual file names in the root directories.
	for (int i = 0; i < conf.get_num_disks(); i++) {
		std::string dir_name = conf.get_disk(i).name;
		native_dir dir(dir_name);
		std::vector<std::string> file_names;
		dir.read_all_files(file_names);
		all_files.insert(file_names.begin(), file_names.end());
	}

	for (std::set<std::string>::const_iterator it = all_files.begin();
			it != all_files.end(); it++) {
		safs_file file(conf, *it);
		if (file.exist()) {
			files.insert(*it);
		}
		else {
			BOOST_LOG_TRIVIAL(error) << boost::format("%1% is corrupted")
				% file.get_name();
		}
	}
	return 0;
}

}
