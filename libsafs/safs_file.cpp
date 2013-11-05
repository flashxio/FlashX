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

#include <limits.h>

#include "native_file.h"
#include "safs_file.h"

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
	size_t min_size = LONG_MAX;
	if (!exist())
		return -1;
	for (unsigned i = 0; i < native_dirs.size(); i++) {
		native_dir dir(native_dirs[i].name);
		std::vector<std::string> local_files;
		dir.read_all_files(local_files);
		native_file f(dir.get_name() + "/" + local_files[0]);
		min_size = min<size_t>(min_size, f.get_size());
	}
	return min_size * native_dirs.size();
}

bool safs_file::create_file(size_t file_size)
{
	size_t size_per_disk = file_size / native_dirs.size();
	if (file_size % native_dirs.size() > 0)
		size_per_disk++;

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
