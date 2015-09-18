/*
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

#include <stddef.h>
#include <assert.h>

#include "common.h"
#include "native_file.h"

ssize_t native_dir::read_all_files(std::vector<std::string> &files) const
{
	DIR *d = opendir(name.c_str());
	if (d == NULL) {
		fprintf(stderr, "can't open %s: %s\n", name.c_str(),
				strerror(errno));
		return -1;
	}
	struct dirent *entryp;
	long path_limit = pathconf(name.c_str(), _PC_NAME_MAX);
	if (path_limit < 0) {
		perror("pathconf");
		path_limit = sizeof(entryp->d_name);
	}
	int len = offsetof(struct dirent, d_name) + path_limit + 1;
	entryp = (struct dirent *) malloc(len);
	struct dirent *result;
	ssize_t num = 0;
	while (true) {
		int ret = readdir_r(d, entryp, &result);
		if (ret != 0) {
			fprintf(stderr, "can't read dir %s: %s\n", name.c_str(),
					strerror(ret));
			num = -1;
			break;
		}
		// End of the directory stream.
		if (result == NULL)
			break;

		assert(result == entryp);
		std::string file = entryp->d_name;
		if (file != "." && file != "..")
			files.push_back(file);
		num++;
	}
	free(entryp);
	closedir(d);
	return num;
}

bool native_dir::create_dir(bool recursive)
{
	if (recursive) {
		std::vector<std::string> strs;
		split_string(name, '/', strs);
		std::string curr_dir;
		for (unsigned i = 0; i < strs.size(); i++) {
			curr_dir += strs[i] + "/";
			int ret = mkdir(curr_dir.c_str(), 0755);
			if (ret < 0 && errno != EEXIST) {
				perror("mkdir");
				return false;
			}
		}
		return true;
	}
	else {
		int ret = mkdir(name.c_str(), 0755);
		return ret == 0;
	}
}

bool native_dir::delete_dir(bool recursive)
{
	if (recursive) {
		std::vector<std::string> all_files;
		read_all_files(all_files);
		for (unsigned i = 0; i < all_files.size(); i++) {
			native_file child_file(name + "/" + all_files[i]);
			if (!child_file.is_dir()) {
				bool ret = child_file.delete_file();
				if (!ret)
					return false;
			}
			else {
				native_dir child(child_file.get_name());
				bool ret = child.delete_dir(true);
				if (!ret)
					return false;
			}
		}
		int ret = rmdir(name.c_str());
		if (ret < 0)
			fprintf(stderr, "rm %s: %s\n", name.c_str(), strerror(errno));
		return ret == 0;
	}
	else {
		int ret = rmdir(name.c_str());
		if (ret < 0)
			fprintf(stderr, "rm %s: %s\n", name.c_str(), strerror(errno));
		return ret == 0;
	}
}
