#ifndef __NATIVE_FILE_H__
#define __NATIVE_FILE_H__

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

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>

#include <string>
#include <vector>

class native_file
{
	std::string file_name;
public:
	native_file(const std::string &file_name) {
		this->file_name = file_name;
	}

	/**
	 * These are file operations on the native file.
	 */
	ssize_t get_size() const {
		struct stat stats;
		if (stat(file_name.c_str(), &stats) < 0) {
			perror("stat");
			return -1;
		}
		return stats.st_size;
	}

	bool exist() const {
		struct stat stats;
		return stat(file_name.c_str(), &stats) == 0;
	}

	bool is_dir() const {
		struct stat stats;
		if (stat(file_name.c_str(), &stats) < 0) {
			perror("stat");
			return false;
		}
		return S_ISDIR(stats.st_mode);
	}

	const std::string &get_name() const {
		return file_name;
	}

	/**
	 * Create/delete a file on the native file system.
	 * If succeed, return 1; otherwise, return 0.
	 */
	bool create_file(size_t size) {
		int fd = open(file_name.c_str(), O_WRONLY | O_CREAT,
				S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
		if (fd < 0) {
			fprintf(stderr, "can't create %s: %s\n", file_name.c_str(),
					strerror(errno));
			return false;
		}
		bool bret = true;
		if (size > 0) {
			int ret = posix_fallocate(fd, 0, size);
			if (ret != 0) {
				fprintf(stderr, "can't allocate %ld bytes for %s, error: %s\n",
						size, file_name.c_str(), strerror(ret));
				bret = false;
			}
		}
		close(fd);
		return bret;
	}

	bool delete_file() {
		int ret = unlink(file_name.c_str());
		if (ret < 0) {
			fprintf(stderr, "can't delete %s\n", file_name.c_str());
			return false;
		}
		return true;
	}
};

class native_dir
{
	std::string name;
public:
	native_dir(const std::string &name) {
		this->name = name;
	}

	bool is_dir() const {
		native_file f(name);
		return f.is_dir();
	}

	bool exist() const {
		native_file f(name);
		return f.exist();
	}

	const std::string &get_name() const {
		return name;
	}

	ssize_t read_all_files(std::vector<std::string> &files) const;

	bool create_dir(bool recursive);
	bool delete_dir(bool recursive);
};

static inline bool file_exist(const std::string &file)
{
	native_file f(file);
	return f.exist();
}

#endif
