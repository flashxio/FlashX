#ifndef __NATIVE_FILE_H__
#define __NATIVE_FILE_H__

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
			return -1;
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
		if (size < 0) {
			fprintf(stderr, "wrong file size\n");
			return false;
		}

		int fd = open(file_name.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR);
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

#endif
