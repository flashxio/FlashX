#ifndef __SAFS_UTIL_H__
#define __SAFS_UTIL_H__

/**
 * Copyright 2019
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
#include <fcntl.h>

#include <string>

#include "io_interface.h"
#include "native_file.h"
#include "safs_file.h"
#include "file_mapper.h"
#include "RAID_config.h"

using namespace safs;

namespace fg { namespace utils {

static config_map::ptr configs;
static const int BUF_SIZE = 1024 * 64 * 4096;

ssize_t complete_read(int fd, char *buf, size_t count);

class data_source
{
public:
	virtual ssize_t get_data(off_t off, size_t size, char *buf) const = 0;
	virtual size_t get_size() const = 0;
};

class file_data_source: public data_source
{
	int fd;
	size_t file_size;
public:
	file_data_source(const std::string &ext_file);
	virtual ssize_t get_data(off_t off, size_t size, char *buf) const override;
	virtual size_t get_size() const override { return file_size; }
};

class synthetic_data_source: public data_source
{
	size_t size;
public:
	synthetic_data_source(size_t size) : size(size) { }

	virtual ssize_t get_data(off_t off, size_t size, char *buf) const override;
	virtual size_t get_size() const override {
		return size;
	}
};

class verify_callback: public callback
{
	char *orig_buf;
	data_source *source;
	size_t verified_bytes;
	file_mapper::const_ptr fmapper;
public:
	verify_callback(data_source *source, file_mapper::const_ptr fmapper);

	~verify_callback() {
		free(orig_buf);
	}

	int invoke(io_request *rqs[], int num);
	size_t get_verified_bytes() const {
		return verified_bytes;
	}
};

void comm_verify_file(int argc, char *argv[]);
void comm_load_file2fs(int argc, char *argv[]);
void comm_load_part_file2fs(int argc, char *argv[]);
void comm_create_file(int argc, char *argv[]);
void print_help();
void comm_help(int argc, char *argv[]);
void comm_list(int argc, char *argv[]);
void comm_delete_file(int argc, char *argv[]);
void comm_export(int argc, char *argv[]);
void comm_show_info(int argc, char *argv[]);
void comm_rename(int argc, char *argv[]);
typedef void (*command_func_t)(int argc, char *argv[]);

struct command
{
	std::string name;
	command_func_t func;
	std::string help_info;
};

int get_num_commands();
const command *get_command(const std::string &name);
void print_help();

} } // End namespace fg::utils

#endif
