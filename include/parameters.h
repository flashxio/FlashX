#ifndef __MY_PARAMETERS_H__
#define __MY_PARAMETERS_H__

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

#include <assert.h>

#include <iostream>
#include <map>
#include <string>

#define USE_GCLOCK

#define PAGE_SIZE 4096
#define LOG_PAGE_SIZE 12

#define MIN_BLOCK_SIZE 512

/**
 * This maintains key-value pairs for options.
 */
class config_map
{
	std::map<std::string, std::string> configs;
public:
	config_map() {
	}

	config_map(const std::string &conf_file);

	void add_options(const char *opts[], int num);
	// Multiple options may exist in the input string and
	// the options are separated by space.
	void add_options(const std::string &opts);

	const std::string &get_option(const std::string &name) const {
		std::map<std::string, std::string>::const_iterator it
			= configs.find(name);
		assert(it != configs.end());
		return it->second;
	}

	bool has_option(const std::string &name) const {
		return configs.find(name) != configs.end();
	}

	bool read_option(const std::string &name, std::string &value) const {
		if (has_option(name)) {
			value = get_option(name);
			return true;
		}
		else
			return false;
	}

	bool read_option_int(const std::string &name, int &value) const {
		if (has_option(name)) {
			value = atoi(get_option(name).c_str());
			return true;
		}
		else
			return false;
	}

	bool read_option_bool(const std::string &name, bool &value) const {
		if (has_option(name)) {
			value = true;
			return true;
		}
		else
			return false;
	}

	const std::map<std::string, std::string> &get_options() const {
		return configs;
	}
};

class sys_parameters
{
	int RAID_block_size;
	int SA_min_cell_size;
	int test_hit_rate;
	int io_depth_per_file;
	int cache_type;
	long cache_size;
	int RAID_mapping_option;
	bool use_virt_aio;
	bool verify_content;
	bool use_flusher;
	bool cache_large_write;
	int vaio_print_freq;
	int numa_num_process_threads;
	int num_nodes;
	bool merge_reqs;
	int max_obj_alloc_size;
	bool writable;
public:
	sys_parameters();

	void init(const std::map<std::string, std::string> &parameters);

	void print();
	void print_help();

	// in pages
	int get_RAID_block_size() const {
		return RAID_block_size;
	}

	int get_SA_min_cell_size() const {
		return SA_min_cell_size;
	}

	int get_test_hit_rate() const {
		return test_hit_rate;
	}

	int get_aio_depth_per_file() const {
		return io_depth_per_file;
	}

	int get_cache_type() const {
		return cache_type;
	}

	long get_cache_size() const {
		return cache_size;
	}

	int get_RAID_mapping_option() const {
		return RAID_mapping_option;
	}

	bool is_use_virt_aio() const {
		return use_virt_aio;
	}

	bool is_verify_content() const {
		return verify_content;
	}

	bool is_use_flusher() const {
		return use_flusher;
	}

	bool is_cache_large_write() const {
		return cache_large_write;
	}

	int get_vaio_print_freq() const {
		return vaio_print_freq;
	}

	int get_numa_num_process_threads() const {
		return numa_num_process_threads;
	}

	int get_num_nodes() const {
		return num_nodes;
	}

	bool is_merge_reqs() const {
		return merge_reqs;
	}

	int get_max_obj_alloc_size() const {
		return max_obj_alloc_size;
	}

	bool is_writable() const {
		return writable;
	}
};

extern sys_parameters params;

#define AIO_DEPTH_PER_FILE params.get_aio_depth_per_file()

/**
 * The size of an I/O message sent to an I/O thread.
 * It is in the number of I/O requests.
 */
const int IO_MSG_SIZE = AIO_DEPTH_PER_FILE;
/**
 * The size of a high-priority I/O queue.
 * It's in the number of I/O messages.
 */
const int IO_QUEUE_SIZE = 10;
const int MAX_FETCH_REQS = 3;
const int AIO_COMPLETE_BUF_SIZE = 8;

/**
 * This number defines the max number of flushes sent to an SSD.
 * The number is set based on the performance result.
 * It varies in different SSDs. It's better that we can estimate at runtime.
 */
const int MAX_NUM_FLUSHES_PER_FILE = 2048;

const int MAX_DISK_CACHED_REQS = 1000;

/**
 * The number of requests issued by user applications
 * in one access.
 */
const int NUM_REQS_BY_USER = 1000;

/**
 * The min size of IO vector allocated for an IO request..
 */
const int MIN_NUM_ALLOC_IOVECS = 16;
const int NUM_EMBEDDED_IOVECS = 1;

const int CELL_SIZE = 16;
const int CELL_MIN_NUM_PAGES = 8;

const int MAX_NUM_DIRTY_CELLS_IN_QUEUE = 1000;
const int DIRTY_PAGES_THRESHOLD = 1;
const int NUM_WRITEBACK_DIRTY_PAGES = 2;
const int MAX_NUM_WRITEBACK = 4;
const int DISCARD_FLUSH_THRESHOLD = 6;

const long MAX_CACHE_SIZE = 512L * 1024 * 1024 * 1024;

const int NUMA_MSG_SIZE = 4096;
const int NUMA_REPLY_CACHE_SIZE = 50;

const int LOCAL_BUF_SIZE = 100;

struct str2int {
	std::string name;
	int value;
};

class str2int_map {
	str2int *maps;
	int num;
public:
	str2int_map(str2int *maps, int num) {
		this->maps = maps;
		this->num = num;
	}

	int map(const std::string &str) {
		for (int i = 0; i < num; i++) {
			if (maps[i].name.compare(str) == 0)
				return i;
		}
		return -1;
	}

	void print(const std::string &str) {
		std::cout<<str;
		for (int i = 0; i < num; i++) {
			std::cout<<maps[i].name;
			if (i < num - 1)
				std::cout<<", ";
		}
		std::cout<<std::endl;
	}
};

#endif
