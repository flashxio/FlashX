#ifndef __TEST_CONFIG_H__
#define __TEST_CONFIG_H__

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

#include <map>
#include <string>

#include "RAID_config.h"
#include "cache_config.h"

enum {
	SINGLE_LARGE_BUF,
	SINGLE_SMALL_BUF,
	MULTI_BUF,
};

enum {
	SEQ_OFFSET,
	RAND_OFFSET,
	RAND_SEQ_OFFSET,
	RAND_PERMUTE,
	HIT_DEFINED,
	USER_FILE_WORKLOAD = -1
};

class test_config
{
	int access_option;
	long num_reqs;
	int entry_size;
	int nthreads;
	int buf_type;
	bool high_prio;
	bool use_aio;
	int workload;
	// All reads
	double read_ratio;
	int num_repeats;
	std::string workload_file;
	bool user_compute;
public:
	test_config() {
		access_option = -1;
		num_reqs = -1;
		entry_size = safs::PAGE_SIZE;
		nthreads = 1;
		buf_type = SINGLE_LARGE_BUF;
		high_prio = false;
		use_aio = true;
		workload = RAND_OFFSET;
		read_ratio = -1;
		num_repeats = 1;
		user_compute = false;
	}

	void init(const std::map<std::string, std::string> &configs);
	void print();
	void print_help();

	int get_access_option() const {
		return access_option;
	}

	long get_num_reqs() const {
		return num_reqs;
	}

	int get_entry_size() const {
		return entry_size;
	}

	int get_nthreads() const {
		return nthreads;
	}

	int get_buf_type() const {
		return buf_type;
	}

	bool is_high_prio() const {
		return high_prio;
	}

	bool is_use_aio() const {
		return use_aio;
	}

	int get_workload() const {
		return workload;
	}

	double get_read_ratio() const {
		return read_ratio;
	}

	int get_num_repeats() const {
		return num_repeats;
	}

	std::string get_workload_file() const {
		return workload_file;
	}

	bool is_user_compute() const {
		return user_compute;
	}
};

extern test_config config;

#endif
