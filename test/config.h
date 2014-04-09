#ifndef __TEST_CONFIG_H__
#define __TEST_CONFIG_H__

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
		entry_size = PAGE_SIZE;
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
