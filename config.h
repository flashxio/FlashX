#ifndef __TEST_CONFIG_H__
#define __TEST_CONFIG_H__

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
	int buf_size;
	bool high_prio;
	bool use_aio;
	int workload;
	// All reads
	double read_ratio;
	int num_repeats;
	std::string workload_file;
public:
	test_config() {
		access_option = -1;
		num_reqs = -1;
		entry_size = PAGE_SIZE;
		nthreads = 1;
		buf_type = SINGLE_LARGE_BUF;
		buf_size = PAGE_SIZE;
		high_prio = false;
		use_aio = true;
		workload = RAND_OFFSET;
		read_ratio = -1;
		num_repeats = 1;
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

	int get_buf_size() const {
		return buf_size;
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
};

extern test_config config;

#endif
