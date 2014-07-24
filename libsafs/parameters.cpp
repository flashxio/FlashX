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

#include<algorithm>
#include <sstream>

#include "parameters.h"
#include "common.h"
#include "RAID_config.h"
#include "cache_config.h"

sys_parameters params;

str2int RAID_options[] = {
	{"RAID0", RAID0},
	{"RAID5", RAID5},
	{"HASH", HASH},
};

str2int cache_types[] = {
	{ "tree", TREE_CACHE } ,
	{ "associative", ASSOCIATIVE_CACHE },
	{ "hash-index", HASH_INDEX_CACHE },
	{ "cuckoo", CUCKOO_CACHE },
	{ "lru2q", LRU2Q_CACHE },
	{ "gclock", GCLOCK_CACHE },
};

sys_parameters::sys_parameters()
{
	RAID_block_size = 64;
	SA_min_cell_size = 12;
	test_hit_rate = -1;		// this disables generating the artificial hit rate.
	io_depth_per_file = 32;
	cache_type = ASSOCIATIVE_CACHE;
	cache_size = 512 * 1024 * 1024;
	RAID_mapping_option = RAID5;
	use_virt_aio = false;
	verify_content = false;
	use_flusher = false;
	cache_large_write = false;
	vaio_print_freq = 1000000;
	numa_num_process_threads = 1;
	num_nodes = 1;
	merge_reqs = false;
	max_obj_alloc_size = 100 * 1024 * 1024;
	writable = false;
	max_num_pending_ios = 1000;
	huge_page_enabled = false;
}

void sys_parameters::init(const std::map<std::string, std::string> &configs)
{
	str2int_map cache_map(cache_types, 
			sizeof(cache_types) / sizeof(cache_types[0]));
	str2int_map RAID_option_map(RAID_options,
			sizeof(RAID_options) / sizeof(RAID_options[0]));
	std::map<std::string, std::string>::const_iterator it;

	it = configs.find("RAID_block_size");
	if (it != configs.end()) {
		RAID_block_size = (int) str2size(it->second) / PAGE_SIZE;
		assert(power2(RAID_block_size));
	}

	it = configs.find("SA_cell_size");
	if(it != configs.end()) {
		SA_min_cell_size = atoi(it->second.c_str());
	}

	it = configs.find("io_depth");
	if(it != configs.end()) {
		io_depth_per_file = atoi(it->second.c_str());
	}

	it = configs.find("hit_percent");
	if(it != configs.end()) {
		test_hit_rate = atoi(it->second.c_str());
	}

	it = configs.find("cache_type");
	if(it != configs.end()) {
		cache_type = cache_map.map(it->second);
		if (cache_type < 0) {
			fprintf(stderr, "can't find the right cache type\n");
			exit(1);
		}
	}

	it = configs.find("cache_size");
	if(it != configs.end()) {
		cache_size = str2size(it->second);
	}

	it = configs.find("RAID_mapping");
	if (it != configs.end()) {
		RAID_mapping_option = RAID_option_map.map(it->second);
		if (RAID_mapping_option < 0) {
			fprintf(stderr, "can't find the right mapping option\n");
			exit(1);
		}
	}

	it = configs.find("virt_aio");
	if (it != configs.end())
		use_virt_aio = true;

	it = configs.find("verify_content");
	if (it != configs.end()) {
		verify_content = true;
	}

	it = configs.find("use_flusher");
	if (it != configs.end()) {
		use_flusher = true;
	}

	it = configs.find("cache_large_write");
	if (it != configs.end()) {
		cache_large_write = true;
	}

	it = configs.find("vaio_print_freq");
	if (it != configs.end()) {
		vaio_print_freq = str2size(it->second);
	}

	it = configs.find("numa_num_process_threads");
	if (it != configs.end()) {
		numa_num_process_threads = str2size(it->second);
	}

	it = configs.find("num_nodes");
	if (it != configs.end()) {
		num_nodes = str2size(it->second);
	}

	it = configs.find("merge_reqs");
	if (it != configs.end()) {
		merge_reqs = true;
	}

	it = configs.find("max_obj_alloc_size");
	if (it != configs.end()) {
		max_obj_alloc_size = str2size(it->second);
	}

	it = configs.find("writable");
	if (it != configs.end()) {
		writable = true;
	}

	it = configs.find("max_num_pending_ios");
	if (it != configs.end()) {
		max_num_pending_ios = str2size(it->second);
	}

	it = configs.find("huge_page");
	if (it != configs.end()) {
		huge_page_enabled = true;
	}
}

void sys_parameters::print()
{
	std::cout << "system parameters: " << std::endl;
	std::cout << "\tRAID_block_size: " << RAID_block_size << std::endl;
	std::cout << "\tSA_cell_size: " << SA_min_cell_size << std::endl;
	std::cout << "\tio_depth:" << io_depth_per_file << std::endl;
	std::cout << "\thit_percent: " << test_hit_rate << std::endl;
	std::cout << "\tcache_type: " << cache_type << std::endl;
	std::cout << "\tcache_size: " << cache_size << std::endl;
	std::cout << "\tRAID_mapping: " << RAID_mapping_option << std::endl;
	std::cout << "\tvirt_aio: " << use_virt_aio << std::endl;
	std::cout << "\tverify_content: " << verify_content << std::endl;
	std::cout << "\tuse_flusher: " << use_flusher << std::endl;
	std::cout << "\tcache_large_write: " << cache_large_write << std::endl;
	std::cout << "\tvaio_print_freq: " << vaio_print_freq << std::endl;
	std::cout << "\tnuma_num_process_threads: " << numa_num_process_threads << std::endl;
	std::cout << "\tnum_nodes: " << num_nodes << std::endl;
	std::cout << "\tmerge_reqs: " << merge_reqs << std::endl;
	std::cout << "\tmax_obj_alloc_size: " << max_obj_alloc_size << std::endl;
	std::cout << "\twritable: " << writable << std::endl;
	std::cout << "\tmax_num_pending_ios: " << max_num_pending_ios << std::endl;
	std::cout << "\thuge_page_enabled: " << huge_page_enabled << std::endl;
}

void sys_parameters::print_help()
{
	str2int_map cache_map(cache_types, 
			sizeof(cache_types) / sizeof(cache_types[0]));
	str2int_map RAID_option_map(RAID_options,
			sizeof(RAID_options) / sizeof(RAID_options[0]));

	std::cout << "system parameters: " << std::endl;
	std::cout << "\tRAID_block_size: x(k, K, m, M, g, G)" << std::endl;
	std::cout << "\tSA_cell_size: the size of a page set" << std::endl;
	std::cout << "\tio_depth: the number of pending I/O requests per file"
		<< std::endl;
	std::cout << "\thit_percent: the artificial cache hit rate (%)" << std::endl;
	cache_map.print("\tcache_type: ");
	std::cout << "\tcache_size: x(k, K, m, M, g, G)" << std::endl;
	RAID_option_map.print("\tRAID_mapping: ");
	std::cout << "\tvirt_aio: enable virtual AIO for debugging and performance evaluation"
		<< std::endl;
	std::cout << "\tverify_content: verify data for testing" << std::endl;
	std::cout << "\tuse_flusher: use flusher in the page cache" << std::endl;
	std::cout << "\tcache_large_write: enable large write in the page cache."
		<< std::endl;
	std::cout << "\tvaio_print_freq: how frequently a virtual SSD print stat info (in us)"
		<< std::endl;
	std::cout << "\tnuma_num_process_threads: the number of request processing threads per node in part_global_cached_io"
		<< std::endl;
	std::cout << "\tnum_nodes: the number of NUMA nodes the test program should run"
		<< std::endl;
	std::cout << "\tmerge_reqs: whether or not merge requests in the cached IO"
		<< std::endl;
	std::cout << "\tmax_obj_alloc_size: the maximal size that an object allocator can use."
		<< std::endl;
	std::cout << "\twritable: indicate whether or not to write data" << std::endl;
	std::cout << "\tmax_num_pending_ios: the max number of pending IOs in an I/O instance"
		<< std::endl;
	std::cout << "\thuge_page_enabled: determine whether we use huge page for large chunk of memory"
		<< std::endl;
}

static void read_config_file(const std::string &conf_file,
		std::map<std::string, std::string> &configs)
{
	FILE *f = fopen(conf_file.c_str(), "r");
	if (f == NULL) {
		perror("fopen");
		assert(0);
	}

	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	while ((read = getline(&line, &len, f)) > 0) {
		std::string str = line;
		if (str.length() == 1)
			continue;
		if (line[0] == '#')
			continue;

		size_t found = str.find("=");
		/* if there isn't `=', I assume it's a file name*/
		if (found == std::string::npos) {
			fprintf(stderr, "wrong format: %s\n", line);
			assert(0);
		}

		std::string value = str.substr(found + 1);
		value.erase(std::remove_if(value.begin(), value.end(), isspace),
				value.end());
		std::string key = str.substr(0, found);
		key.erase(std::remove_if(key.begin(), key.end(), isspace),
				key.end());
		bool res = configs.insert(std::pair<std::string, std::string>(key,
					value)).second;
		if (!res)
			configs[key] = value;
	}
	fclose(f);
}

config_map::config_map(const std::string &conf_file)
{
	read_config_file(conf_file, configs);
}

/**
 * All options should have the following format:
 *		key=value.
 * All options that don't have the format are ignored.
 */
void config_map::add_options(const char *argv[], int argc)
{
	for (int i = 0; i < argc; i++) {
		std::string str = argv[i];

		size_t found = str.find("=");
		if (found == std::string::npos) {
			continue;
		}

		std::string value = str.substr(found + 1);
		value.erase(std::remove_if(value.begin(), value.end(), isspace),
				value.end());
		std::string key = str.substr(0, found);
		key.erase(std::remove_if(key.begin(), key.end(), isspace),
				key.end());
		bool res = configs.insert(std::pair<std::string, std::string>(key,
					value)).second;
		if (!res)
			configs[key] = value;
	}
}

size_t split_string(const std::string &str, const std::string &delimiter,
		std::vector<std::string> &splits)
{
	std::stringstream ss(str);
	size_t num = 0;
	while (ss.good()) {
		std::string s;
		ss >> s;
		splits.push_back(s);
		num++;
	}
	return num;
}

void config_map::add_options(const std::string &opts)
{
	std::vector<std::string> parts;
	split_string(opts, " ", parts);
	std::vector<const char *> part_chars;
	for (size_t i = 0; i < parts.size(); i++) {
		part_chars.push_back(parts[i].c_str());
	}
	config_map::add_options(part_chars.data(), part_chars.size());
}
