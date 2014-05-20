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

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include <string.h>
#include <errno.h>
#include <pthread.h>
#include <assert.h>
#include <numa.h>
#include <numaif.h>
#ifdef PROFILER
#include <google/profiler.h>
#endif

#include <vector>
#include <set>
#include <iostream>
#include <string>
#include <deque>
#include <algorithm>

#define NUM_THREADS 1024

#include "workload.h"
#include "thread_private.h"
#include "RAID_config.h"
#include "io_interface.h"
#include "cache_config.h"
#include "config.h"
#include "debugger.h"

//#define USE_PROCESS

#define GC_SIZE 10000

enum {
	NORMAL,
	DIRECT,
	MMAP,
};

struct timeval global_start;
char static_buf[PAGE_SIZE * 8] __attribute__((aligned(PAGE_SIZE)));

thread_private *threads[NUM_THREADS];

str2int access_methods[] = {
	{ "normal", READ_ACCESS },
	{ "direct", DIRECT_ACCESS },
	{ "aio", AIO_ACCESS },
	{ "remote", REMOTE_ACCESS },
	{ "global_cache", GLOBAL_CACHE_ACCESS },
	{ "parted_global", PART_GLOBAL_ACCESS },
};

str2int workloads[] = {
	{ "SEQ", SEQ_OFFSET },
	{ "RAND", RAND_OFFSET },
	{ "RAND_PERMUTE", RAND_PERMUTE },
	{ "HIT_DEFINED", HIT_DEFINED },
	{ "user_file", USER_FILE_WORKLOAD },
};

str2int req_buf_types[] = {
	{ "SINGLE_LARGE", SINGLE_LARGE_BUF },
	{ "SINGLE_SMALL", SINGLE_SMALL_BUF },
	{ "MULTI", MULTI_BUF },
};

#ifdef PROFILER
std::string prof_file = "rand-read.prof";
#endif

test_config config;

void test_config::init(const std::map<std::string, std::string> &configs)
{
	str2int_map access_map(access_methods,
			sizeof(access_methods) / sizeof(access_methods[0]));
	str2int_map workload_map(workloads, 
			sizeof(workloads) / sizeof(workloads[0]));
	str2int_map buf_type_map(req_buf_types,
			sizeof(req_buf_types) / sizeof(req_buf_types[0]));

	std::map<std::string, std::string>::const_iterator it;

	it = configs.find("option");
	if (it != configs.end()) {
		access_option = access_map.map(it->second);
		if (access_option < 0) {
			fprintf(stderr, "can't find the right access option\n");
			exit(1);
		}
	}

	it = configs.find("num_reqs");
	if (it != configs.end()) {
		num_reqs = atoi(it->second.c_str());
	}

	it = configs.find("threads");
	if (it != configs.end()) {
		nthreads = atoi(it->second.c_str());
	}

	it = configs.find("read_percent");
	if (it != configs.end()) {
		read_ratio = (((double) atoi(it->second.c_str())) / 100);
	}

	it = configs.find("repeats");
	if (it != configs.end()) {
		num_repeats = atoi(it->second.c_str());
	}

	it = configs.find("entry_size");
	if (it != configs.end()) {
		entry_size = (int) str2size(it->second);
		workload_gen::set_default_entry_size(entry_size);
	}

	it = configs.find("workload");
	if (it != configs.end()) {
		workload = workload_map.map(it->second);
		if (workload == -1) {
			workload_file = it->second;
		}
	}

	it = configs.find("access");
	if (it != configs.end()) {
		if(it->second.compare("read") == 0)
			workload_gen::set_default_access_method(READ);
		else if(it->second.compare("write") == 0)
			workload_gen::set_default_access_method(WRITE);
		else {
			fprintf(stderr, "wrong default access method\n");
			exit(1);
		}
	}

	it = configs.find("high_prio");
	if (it != configs.end()) {
		high_prio = true;
	}

	it = configs.find("buf_type");
	if (it != configs.end()) {
		buf_type = buf_type_map.map(it->second);
	}

	it = configs.find("sync");
	if (it != configs.end()) {
		use_aio = false;
	}

	it = configs.find("user_compute");
	if (it != configs.end()) {
		user_compute = true;
	}

#ifdef PROFILER
	it = configs.find("prof");
	if (it != configs.end()) {
		prof_file = it->second;
	}
#endif
}

void test_config::print()
{
	printf("the configuration of the test program\n");
	printf("\toption: %d\n", access_option);
	printf("\tnum_reqs: %ld\n", num_reqs);
	printf("\tthreads: %d\n", nthreads);
	printf("\tread_ratio: %f\n", read_ratio);
	printf("\trepeats: %d\n", num_repeats);
	printf("\tentry_size: %d\n", entry_size);
	printf("\tworkload: %d\n", workload);
	printf("\thigh_prio: %d\n", high_prio);
	printf("\tbuf_type: %d\n", buf_type);
	printf("\tsync: %d\n", !use_aio);
	printf("\tuser_compute: %d\n", user_compute);
}

void test_config::print_help()
{
	str2int_map access_map(access_methods,
			sizeof(access_methods) / sizeof(access_methods[0]));
	str2int_map workload_map(workloads, 
			sizeof(workloads) / sizeof(workloads[0]));
	str2int_map buf_type_map(req_buf_types,
			sizeof(req_buf_types) / sizeof(req_buf_types[0]));

	printf("test options:\n");
	access_map.print("access options: ");
	printf("\tnum_reqs: the number of requests generated in a workload\n");
	printf("\tread_ratio: the read percentage of a synthetic workload\n");
	printf("\trepeats: the number of repeats in a random permutation workload\n");
	printf("\tthreads: the number of test threads\n");
	printf("\tentry_size: the size of each access\n");
	workload_map.print("\tworkloads: ");
	printf("\thigh_prio: run the test program in a higher OS priority\n");
	buf_type_map.print("\tbuf types: ");
	printf("\tsync: whether to use sync or async\n");
	printf("\troot_conf: a config file to specify the root paths of the RAID\n");
	printf("\tuser_compute: whether to use user_compute\n");
}

void int_handler(int sig_num)
{
#ifdef PROFILER
	if (!prof_file.empty())
		ProfilerStop();
#endif
#ifdef STATISTICS
	for (int i = 0; i < config.get_nthreads(); i++) {
		if (threads[i])
			threads[i]->print_stat();
	}
	print_io_thread_stat();
#endif
	exit(0);
}

const long TEST_DATA_SIZE = 15L * 1024 * 1024 * 4096;

class debug_workload_gens: public debug_task
{
	std::vector<workload_gen *> workloads;
public:
	debug_workload_gens(const std::vector<workload_gen *> &workloads) {
		this->workloads = workloads;
	}

	void run() {
		for (unsigned i = 0; i < workloads.size(); i++)
			workloads[i]->print_state();
	}
};

int main(int argc, char *argv[])
{
	int ret = 0;
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;

	if (argc < 3) {
		fprintf(stderr, "there are %d argments\n", argc);
		fprintf(stderr, "test_rand_io conf_file data_file [data_files ...] [conf_key=conf_value]\n");

		config.print_help();
		params.print_help();
		exit(1);
	}
	std::string conf_file = argv[1];

	std::vector<std::string> data_files;
	data_files.push_back(argv[2]);
	for (int i = 3; i < argc; i++) {
		// This specifies an configuration option for the test program.
		if (strchr(argv[i], '=') != NULL)
			break;

		data_files.push_back(argv[i]);
	}

	signal(SIGINT, int_handler);
	// The file that contains all data files.

	config_map configs(conf_file);
	configs.add_options((const char **) argv + 3, argc - 3);
	config.init(configs.get_options());
	config.print();

	printf("use a different random sequence\n");
	srandom(time(NULL));

	if (config.get_nthreads() > NUM_THREADS) {
		fprintf(stderr, "too many threads\n");
		exit(1);
	}

	int remainings = config.get_num_reqs() % config.get_nthreads();
	int shift = 0;
	long start;
	long end = 0;

	assert(config.get_nthreads() % params.get_num_nodes() == 0);
	init_io_system(configs);
	std::vector<file_io_factory::shared_ptr> factories;
	for (unsigned i = 0; i < data_files.size(); i++) {
		file_io_factory::shared_ptr factory = create_io_factory(data_files[i],
				config.get_access_option());
		assert(factory);
		factories.push_back(factory);
	}
	std::vector<int> node_id_array;
	for (int i = 0; i < params.get_num_nodes(); i++)
		node_id_array.push_back(i);
	assert(config.get_nthreads() % node_id_array.size() == 0);
	assert(config.get_nthreads() % data_files.size() == 0);
	int nthread_per_node = config.get_nthreads() / node_id_array.size();
	int nthread_per_file = config.get_nthreads() / data_files.size();
	std::vector<workload_gen *> workload_gens;
	for (int i = 0; i < config.get_nthreads(); i++) {
		int node_idx = i / nthread_per_node;
		int file_idx = i / nthread_per_file;
		int node_id = node_id_array[node_idx];
		file_io_factory::shared_ptr factory = factories[file_idx];
		/*
		 * we still assign each thread a range regardless of the number
		 * of threads. read_private will choose the right file descriptor
		 * according to the offset.
		 */
		start = end;
		end = start + ((long) config.get_num_reqs() / config.get_nthreads()
				+ (shift < remainings)) * PAGE_SIZE / config.get_entry_size();
		if (remainings != shift)
			shift++;
#ifdef DEBUG
		printf("thread %d starts %ld ends %ld\n", i, start, end);
#endif

		workload_gen *gen;
		switch (config.get_workload()) {
			case SEQ_OFFSET:
				gen = new seq_workload(start, end, config.get_entry_size());
				break;
			case RAND_OFFSET:
				assert(config.get_read_ratio() >= 0);
				gen = new rand_workload(start, end, config.get_entry_size(),
						end - start, (int) (config.get_read_ratio() * 100));
				break;
			case RAND_PERMUTE:
				assert(config.get_read_ratio() >= 0);
				gen = new global_rand_permute_workload(config.get_entry_size(),
						(((long) config.get_num_reqs()) * PAGE_SIZE) / config.get_entry_size(),
						config.get_num_repeats(), config.get_read_ratio());
				break;
			case -1:
				{
					static long length = 0;
					static workload_t *workloads = NULL;
					if (workloads == NULL)
						workloads = load_file_workload(config.get_workload_file(), length);
					long num_reqs = length;
					if (config.get_num_reqs() >= 0)
						num_reqs = min(config.get_num_reqs(), num_reqs);
					gen = new file_workload(workloads, num_reqs, i, config.get_nthreads(),
							(int) (config.get_read_ratio() * 100));
					break;
				}
			default:
				fprintf(stderr, "unsupported workload\n");
				exit(1);
		}
		workload_gens.push_back(gen);

		threads[i] = new thread_private(node_id, i, config.get_entry_size(), factory, gen);
	}
	debug.register_task(new debug_workload_gens(workload_gens));

	if (config.is_high_prio()) {
		ret = setpriority(PRIO_PROCESS, getpid(), -20);
		if (ret < 0) {
			perror("setpriority");
			exit(1);
		}
	}

	gettimeofday(&start_time, NULL);
	global_start = start_time;
#ifdef PROFILER
	if (!prof_file.empty())
		ProfilerStart(prof_file.c_str());
#endif
	for (int i = 0; i < config.get_nthreads(); i++) {
		threads[i]->start();
	}

	for (int i = 0; i < config.get_nthreads(); i++) {
		threads[i]->join();
		read_bytes += threads[i]->get_read_bytes();
	}
#ifdef PROFILER
	if (!prof_file.empty())
		ProfilerStop();
#endif
	gettimeofday(&end_time, NULL);
	printf("read %ld bytes, takes %f seconds\n",
			read_bytes, end_time.tv_sec - start_time.tv_sec
			+ ((float)(end_time.tv_usec - start_time.tv_usec))/1000000);

#ifdef STATISTICS
	for (int i = 0; i < config.get_nthreads(); i++) {
		threads[i]->print_stat();
	}
#endif
	for (unsigned i = 0; i < workload_gens.size(); i++)
		delete workload_gens[i];
	for (int i = 0; i < config.get_nthreads(); i++)
		delete threads[i];
#ifdef STATISTICS
	print_io_thread_stat();
#endif
	factories.clear();
	destroy_io_system();
}
