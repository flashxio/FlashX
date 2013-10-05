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
	{ "RAND_SEQ", RAND_SEQ_OFFSET },
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

void test_config::init(const std::map<std::string, std::string> configs)
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

	it = configs.find("pages");
	if (it != configs.end()) {
		npages = atoi(it->second.c_str());
	}

	it = configs.find("threads");
	if (it != configs.end()) {
		nthreads = atoi(it->second.c_str());
	}

	it = configs.find("num_nodes");
	if (it != configs.end()) {
		num_nodes = str2size(it->second);
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

	it = configs.find("verify_content");
	if (it != configs.end()) {
		verify_read_content = true;
	}

	it = configs.find("high_prio");
	if (it != configs.end()) {
		high_prio = true;
	}

	it = configs.find("buf_type");
	if (it != configs.end()) {
		buf_type = buf_type_map.map(it->second);
	}

	it = configs.find("buf_size");
	if (it != configs.end()) {
		buf_size = (int) str2size(it->second);
	}

	it = configs.find("sync");
	if (it != configs.end()) {
		use_aio = false;
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
	printf("\tpages: %ld\n", npages);
	printf("\tthreads: %d\n", nthreads);
	printf("\tnum_nodes: %d\n", num_nodes);
	printf("\tread_ratio: %f\n", read_ratio);
	printf("\trepeats: %d\n", num_repeats);
	printf("\tentry_size: %d\n", entry_size);
	printf("\tworkload: %d\n", workload);
	printf("\tverify_content: %d\n", verify_read_content);
	printf("\thigh_prio: %d\n", high_prio);
	printf("\tbuf_type: %d\n", buf_type);
	printf("\tbuf_size: %d\n", buf_size);
	printf("\tsync: %d\n", !use_aio);
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
	printf("\tpages: the number of pages to access in a synthetic workload\n");
	printf("\tread_ratio: the read percentage of a synthetic workload\n");
	printf("\trepeats: the number of repeats in a random permutation workload\n");
	printf("\tthreads: the number of test threads\n");
	printf("\tnum_nodes: the number of NUMA nodes the test program should run\n");
	printf("\tentry_size: the size of each access\n");
	workload_map.print("workloads: ");
	printf("\tverify_content: to verify the data read from the disks\n");
	printf("\thigh_prio: run the test program in a higher OS priority\n");
	buf_type_map.print("buf types: ");
	printf("\tbuf_size: the buffer size for each access\n");
	printf("\tsync: whether to use sync or async\n");
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

void read_config_file(const std::string &conf_file,
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
		configs.insert(std::pair<std::string, std::string>(key, value));
	}
	fclose(f);
}

void parse_args(int argc, char *argv[],
		std::map<std::string, std::string> &configs)
{
	for (int i = 0; i < argc; i++) {
		std::string str = argv[i];

		size_t found = str.find("=");
		/* if there isn't `=', I assume it's a file name*/
		if (found == std::string::npos) {
			fprintf(stderr, "wrong format: %s\n", argv[i]);
			assert(0);
		}

		std::string value = str.substr(found + 1);
		value.erase(std::remove_if(value.begin(), value.end(), isspace),
				value.end());
		std::string key = str.substr(0, found);
		key.erase(std::remove_if(key.begin(), key.end(), isspace),
				key.end());
		configs.insert(std::pair<std::string, std::string>(key, value));
	}
}

int main(int argc, char *argv[])
{
	int ret = 0;
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;

	if (argc < 3) {
		fprintf(stderr, "there are %d argments\n", argc);
		fprintf(stderr, "read conf_file data_file [conf_key=conf_value]\n");

		config.print_help();
		sys_params.print_help();
		exit(1);
	}
	std::string conf_file = argv[1];
	std::string file_file = argv[2];

	signal(SIGINT, int_handler);
	// The file that contains all data files.

	std::map<std::string, std::string> configs;
	read_config_file(conf_file, configs);
	parse_args(argc - 3, argv + 3, configs);
	sys_params.init(configs);
	config.init(configs);

	sys_params.print();
	config.print();

	printf("use a different random sequence\n");
	srandom(time(NULL));

	int flags = O_RDWR;
	if (config.get_access_option() != READ_ACCESS) {
		printf("file is opened with direct I/O\n");
		flags |= O_DIRECT;
	}

	if (config.get_nthreads() > NUM_THREADS) {
		fprintf(stderr, "too many threads\n");
		exit(1);
	}

	int remainings = config.get_npages() % config.get_nthreads();
	int shift = 0;
	long start;
	long end = 0;

	RAID_config raid_conf(file_file, sys_params.get_RAID_mapping_option(),
			sys_params.get_RAID_block_size());

	std::set<int> node_ids = raid_conf.get_node_ids();
	// In this way, we can guarantee that the cache is created
	// on the nodes with the data files.
	for (int i = 0; i < config.get_num_nodes()
			&& node_ids.size() < (unsigned) config.get_num_nodes(); i++)
		node_ids.insert(i);
	std::vector<int> node_id_array;
	// We only get a specified number of nodes.
	for (std::set<int>::const_iterator it = node_ids.begin();
			it != node_ids.end()
			&& (int) node_id_array.size() < config.get_num_nodes(); it++)
		node_id_array.push_back(*it);

	assert(config.get_nthreads() % config.get_num_nodes() == 0);
	assert(node_id_array.size() == (unsigned) config.get_num_nodes());
	printf("There are %ld nodes\n", node_id_array.size());

	cache_config *cache_conf = NULL;
	if (config.get_access_option() == GLOBAL_CACHE_ACCESS)
		cache_conf = new even_cache_config(sys_params.get_cache_size(),
				sys_params.get_cache_type(), node_id_array);
	else if (config.get_access_option() == PART_GLOBAL_ACCESS) {
		assert(config.get_num_nodes() == 4);
		cache_conf = new even_cache_config(sys_params.get_cache_size(),
				sys_params.get_cache_type(), node_id_array);
	}

	init_io_system(raid_conf, node_id_array);

	file_io_factory *factory = create_io_factory(raid_conf, node_id_array,
			config.get_access_option(), sys_params.get_aio_depth_per_file(),
			cache_conf);
	int nthread_per_node = config.get_nthreads() / node_id_array.size();
	std::vector<workload_gen *> workload_gens;
	for (unsigned i = 0; i < node_id_array.size(); i++) {
		int node_id = node_id_array[i];
		for (int j = 0; j < nthread_per_node; j++) {
			/*
			 * we still assign each thread a range regardless of the number
			 * of threads. read_private will choose the right file descriptor
			 * according to the offset.
			 */
			start = end;
			end = start + ((long) config.get_npages() / config.get_nthreads()
					+ (shift < remainings)) * PAGE_SIZE / config.get_entry_size();
			if (remainings != shift)
				shift++;
#ifdef DEBUG
			printf("thread %d starts %ld ends %ld\n", j, start, end);
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
				case RAND_SEQ_OFFSET:
					gen = new rand_seq_workload(start, end, config.get_entry_size(),
							1024 * 1024 * 4);
					break;
				case RAND_PERMUTE:
					assert(config.get_read_ratio() >= 0);
					gen = new global_rand_permute_workload(config.get_entry_size(),
							(((long) config.get_npages()) * PAGE_SIZE) / config.get_entry_size(),
							config.get_num_repeats(), config.get_read_ratio());
					break;
				case -1:
					{
						static long length = 0;
						static workload_t *workloads = NULL;
						if (workloads == NULL)
							workloads = load_file_workload(config.get_workload_file(), length);
						gen = new file_workload(workloads, length, j, config.get_nthreads(),
								(int) (config.get_read_ratio() * 100));
						break;
					}
				default:
					fprintf(stderr, "unsupported workload\n");
					exit(1);
			}
			workload_gens.push_back(gen);

			int idx = i * nthread_per_node + j;
			threads[idx] = new thread_private(node_id, idx, config.get_entry_size(), factory, gen);
		}
	}

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
	print_io_thread_stat();
#endif
	for (unsigned i = 0; i < workload_gens.size(); i++)
		delete workload_gens[i];
	for (int i = 0; i < config.get_nthreads(); i++)
		delete threads[i];
	destroy_io_factory(factory);
	delete cache_conf;
}
