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

#define NUM_THREADS 1024

#include "workload.h"
#include "thread_private.h"
#include "RAID_config.h"
#include "io_interface.h"
#include "cache_config.h"

//#define USE_PROCESS

#define GC_SIZE 10000

enum {
	NORMAL,
	DIRECT,
	MMAP,
};

long npages;
int entry_size = 128;
int nthreads = 1;
int buf_type = SINGLE_LARGE_BUF;
int buf_size = PAGE_SIZE;
struct timeval global_start;
char static_buf[PAGE_SIZE * 8] __attribute__((aligned(PAGE_SIZE)));
bool verify_read_content = false;
bool high_prio = false;

thread_private *threads[NUM_THREADS];

struct str2int {
	std::string name;
	int value;
};

str2int access_methods[] = {
	{ "normal", READ_ACCESS },
	{ "direct", DIRECT_ACCESS },
	{ "aio", AIO_ACCESS },
	{ "remote", REMOTE_ACCESS },
	{ "global_cache", GLOBAL_CACHE_ACCESS },
	{ "parted_global", PART_GLOBAL_ACCESS },
};

str2int cache_types[] = {
	{ "tree", TREE_CACHE } ,
	{ "associative", ASSOCIATIVE_CACHE },
	{ "hash-index", HASH_INDEX_CACHE },
	{ "cuckoo", CUCKOO_CACHE },
	{ "lru2q", LRU2Q_CACHE },
	{ "gclock", GCLOCK_CACHE },
};

enum {
	SEQ_OFFSET,
	RAND_OFFSET,
	RAND_SEQ_OFFSET,
	RAND_PERMUTE,
	HIT_DEFINED,
	USER_FILE_WORKLOAD = -1
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

str2int RAID_options[] = {
	{"RAID0", RAID0},
	{"RAID5", RAID5},
	{"HASH", HASH},
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

#ifdef PROFILER
std::string prof_file = "rand-read.prof";
#endif

void int_handler(int sig_num)
{
#ifdef PROFILER
	if (!prof_file.empty())
		ProfilerStop();
#endif
	exit(0);
}

int main(int argc, char *argv[])
{
	long cache_size = 512 * 1024 * 1024;
	int access_option = -1;
	int RAID_mapping_option = RAID5;
	int ret = 0;
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;
	int num_nodes = 1;
	int cache_type = -1;
	int workload = RAND_OFFSET;
	int RAID_block_size = sys_parameters::get_default_RAID_block_size();
	int SA_min_cell_size = sys_parameters::get_default_SA_min_cell_size();
	// No cache hits.
	double hit_ratio = 0;
	// All reads
	double read_ratio = 1;
	int num_repeats = 1;
	std::string workload_file;
	str2int_map access_map(access_methods,
			sizeof(access_methods) / sizeof(access_methods[0]));
	str2int_map workload_map(workloads, 
			sizeof(workloads) / sizeof(workloads[0]));
	str2int_map cache_map(cache_types, 
			sizeof(cache_types) / sizeof(cache_types[0]));
	str2int_map buf_type_map(req_buf_types,
			sizeof(req_buf_types) / sizeof(req_buf_types[0]));
	str2int_map RAID_option_map(RAID_options,
			sizeof(RAID_options) / sizeof(RAID_options[0]));

	if (argc < 5) {
		fprintf(stderr, "there are %d argments\n", argc);
		fprintf(stderr, "read files option pages threads cache_size entry_size workload cache_type num_nodes verify_content high_prio multibuf buf_size hit_percent read_percent repeats RAID_mapping RAID_block_size SA_cell_size\n");
		access_map.print("available access options: ");
		workload_map.print("available workloads: ");
		cache_map.print("available cache types: ");
		buf_type_map.print("available buf types: ");
		exit(1);
	}

	signal(SIGINT, int_handler);
	// The file that contains all data files.
	std::string file_file;

	for (int i = 1; i < argc; i++) {
		std::string str = argv[i];
		size_t found = str.find("=");
		/* if there isn't `=', I assume it's a file name*/
		if (found == std::string::npos) {
			file_file = str;
			continue;
		}

		std::string value = str.substr(found + 1);
		std::string key = str.substr(0, found);
		if (key.compare("option") == 0) {
			access_option = access_map.map(value);
			if (access_option < 0) {
				fprintf(stderr, "can't find the right access option\n");
				exit(1);
			}
		}
		else if(key.compare("cache_type") == 0) {
			cache_type = cache_map.map(value);
		}
		else if(key.compare("pages") == 0) {
			npages = atoi(value.c_str());
		}
		else if(key.compare("threads") == 0) {
			nthreads = atoi(value.c_str());
		}
		else if(key.compare("cache_size") == 0) {
			cache_size = str2size(value);
		}
		else if(key.compare("num_nodes") == 0) {
			num_nodes = str2size(value);
		}
		else if(key.compare("hit_percent") == 0) {
			hit_ratio = (((double) atoi(value.c_str())) / 100);
		}
		else if(key.compare("read_percent") == 0) {
			read_ratio = (((double) atoi(value.c_str())) / 100);
		}
		else if (key.compare("repeats") == 0) {
			num_repeats = atoi(value.c_str());
		}
		else if(key.compare("entry_size") == 0) {
			entry_size = (int) str2size(value);
			workload_gen::set_default_entry_size(entry_size);
		}
		else if(key.compare("workload") == 0) {
			workload = workload_map.map(value);
			if (workload == -1) {
				workload_file = value;
			}
		}
		else if(key.compare("access") == 0) {
			if(value.compare("read") == 0)
				workload_gen::set_default_access_method(READ);
			else if(value.compare("write") == 0)
				workload_gen::set_default_access_method(WRITE);
			else {
				fprintf(stderr, "wrong default access method\n");
				exit(1);
			}
		}
		else if(key.compare("verify_content") == 0) {
			verify_read_content = true;
		}
		else if(key.compare("high_prio") == 0) {
			high_prio = true;
		}
		else if (key.compare("buf_type") == 0) {
			buf_type = buf_type_map.map(value);
		}
		else if (key.compare("buf_size") == 0) {
			buf_size = (int) str2size(value);
		}
		else if (key.compare("RAID_mapping") == 0) {
			RAID_mapping_option = RAID_option_map.map(value);
			if (RAID_mapping_option < 0) {
				fprintf(stderr, "can't find the right mapping option\n");
				exit(1);
			}
		}
		else if(key.compare("RAID_block_size") == 0) {
			RAID_block_size = (int) str2size(value);
		}
		else if(key.compare("SA_cell_size") == 0) {
			SA_min_cell_size = atoi(value.c_str());
		}
#ifdef PROFILER
		else if(key.compare("prof") == 0) {
			prof_file = value;
		}
#endif
		else {
			fprintf(stderr, "wrong option\n");
			exit(1);
		}
	}
	printf("access: %d, npages: %ld, nthreads: %d, cache_size: %ld, cache_type: %d, entry_size: %d, workload: %d, num_nodes: %d, verify_content: %d, high_prio: %d, hit_ratio: %f, read_ratio: %f, repeats: %d, RAID_mapping: %d, RAID block size: %d, SA_cell_size: %d\n",
			access_option, npages, nthreads, cache_size, cache_type, entry_size, workload, num_nodes, verify_read_content, high_prio, hit_ratio, read_ratio, num_repeats, RAID_mapping_option, RAID_block_size, SA_min_cell_size);
	params.init(RAID_block_size, SA_min_cell_size, (int) (hit_ratio * 100));

	int flags = O_RDWR;
	if (access_option != READ_ACCESS) {
		printf("file is opened with direct I/O\n");
		flags |= O_DIRECT;
	}

	if (nthreads > NUM_THREADS) {
		fprintf(stderr, "too many threads\n");
		exit(1);
	}

	int remainings = npages % nthreads;
	int shift = 0;
	long start;
	long end = 0;

	RAID_config raid_conf(file_file, RAID_mapping_option, RAID_block_size);

	std::set<int> node_ids = raid_conf.get_node_ids();
	// In this way, we can guarantee that the cache is created
	// on the nodes with the data files.
	for (int i = 0; i < num_nodes
			&& node_ids.size() < (unsigned) num_nodes; i++)
		node_ids.insert(i);
	std::vector<int> node_id_array;
	// We only get a specified number of nodes.
	for (std::set<int>::const_iterator it = node_ids.begin();
			it != node_ids.end() && (int) node_id_array.size() < num_nodes; it++)
		node_id_array.push_back(*it);

	assert(nthreads % num_nodes == 0);
	assert(node_id_array.size() == (unsigned) num_nodes);
	printf("There are %ld nodes\n", node_id_array.size());

	cache_config *cache_conf = NULL;
	if (access_option == GLOBAL_CACHE_ACCESS)
		cache_conf = new even_cache_config(cache_size, cache_type,
				node_id_array);
	else if (access_option == PART_GLOBAL_ACCESS) {
		assert(num_nodes == 4);
		cache_conf = new even_cache_config(cache_size, cache_type,
				node_id_array);
	}

	init_io_system(raid_conf, node_id_array);

	int io_depth = AIO_DEPTH_PER_FILE / nthreads;
	if (io_depth == 0)
		io_depth = 1;
	file_io_factory *factory = create_io_factory(raid_conf, node_id_array,
			access_option, io_depth, cache_conf);
	int nthread_per_node = nthreads / node_id_array.size();
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
			end = start + ((long) npages / nthreads + (shift < remainings))
				* PAGE_SIZE / entry_size;
			if (remainings != shift)
				shift++;
#ifdef DEBUG
			printf("thread %d starts %ld ends %ld\n", j, start, end);
#endif

			workload_gen *gen;
			switch (workload) {
				case SEQ_OFFSET:
					gen = new seq_workload(start, end, entry_size);
					break;
				case RAND_OFFSET:
					gen = new rand_workload(start, end, entry_size);
					break;
				case RAND_SEQ_OFFSET:
					gen = new rand_seq_workload(start, end, entry_size,
							1024 * 1024 * 4);
					break;
				case RAND_PERMUTE:
					gen = new global_rand_permute_workload(entry_size,
							(((long) npages) * PAGE_SIZE) / entry_size,
							num_repeats, read_ratio);
					break;
				case HIT_DEFINED:
					gen = new cache_hit_defined_workload(entry_size,
							(((long) npages) * PAGE_SIZE) / entry_size,
							cache_size, hit_ratio, read_ratio);
					break;
				case -1:
					{
						static long length = 0;
						static workload_t *workloads = NULL;
						if (workloads == NULL)
							workloads = load_file_workload(workload_file, length);
						gen = new file_workload(workloads, length, j, nthreads);
						break;
					}
				default:
					fprintf(stderr, "unsupported workload\n");
					exit(1);
			}
			workload_gens.push_back(gen);

			int idx = i * nthread_per_node + j;
			threads[idx] = new thread_private(node_id, idx, entry_size, factory, gen);
		}
	}

	if (high_prio) {
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
	for (int i = 0; i < nthreads; i++) {
		ret = threads[i]->start_thread();
		if (ret) {
			perror("pthread_create");
			exit(1);
		}
	}

	for (int i = 0; i < nthreads; i++) {
		threads[i]->wait_thread_end();
		if (ret) {
			perror("pthread_join");
			exit(1);
		}
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
	for (int i = 0; i < nthreads; i++) {
		threads[i]->print_stat();
	}
	print_io_thread_stat();
#endif
	for (unsigned i = 0; i < workload_gens.size(); i++)
		delete workload_gens[i];
	for (int i = 0; i < nthreads; i++)
		delete threads[i];
	destroy_io_factory(factory);
	delete cache_conf;
}
