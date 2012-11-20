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

#include <iostream>
#include <string>
#include <deque>

#define NUM_THREADS 1024

#include "cache.h"
#include "associative_cache.h"
#include "workload.h"
#include "global_cached_private.h"
#include "aio_private.h"
#include "direct_private.h"
#include "part_global_cached_private.h"
#include "read_private.h"
#include "thread_private.h"
#include "remote_access.h"
#include "disk_read_thread.h"

//#define USE_PROCESS

#define GC_SIZE 10000

enum {
	NORMAL,
	DIRECT,
	MMAP,
};

long npages;
int nthreads = 1;
struct timeval global_start;
char static_buf[PAGE_SIZE * 8] __attribute__((aligned(PAGE_SIZE)));
bool verify_read_content = false;
bool high_prio = false;

thread_private *threads[NUM_THREADS];

thread_private *get_thread(int idx)
{
	thread_private *thread = threads[idx];
	assert(idx == thread->get_idx());
	return thread;
}

struct str2int {
	std::string name;
	int value;
};

enum {
	READ_ACCESS,
	DIRECT_ACCESS,
#ifdef ENABLE_AIO
	AIO_ACCESS,
#endif
	GLOBAL_CACHE_ACCESS,
	PART_GLOBAL_ACCESS,
};

str2int access_methods[] = {
	{ "normal", READ_ACCESS },
	{ "direct", DIRECT_ACCESS },
#ifdef ENABLE_AIO
	{ "aio", AIO_ACCESS },
#endif
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
	RAND_PERMUTE,
	LOCAL_RAND_PERMUTE,
	STRIDE,
	BALANCED,
	RAID_RAND,
	USER_FILE_WORKLOAD = -1
};

str2int workloads[] = {
	{ "SEQ", SEQ_OFFSET },
	{ "RAND", RAND_OFFSET },
	{ "RAND_PERMUTE", RAND_PERMUTE },
	{ "LOCAL_RAND_PERMUTE", LOCAL_RAND_PERMUTE },
	{ "STRIDE", STRIDE },
	{ "BALANCED", BALANCED },
	{ "RAID_RAND", RAID_RAND},
	{ "user_file", USER_FILE_WORKLOAD },
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

long str2size(std::string str)
{
	int len = str.length();
	long multiply = 1;
	if (str[len - 1] == 'M' || str[len - 1] == 'm') {
		multiply *= 1024 * 1024;
		str[len - 1] = 0;
	}
	else if (str[len - 1] == 'K' || str[len - 1] == 'k') {
		multiply *= 1024;
		str[len - 1] = 0;
	}
	else if (str[len - 1] == 'G' || str[len - 1] == 'g') {
		multiply *= 1024 * 1024 * 1024;
		str[len - 1] = 0;
	}
	return atol(str.c_str()) * multiply;
}

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
	bool preload = false;
	long cache_size = 512 * 1024 * 1024;
	int entry_size = 128;
	int access_option = -1;
	int ret;
	int i, j;
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;
	int num_files = 0;
	int num_nodes = 1;
	int cache_type = -1;
	std::string file_names[NUM_THREADS];
	int workload = RAND_OFFSET;
	std::string workload_file;
	str2int_map access_map(access_methods,
			sizeof(access_methods) / sizeof(access_methods[0]));
	str2int_map workload_map(workloads, 
			sizeof(workloads) / sizeof(workloads[0]));
	str2int_map cache_map(cache_types, 
			sizeof(cache_types) / sizeof(cache_types[0]));

	if (argc < 5) {
		fprintf(stderr, "there are %d argments\n", argc);
		fprintf(stderr, "read files option pages threads cache_size entry_size preload workload cache_type num_nodes verify_content high_prio\n");
		access_map.print("available access options: ");
		workload_map.print("available workloads: ");
		cache_map.print("available cache types: ");
		exit(1);
	}

	signal(SIGINT, int_handler);

	for (int i = 1; i < argc; i++) {
		std::string str = argv[i];
		size_t found = str.find("=");
		/* if there isn't `=', I assume it's a file name*/
		if (found == std::string::npos) {
			file_names[num_files++] = str;
			continue;
		}

		std::string value = str.substr(found + 1);
		std::string key = str.substr(0, found);
		if (key.compare("option") == 0) {
			access_option = access_map.map(value);
			if (access_option < 0) {
				fprintf(stderr, "can't find the right option\n");
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
		else if(key.compare("preload") == 0) {
			preload = true;
		}
		else if(key.compare("verify_content") == 0) {
			verify_read_content = true;
		}
		else if(key.compare("high_prio") == 0) {
			high_prio = true;
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
	printf("access: %d, npages: %ld, nthreads: %d, cache_size: %ld, cache_type: %d, entry_size: %d, workload: %d, num_nodes: %d, verify_content: %d, high_prio: %d\n",
			access_option, npages, nthreads, cache_size, cache_type, entry_size, workload, num_nodes, verify_read_content, high_prio);

	int num_entries = (int) (((long) npages) * PAGE_SIZE / entry_size);

	if (nthreads > NUM_THREADS) {
		fprintf(stderr, "too many threads\n");
		exit(1);
	}

	int remainings = npages % nthreads;
	int shift = 0;
	long start;
	long end = 0;
	const char *cnames[num_files];
	int num;
	disk_read_thread **read_threads = new disk_read_thread*[num_files];
	for (int k = 0; k < num_files; k++) {
		cnames[k] = file_names[k].c_str();
		if (access_option == GLOBAL_CACHE_ACCESS
				|| access_option == PART_GLOBAL_ACCESS)
			read_threads[k] = new disk_read_thread(cnames[k],
					npages * PAGE_SIZE);
	}

	num = num_files;
	/* initialize the threads' private data. */
	for (j = 0; j < nthreads; j++) {
		switch (access_option) {
			case READ_ACCESS:
				threads[j] = new thread_private(j, entry_size,
						new buffered_io(cnames, num, npages * PAGE_SIZE));
				break;
			case DIRECT_ACCESS:
				threads[j] = new thread_private(j, entry_size,
						new direct_io(cnames, num, npages * PAGE_SIZE));
				break;
#if ENABLE_AIO
			case AIO_ACCESS:
				threads[j] = new thread_private(j, entry_size,
						new async_io(cnames, num, npages * PAGE_SIZE));
				break;
#endif
			case GLOBAL_CACHE_ACCESS:
				{
//					io_interface *underlying = new async_io(cnames, num, npages * PAGE_SIZE);
//					io_interface *underlying = new remote_disk_access(
//							read_threads, num_files);
					io_interface *underlying = new direct_io(cnames, num,
							npages * PAGE_SIZE);
					global_cached_io *io = new global_cached_io(underlying,
							cache_size, cache_type);
					if (preload)
						io->preload(0, npages * PAGE_SIZE);
					threads[j] = new thread_private(j, entry_size, io);
				}
				break;
			case PART_GLOBAL_ACCESS:
				{
					io_interface *underlying = new direct_io(cnames, num,
							npages * PAGE_SIZE);
					threads[j] = new thread_private(j, entry_size,
							new part_global_cached_io(num_nodes, underlying,
								j, cache_size, cache_type));
				}
				break;
			default:
				fprintf(stderr, "wrong access option\n");
				exit(1);
		}
		
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
		printf("thread %d starts %ld ends %ld\n", j, start, end);

		workload_gen *gen;
		workload_chunk *chunk = NULL;
		switch (workload) {
			case SEQ_OFFSET:
				gen = new seq_workload(start, end, entry_size);
				break;
			case RAND_OFFSET:
				gen = new rand_workload(start, end, entry_size);
				break;
			case RAND_PERMUTE:
				gen = new global_rand_permute_workload(num_entries,
						entry_size, start, end);
				break;
			case LOCAL_RAND_PERMUTE:
				gen = new local_rand_permute_workload(start, end, entry_size);
				break;
			case STRIDE:
				gen = new stride_workload(start, end, entry_size);
				break;
			case BALANCED:
				if (chunk == NULL) {
					chunk = new stride_workload_chunk(0,
							(long) npages * PAGE_SIZE / entry_size, entry_size);
				}
				gen = new balanced_workload(chunk);
				break;
			case RAID_RAND:
				/* In this workload each thread starts */
				gen = new RAID0_rand_permute_workload(npages,
						entry_size, nthreads, j);
				break;
			case -1:
				gen = new file_workload(workload_file, nthreads);
				break;
			default:
				fprintf(stderr, "unsupported workload\n");
				exit(1);
		}
		threads[j]->set_workload(gen);
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
	for (i = 0; i < nthreads; i++) {
		ret = threads[i]->start_thread();
		if (ret) {
			perror("pthread_create");
			exit(1);
		}
	}

	for (i = 0; i < nthreads; i++) {
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
	printf("there are %d cells\n", avail_cells);
	printf("there are %d waits for unused\n", num_wait_unused);
	printf("there are %d lock contentions\n", lock_contentions);
#endif
}

const rand_permute *global_rand_permute_workload::permute;
