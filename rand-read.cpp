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
int entry_size = 128;
int nthreads = 1;
int buf_type = SINGLE_LARGE_BUF;
int buf_size = PAGE_SIZE;
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
	REMOTE_ACCESS,
	GLOBAL_CACHE_ACCESS,
	PART_GLOBAL_ACCESS,
};

str2int access_methods[] = {
	{ "normal", READ_ACCESS },
	{ "direct", DIRECT_ACCESS },
#ifdef ENABLE_AIO
	{ "aio", AIO_ACCESS },
#endif
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
	RAND_PERMUTE,
	HIT_DEFINED,
	USER_FILE_WORKLOAD = -1
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
	int access_option = -1;
	int ret;
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;
	int num_nodes = 1;
	int cache_type = -1;
	int workload = RAND_OFFSET;
	std::string workload_file;
	str2int_map access_map(access_methods,
			sizeof(access_methods) / sizeof(access_methods[0]));
	str2int_map workload_map(workloads, 
			sizeof(workloads) / sizeof(workloads[0]));
	str2int_map cache_map(cache_types, 
			sizeof(cache_types) / sizeof(cache_types[0]));
	str2int_map buf_type_map(req_buf_types,
			sizeof(req_buf_types) / sizeof(req_buf_types[0]));

	if (argc < 5) {
		fprintf(stderr, "there are %d argments\n", argc);
		fprintf(stderr, "read files option pages threads cache_size entry_size preload workload cache_type num_nodes verify_content high_prio multibuf buf_size\n");
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
		else if (key.compare("buf_type") == 0) {
			buf_type = buf_type_map.map(value);
		}
		else if (key.compare("buf_size") == 0) {
			buf_size = (int) str2size(value);
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

	std::vector<file_info> files;
	int num_files = retrieve_data_files(file_file, files);
	NUMA_access_mapper *NUMA_mapper = new NUMA_access_mapper(files);

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
	std::set<int> node_ids;
	for (int k = 0; k < num_files; k++) {
		cnames[k] = files[k].name.c_str();
		if (access_option == GLOBAL_CACHE_ACCESS
				|| access_option == PART_GLOBAL_ACCESS
				|| access_option == REMOTE_ACCESS)
			read_threads[k] = new disk_read_thread(cnames[k],
					npages * PAGE_SIZE / num_files, files[k].node_id);
		node_ids.insert(files[k].node_id);
	}

	// In this way, we can guarantee that the cache is created
	// on the nodes with the data files.
	for (int i = 0; i < num_nodes
			&& node_ids.size() < (unsigned) num_nodes; i++)
		node_ids.insert(i);
	std::vector<int> node_id_array(node_ids.begin(), node_ids.end());

	num = num_files;
	assert(nthreads % num_nodes == 0);
	assert(node_id_array.size() >= (unsigned) num_nodes);
	int nthreads_per_node = nthreads / num_nodes;
	for (int i = 0, j = 0; i < num_nodes; i++) {
		int node_id = node_id_array[i];
		numa_run_on_node(node_id);
		/* initialize the threads' private data. */
		for (int k = 0; k < nthreads_per_node; k++, j++) {
			switch (access_option) {
				case READ_ACCESS:
					threads[j] = new thread_private(j, entry_size,
							new buffered_io(cnames, num, npages * PAGE_SIZE, node_id));
					break;
				case DIRECT_ACCESS:
					threads[j] = new thread_private(j, entry_size,
							new direct_io(cnames, num, npages * PAGE_SIZE, node_id));
					break;
#if ENABLE_AIO
				case AIO_ACCESS:
					{
						int depth_per_file = AIO_DEPTH_PER_FILE / nthreads;
						if (depth_per_file == 0)
							depth_per_file = 1;
						threads[j] = new thread_private(j, entry_size,
								new async_io(cnames, num, npages * PAGE_SIZE,
									depth_per_file, node_id));
					}
					break;
#endif
				case REMOTE_ACCESS:
					threads[j] = new thread_private(j, entry_size,
							new remote_disk_access(read_threads, num_files, node_id));
					break;
				case GLOBAL_CACHE_ACCESS:
					{
						io_interface *underlying = new remote_disk_access(
								read_threads, num_files, node_id);
						global_cached_io *io = new global_cached_io(underlying,
								cache_size, cache_type);
						if (preload)
							io->preload(0, npages * PAGE_SIZE);
						threads[j] = new thread_private(j, entry_size, io);
					}
					break;
				case PART_GLOBAL_ACCESS:
					{
						io_interface *underlying = new remote_disk_access(
								read_threads, num_files, node_id);
						threads[j] = new thread_private(j, entry_size,
								new part_global_cached_io(num_nodes, underlying,
									j, cache_size, cache_type, NUMA_mapper));
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
			switch (workload) {
				case SEQ_OFFSET:
					gen = new seq_workload(start, end, entry_size);
					break;
				case RAND_OFFSET:
					gen = new rand_workload(start, end, entry_size);
					break;
				case RAND_PERMUTE:
					gen = new global_rand_permute_workload(entry_size,
							start, end);
					break;
				case HIT_DEFINED:
					gen = new cache_hit_defined_workload(entry_size, start,
							end, cache_size, 0);
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
			threads[j]->set_workload(gen);
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
	for (int i = 0; i < num_files; i++) {
		disk_read_thread *t = read_threads[i];
		if (t)
			printf("queue on file %s wait for requests for %d times, is full for %d times, and %d accesses and %d io waits\n",
					cnames[i], read_threads[i]->get_queue()->get_num_empty(),
					read_threads[i]->get_queue()->get_num_empty(), read_threads[i]->get_num_accesses(),
					read_threads[i]->get_num_iowait());
	}
	for (int i = 0; i < nthreads; i++) {
		threads[i]->print_stat();
	}
	printf("there are %d cells\n", avail_cells);
	printf("there are %d waits for unused\n", num_wait_unused);
	printf("there are %d lock contentions\n", lock_contentions);
#endif
}
