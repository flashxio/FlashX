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
#include <sys/syscall.h>
#define gettid() syscall(__NR_gettid)

#include <iostream>
#include <string>
#include <deque>

#define NUM_THREADS 1024

#include "cache.h"
#include "tree_cache.h"
#include "associative_cache.h"
#include "cuckoo_cache.h"
#include "workload.h"
#include "LRU2Q.h"
#include "global_cached_private.h"
#include "aio_private.h"
#include "direct_private.h"
#include "mmap_private.h"
#include "part_cached_private.h"
#include "part_global_cached_private.h"
#include "read_private.h"
#include "thread_private.h"

//#define USE_PROCESS

#define GC_SIZE 10000
#define BULK_SIZE 1000

enum {
	NORMAL,
	DIRECT,
	MMAP,
};

void *page::data_start;

long npages;
int nthreads = 1;
struct timeval global_start;
char static_buf[PAGE_SIZE * 8] __attribute__((aligned(PAGE_SIZE)));
int access_method = READ;
bool verify_read_content = false;
bool high_prio = false;

thread_private *threads[NUM_THREADS];

void *rand_read(void *arg)
{
	ssize_t ret = -1;
	thread_private *priv = threads[(long) arg];
	int entry_size = priv->get_entry_size();

#if NCPUS > 0
	cpu_set_t cpuset;
	pthread_t thread = pthread_self();
	CPU_ZERO(&cpuset);
	int cpu_num = priv->idx % NCPUS;
	CPU_SET(cpu_num, &cpuset);
	ret = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
	if (ret != 0) {
		perror("pthread_setaffinity_np");
		exit(1);
	}
	printf("attach thread %d to CPU %d\n", priv->idx, cpu_num);
#endif

#if NUM_NODES > 1
	struct bitmask *nodemask = numa_allocate_cpumask();
	numa_bitmask_clearall(nodemask);
	int node_num = priv->idx / (nthreads / NUM_NODES);
	printf("thread %d is associated to node %d\n", priv->idx, node_num);
	numa_bitmask_setbit(nodemask, node_num);
	numa_bind(nodemask);
	numa_set_strict(1);
	numa_set_bind_policy(1);
#endif

	printf("pid: %d, tid: %ld\n", getpid(), gettid());
	priv->thread_init();
	rand_buf *buf = priv->buf;

	gettimeofday(&priv->start_time, NULL);
	while (priv->gen->has_next()) {
		if (priv->support_bulk()) {
			io_request reqs[BULK_SIZE];
			int i;
//			io_request *reqs = gc->allocate_obj(BULK_SIZE);
			for (i = 0; i < BULK_SIZE
					&& priv->gen->has_next() && !buf->is_full(); i++) {
				reqs[i].init(buf->next_entry(),
						priv->gen->next_offset(), entry_size, READ, priv);
			}
			// TODO right now it only support read.
			ret = priv->access(reqs, i, READ);
			if (ret < 0) {
				perror("access_vector");
				exit(1);
			}
		}
		else {
			char *entry = buf->next_entry();
			off_t off = priv->gen->next_offset();

			/*
			 * generate the data for writing the file,
			 * so the data in the file isn't changed.
			 */
			if (access_method == WRITE) {
				unsigned long *p = (unsigned long *) entry;
				long start = off / sizeof(long);
				for (unsigned int i = 0; i < entry_size / sizeof(*p); i++)
					p[i] = start++;
			}

			ret = priv->access(entry, off, entry_size, access_method);
			if (ret > 0) {
//				assert(ret == buf->get_entry_size());
				if (access_method == READ && verify_read_content) {
					if (*(unsigned long *) entry != off / sizeof(long))
						printf("entry: %ld, off: %ld\n", *(unsigned long *) entry, off / sizeof(long));
					assert(*(unsigned long *) entry == off / sizeof(long));
				}
				if (ret > 0)
					priv->read_bytes += ret;
				else
					break;
			}
			buf->free_entry(entry);
			if (ret < 0) {
				perror("access");
				exit(1);
			}
		}
	}
	priv->cleanup();
	printf("thread %d exits\n", priv->idx);
	gettimeofday(&priv->end_time, NULL);
	
#ifdef USE_PROCESS
	exit(priv->read_bytes);
#else
	pthread_exit((void *) priv->read_bytes);
#endif
}

int process_create(pid_t *pid, void (*func)(void *), void *priv)
{
	pid_t id = fork();

	if (id < 0)
		return -1;

	if (id == 0) {	// child
		func(priv);
		exit(0);
	}

	if (id > 0)
		*pid = id;
	return 0;
}

int process_join(pid_t pid)
{
	int status;
	pid_t ret = waitpid(pid, &status, 0);
	return ret < 0 ? ret : 0;
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
	MMAP_ACCESS,
	LOCAL_CACHE_ACCESS,
	GLOBAL_CACHE_ACCESS,
	PART_GLOBAL_ACCESS,
};

str2int access_methods[] = {
	{ "normal", READ_ACCESS },
	{ "direct", DIRECT_ACCESS },
#ifdef ENABLE_AIO
	{ "aio", AIO_ACCESS },
#endif
	{ "mmap", MMAP_ACCESS },
	{ "local_cache", LOCAL_CACHE_ACCESS },
	{ "global_cache", GLOBAL_CACHE_ACCESS },
	{ "parted_global", PART_GLOBAL_ACCESS },
};

str2int cache_types[] = {
	{ "tree", TREE_CACHE } ,
	{ "associative", ASSOCIATIVE_CACHE },
	{ "hash-index", HASH_INDEX_CACHE },
	{ "cuckoo", CUCKOO_CACHE },
	{ "lru2q", LRU2Q_CACHE },
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
#ifdef PROFILER
	std::string prof_file;
#endif

	if (argc < 5) {
		fprintf(stderr, "there are %d argments\n", argc);
		fprintf(stderr, "read files option pages threads cache_size entry_size preload workload cache_type num_nodes verify_content high_prio\n");
		access_map.print("available access options: ");
		workload_map.print("available workloads: ");
		cache_map.print("available cache types: ");
		exit(1);
	}

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
		}
		else if(key.compare("workload") == 0) {
			workload = workload_map.map(value);
			if (workload == -1) {
				workload_file = value;
			}
		}
		else if(key.compare("access") == 0) {
			if(value.compare("read") == 0)
				access_method = READ;
			else if(value.compare("write") == 0)
				access_method = WRITE;
			else {
				fprintf(stderr, "wrong access method\n");
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

	int num_entries = npages * (PAGE_SIZE / entry_size);

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
	for (int k = 0; k < num_files; k++)
		cnames[k] = file_names[k].c_str();
	num = num_files;
	memory_manager *manager = new memory_manager(cache_size);
	/* initialize the threads' private data. */
	for (j = 0; j < nthreads; j++) {
		switch (access_option) {
			case READ_ACCESS:
				threads[j] = new read_private(cnames, num, npages * PAGE_SIZE, j, entry_size);
				break;
			case DIRECT_ACCESS:
				threads[j] = new direct_private(cnames, num, npages * PAGE_SIZE, j, entry_size);
				break;
#if ENABLE_AIO
			case AIO_ACCESS:
				threads[j] = new aio_private(cnames, num, npages * PAGE_SIZE, j, entry_size);
				break;
#endif
			case MMAP_ACCESS:
				/* TODO for now mmap doesn't support accessing multiple files. */
				threads[j] = new mmap_private(cnames[0], j, entry_size);
				break;
			case LOCAL_CACHE_ACCESS:
				threads[j] = new part_cached_private(cnames, num, npages * PAGE_SIZE,
						j, cache_size / nthreads, entry_size);
				break;
			case GLOBAL_CACHE_ACCESS:
				threads[j] = new global_cached_private(cnames, num, npages * PAGE_SIZE, j,
						cache_size, entry_size, cache_type, manager);
				if (preload)
					((global_cached_private *) threads[j])->preload(0, npages * PAGE_SIZE);
				break;
			case PART_GLOBAL_ACCESS:
				threads[j] = new part_global_cached_private(num_nodes, cnames, num,
						npages * PAGE_SIZE, j, cache_size, entry_size, cache_type, manager);
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
		end = start + ((long) npages / nthreads + (shift < remainings)) * PAGE_SIZE / entry_size;
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
				gen = new RAID0_rand_permute_workload(npages, entry_size, nthreads, j);
				break;
			case -1:
				gen = new file_workload(workload_file, nthreads);
				break;
			default:
				fprintf(stderr, "unsupported workload\n");
				exit(1);
		}
		threads[j]->gen = gen;
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
#ifdef USE_PROCESS
		ret = process_create(&threads[i]->id, rand_read, (void *) i);
#else
		ret = pthread_create(&threads[i]->id, NULL, rand_read, (void *) (long) i);
#endif
		if (ret) {
			perror("pthread_create");
			exit(1);
		}
	}

	for (i = 0; i < nthreads; i++) {
		ssize_t size;
#ifdef USE_PROCESS
		ret = process_join(threads[i]->id);
#else
		ret = pthread_join(threads[i]->id, (void **) &size);
#endif
		if (ret) {
			perror("pthread_join");
			exit(1);
		}
		read_bytes += size;
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
	printf("%d keys are evicted from the hash table because of conflicts\n", removed_indices);
	printf("there are %d lock contentions\n", lock_contentions);
#endif
	printf("middle evicts: %d, end evicts: %d\n", middle_evicts, end_evicts);
}

const rand_permute *global_rand_permute_workload::permute;
