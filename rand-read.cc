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

#include <iostream>
#include <string>

#include "cache.h"

//#define USE_PROCESS

#define NUM_PAGES 16384
#define NUM_THREADS 32
enum {
	NORMAL,
	DIRECT,
	MMAP,
};

int cache_hits;
off_t *offset;
int npages;
int nthreads = 1;
struct timeval global_start;
char static_buf[PAGE_SIZE * 8] __attribute__((aligned(PAGE_SIZE)));
volatile int first[NUM_THREADS];

void permute_offset(off_t *offset, int num)
{
	int i;
	for (i = num - 1; i >= 1; i--) {
		int j = random() % i;
		off_t tmp = offset[j];
		offset[j] = offset[i];
		offset[i] = tmp;
	}
}

float time_diff(struct timeval time1, struct timeval time2)
{
	return time2.tv_sec - time1.tv_sec
			+ ((float)(time2.tv_usec - time1.tv_usec))/1000000;
}

class rand_buf
{
	/* where the data read from the disk is stored */
	char *buf;
	/* shows the locations in the array where data has be to stored.*/
	off_t *buf_offset;
	int entry_size;
	int num_entries;

	int current;
public:
	rand_buf(int buf_size, int entry_size) {
		this->entry_size = entry_size;
		num_entries = buf_size / entry_size;
		buf = (char *) valloc(buf_size);
		buf_offset = (off_t *) malloc(sizeof (*buf_offset) * num_entries);

		if (buf == NULL){
			fprintf(stderr, "can't allocate buffer\n");
			exit(1);
		}
		/* trigger page faults and bring pages to memory. */
		for (int i = 0; i < buf_size / PAGE_SIZE; i++)
			buf[i * PAGE_SIZE] = 0;

		for (int i = 0; i < num_entries; i++)
			buf_offset[i] = i * entry_size;
		permute_offset(buf_offset, num_entries);

		current = 0;
	}

	~rand_buf() {
		free(buf);
		free(buf_offset);
	}

	char *next_entry() {
		int off = buf_offset[current];
		current = (current + 1) % num_entries;;
		return &buf[off];
	}

	int get_entry_size() {
		return entry_size;
	}
};

/* this data structure stores the thread-private info. */
class thread_private
{
public:
#ifdef USE_PROCESS
	pid_t id;
#else
	pthread_t id;
#endif
	/* the location in the thread descriptor array. */
	int idx;
	rand_buf buf;

	/* the range in the file where we need to read data. */
	int start_i;	// the first entry in the range
	int end_i;		// the last entry in the range

	virtual ssize_t access(char *, off_t, ssize_t) = 0;
	virtual int thread_init() = 0;

	thread_private(int idx, int entry_size): buf(NUM_PAGES / nthreads * PAGE_SIZE, entry_size) {
		this->idx = idx;
	}
};

class read_private: public thread_private
{
	const char *file_name;
	int fd;
	int flags;
public:
	read_private(const char *name, int idx, int entry_size,
			int flags = O_RDONLY): thread_private(idx, entry_size), file_name(name) {
		this->flags = flags;
	}

	int thread_init() {
		int ret;

		fd = open(file_name, flags);
		if (fd < 0) {
			perror("open");
			exit (1);
		}
//		ret = posix_fadvise(fd, (off_t) start_i * buf.get_entry_size(),
//				(off_t) (end_i - start_i) * buf.get_entry_size(), POSIX_FADV_RANDOM);
		ret = posix_fadvise(fd, 0, 0, POSIX_FADV_RANDOM);
		if (ret < 0) {
			perror("posix_fadvise");
			exit(1);
		}
		return 0;
	}

	ssize_t access(char *buf, off_t offset, ssize_t size) {
		off_t ret = lseek(fd, offset, SEEK_SET);
		if (ret < 0) {
			perror("lseek");
			printf("%ld\n", offset);
			exit(1);
		}
		ret = read(fd, buf, size);
		return ret;
	}
};

class direct_private: public read_private
{
public:
	direct_private(const char *name, int idx, int entry_size): read_private(name, idx,
			entry_size, O_DIRECT | O_RDONLY) {}
};

class mmap_private: public thread_private
{
	void *addr;

public:
	mmap_private(const char *new_name, int idx,
			int entry_size): thread_private(idx, entry_size) {
		static void *addr = NULL;
		static const char *file_name = NULL;
		/* if we are mapping to a different file, do the real mapping. */
		if (file_name == NULL || strcmp(file_name, new_name)) {
			int fd = open(new_name, O_RDONLY);
			int ret;

			if (fd < 0) {
				perror("open");
				exit (1);
			}
			addr = mmap(NULL, ((ssize_t) npages) * PAGE_SIZE,
					PROT_READ, MAP_PRIVATE, fd, 0);
			if (addr == NULL) {
				perror("mmap");
				exit(1);
			}
			ret = madvise(addr, ((ssize_t) npages) * PAGE_SIZE, MADV_RANDOM);
			if (ret < 0) {
				perror("madvise");
				exit(1);
			}
			file_name = new_name;
		}
		this->addr = addr;
	}

	int thread_init() {
		return 0;
	}

	ssize_t access(char *buf, off_t offset, ssize_t size) {
		int i;
		char *page = (char *) addr + offset;
		/* I try to avoid gcc optimization eliminating the code below,
		 * and it works. */
		first[idx] = page[0];
		__asm__ __volatile__("" : : "g"(&first[idx]));
		return size;
	}
};

class part_cached_private: public direct_private
{
	page_cache *cache;

public:
	part_cached_private(const char *name, int idx, long cache_size,
			int entry_size): direct_private(name, idx, entry_size) {
		cache = new tree_cache(cache_size);
	}

	ssize_t access(char *buf, off_t offset, ssize_t size) {
		ssize_t ret = cache->get_from_cache(buf, offset, size);
		if (ret < 0) {
			struct page *p = cache->get_empty_page(ROUND_PAGE(offset));
			ret = read_private::access((char *) p->get_data(),
					ROUND_PAGE(offset), PAGE_SIZE);
			if (ret < 0) {
				perror("read");
				exit(1);
			}
			offset -= ROUND_PAGE(offset);
			/* I assume the data I read never crosses the page boundary */
			memcpy(buf, (char *) p->get_data() + offset, size);
			ret = size;
		}
		else
			cache_hits++;
		return ret;
	}
};

class global_cached_private: public direct_private
{
	page_cache *cache;
public:
	global_cached_private(const char *name, int idx,
			int entry_size): direct_private(name, idx, entry_size) { }

	ssize_t access(char *buf, off_t offset, ssize_t size) {
		return 0;
	}
};

thread_private *threads[NUM_THREADS];

void *rand_read(void *arg)
{
	ssize_t ret;
	int i, j, start_i, end_i;
	ssize_t read_bytes = 0;
	struct timeval start_time, end_time;
	thread_private *priv = threads[(long) arg];
	rand_buf *buf;
	int fd;

	priv->thread_init();
	start_i = priv->start_i;
	end_i = priv->end_i;
	buf = &priv->buf;

	gettimeofday(&start_time, NULL);
	for (i = start_i, j = 0; i < end_i; i++, j++) {
		ret = priv->access(buf->next_entry(),
				offset[i], buf->get_entry_size());
		assert(ret == buf->get_entry_size());
		if (ret > 0)
			read_bytes += ret;
		else
			break;
	}
	if (ret < 0) {
		perror("read");
		exit(1);
	}
	gettimeofday(&end_time, NULL);
	printf("read %ld bytes, start at %f seconds, takes %f seconds\n",
			read_bytes, time_diff(global_start, start_time),
			time_diff(start_time, end_time));
	
#ifdef USE_PROCESS
	exit(read_bytes);
#else
	pthread_exit((void *) read_bytes);
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

enum {
	READ_ACCESS,
	DIRECT_ACCESS,
	MMAP_ACCESS,
	LOCAL_CACHE_ACCESS,
	GLOBAL_CACHE_ACCESS
};

struct access_method {
	std::string name;
	int access;
} methods[] = {
	"normal", READ_ACCESS,
	"direct", DIRECT_ACCESS,
	"mmap", MMAP_ACCESS,
	"local_cache", LOCAL_CACHE_ACCESS,
	"global_cache", GLOBAL_CACHE_ACCESS
};

void print_methods() {
	std::cout<<"available methods: ";
	int num_methods = sizeof(methods) / sizeof(methods[0]);
	for (int i = 0; i < num_methods; i++) {
		std::cout<<methods[i].name;
		if (i < num_methods - 1)
			std::cout<<", ";
	}
	std::cout<<std::endl;
}

int map_method(const std::string &method)
{
	int num_methods = sizeof(methods) / sizeof(methods[0]);
	for (int i = 0; i < num_methods; i++) {
		if (methods[i].name.compare(method) == 0)
			return i;
	}
	return -1;
}

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
	long cache_size = 512 * 1024 * 1024;
	int entry_size = 128;
	int access_option;
	int ret;
	int i, j;
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;
	int num_files = 0;
	std::string file_names[NUM_THREADS];

	if (argc < 5) {
		fprintf(stderr, "read files option pages threads cache_size entry_size\n");
		print_methods();
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
			access_option = map_method(value);
			if (access_option < 0) {
				fprintf(stderr, "can't find the right option\n");
				exit(1);
			}
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
		else if(key.compare("entry_size") == 0) {
			entry_size = (int) str2size(value);
		}
		else {
			fprintf(stderr, "wrong option\n");
			exit(1);
		}
	}

	int num_entries = npages * (PAGE_SIZE / entry_size);
	offset = (off_t *) malloc(sizeof(*offset) * num_entries);
	for(i = 0; i < num_entries; i++) {
		offset[i] = ((off_t) i) * entry_size;
		if (offset[i] < 0) {
			printf("offset[%d]: %ld\n", i, offset[i]);
			exit(1);
		}
	}
	permute_offset(offset, num_entries);

	if (nthreads > NUM_THREADS) {
		fprintf(stderr, "too many threads\n");
		exit(1);
	}
	if (num_files > 1 && num_files != nthreads) {
		fprintf(stderr, "if there are multiple files, \
				the number of files must be the same as the number of threads\n");
		exit(1);
	}

	/* initialize the threads' private data. */
	for (j = 0; j < nthreads; j++) {
		const char *file_name;
		if (num_files > 1) {
			file_name = file_names[j].c_str();
		}
		else {
			file_name = file_names[0].c_str();
		}
		switch (access_option) {
			case READ_ACCESS:
				threads[j] = new read_private(file_name, j, entry_size);
				break;
			case DIRECT_ACCESS:
				threads[j] = new direct_private(file_name, j, entry_size);
				break;
			case MMAP_ACCESS:
				threads[j] = new mmap_private(file_name, j, entry_size);
				break;
			case LOCAL_CACHE_ACCESS:
				threads[j] = new part_cached_private(file_name, j, cache_size / nthreads, entry_size);
				break;
			case GLOBAL_CACHE_ACCESS:
				threads[j] = new global_cached_private(file_name, j, entry_size);
				break;
			default:
				fprintf(stderr, "wrong access option\n");
				exit(1);
		}
		
		if (num_files > 1) {
			threads[j]->start_i = 0;
			threads[j]->end_i = npages * PAGE_SIZE / entry_size;
		}
		else {
			threads[j]->start_i = npages / nthreads * PAGE_SIZE / entry_size * j;
			threads[j]->end_i = threads[j]->start_i + npages / nthreads * PAGE_SIZE / entry_size;
		}
	}

	ret = setpriority(PRIO_PROCESS, getpid(), -20);
	if (ret < 0) {
		perror("setpriority");
		exit(1);
	}

	gettimeofday(&start_time, NULL);
	global_start = start_time;
	for (i = 0; i < nthreads; i++) {
#ifdef USE_PROCESS
		ret = process_create(&threads[i]->id, rand_read, (void *) i);
#else
		ret = pthread_create(&threads[i]->id, NULL, rand_read, (void *) i);
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
	gettimeofday(&end_time, NULL);
	printf("read %ld bytes, takes %f seconds\n",
			read_bytes, end_time.tv_sec - start_time.tv_sec
			+ ((float)(end_time.tv_usec - start_time.tv_usec))/1000000);
	printf("there are %d cache hits\n", cache_hits);
}
