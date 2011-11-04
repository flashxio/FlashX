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
int flags = O_RDONLY;
int npages;
int nthreads;
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
	/* where the data read from the disk is stored */
	char *buf;
	/* shows the locations in the array where data has be to stored.*/
	off_t *buf_offset;

	/* the range in the file where we need to read data. */
	int start_i;
	int end_i;

	virtual ssize_t access(char *, off_t, ssize_t) = 0;
	virtual int thread_init() = 0;

	thread_private(int idx) {
		this->idx = idx;
		buf = (char *) valloc(PAGE_SIZE * (NUM_PAGES / nthreads));
		buf_offset = (off_t *) malloc(sizeof (*buf_offset) * NUM_PAGES / nthreads);

		if (buf == NULL){
			fprintf(stderr, "can't allocate buffer\n");
			exit(1);
		}
		/* trigger page faults and bring pages to memory. */
		for (int i = 0; i < NUM_PAGES / nthreads; i++)
			buf[i * PAGE_SIZE] = 0;

		for (int i = 0; i < NUM_PAGES / nthreads; i++)
			buf_offset[i] = i;
		permute_offset(buf_offset, NUM_PAGES / nthreads);
	}
};

class read_private: public thread_private
{
	char *file_name;
	int fd;
public:
	read_private(char *name, int idx): thread_private(idx) {
		file_name = name;
	}

	int thread_init() {
		int ret;

		fd = open(file_name, flags);
		if (fd < 0) {
			perror("open");
			exit (1);
		}
		ret = posix_fadvise(fd, start_i,
				end_i, POSIX_FADV_RANDOM);
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

class mmap_private: public thread_private
{
	void *addr;

public:
	mmap_private(char *new_name, int idx): thread_private(idx) {
		static void *addr = NULL;
		static char *file_name = NULL;
		/* if we are mapping to a different file, do the real mapping. */
		if (file_name == NULL || strcmp(file_name, new_name)) {
			int fd = open(new_name, flags);
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

class part_cached_private: public read_private
{
	page_cache *cache;

public:
	part_cached_private(char *name, int idx): read_private(name, idx) { }

	ssize_t access(char *buf, off_t offset, ssize_t size) {
		ssize_t ret = cache->get_from_cache(buf, offset, size);
		if (ret < 0) {
			struct page *p = cache->get_empty_page(offset);
			ret = read_private::access((char *) p->get_data(),
					ROUND_PAGE(offset), PAGE_SIZE);
			if (ret < 0) {
				perror("read");
				exit(1);
			}
			p->set_offset(offset);
			offset -= ROUND_PAGE(offset);
			/* I assume the data I read never crosses the page boundary */
			memcpy(buf, (char *) p->get_data() + offset, size);
		}
		else
			cache_hits++;
		return ret;
	}
};

class global_cached_private: public read_private
{
	page_cache *cache;
public:
	global_cached_private(char *name, int idx): read_private(name, idx) { }

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
	char *buf;
	off_t *buf_offset;
	int fd;

	priv->thread_init();
	start_i = priv->start_i;
	end_i = priv->end_i;
	buf = priv->buf;
	buf_offset = priv->buf_offset;

	gettimeofday(&start_time, NULL);
	for (i = start_i, j = 0; i < end_i; i++, j++) {
		if (j == NUM_PAGES / nthreads)
			j = 0;
		ret = priv->access(buf + buf_offset[j] * PAGE_SIZE,
				offset[i], PAGE_SIZE);
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
	MMAP_ACCESS,
	LOCAL_CACHE_ACCESS,
	GLOBAL_CACHE_ACCESS
};

int main(int argc, char *argv[])
{
	int cache_size = 512 * 1024 * 1024;
	int access_options;
	int ret;
	int i, j;
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;
	int is_mmap = 0;
	int num_files = 0;
	char *file_names[NUM_THREADS];

	if (argc < 5) {
		fprintf(stderr, "read files option num_pages num_threads\n");
		exit(1);
	}

	switch (atoi(argv[argc - 3])) {
		case NORMAL:
			access_options = READ_ACCESS;
			break;
		case DIRECT:
			flags |= O_DIRECT;
			access_options = READ_ACCESS;
			break;
		case MMAP:
			is_mmap = 1;
			access_options = MMAP_ACCESS;
			break;
		default:
			fprintf(stderr, "wrong option\n");
			exit(1);
	}
	num_files = argc - 4;
	for (i = 0; i < num_files; i++) {
		file_names[i] = argv[1 + i];
	}

	npages = atoi(argv[argc - 2]);
	offset = (off_t *) malloc(sizeof(*offset) * npages);
	for(i = 0; i < npages; i++) {
		offset[i] = ((off_t) i) * 4096L;
		if (offset[i] < 0) {
			printf("offset[%d]: %ld\n", i, offset[i]);
			exit(1);
		}
	}
	permute_offset(offset, npages);

	nthreads = atoi(argv[argc - 1]);
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
		char *file_name;
		if (num_files > 1) {
			file_name = file_names[j];
		}
		else {
			file_name = file_names[0];
		}
		switch (access_options) {
			case READ_ACCESS:
				threads[j] = new read_private(file_name, j);
				break;
			case MMAP_ACCESS:
				threads[j] = new mmap_private(file_name, j);
				break;
			case LOCAL_CACHE_ACCESS:
				threads[j] = new part_cached_private(file_name, j);
				break;
			case GLOBAL_CACHE_ACCESS:
				threads[j] = new global_cached_private(file_name, j);
				break;
			default:
				fprintf(stderr, "wrong access option\n");
				exit(1);
		}
		
		if (num_files > 1) {
			threads[j]->start_i = 0;
			threads[j]->end_i = npages;
		}
		else {
			threads[j]->start_i = npages / nthreads * j;
			threads[j]->end_i = threads[j]->start_i + npages / nthreads;
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
