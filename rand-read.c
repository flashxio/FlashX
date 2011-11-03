#define _GNU_SOURCE
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include <string.h>
#include <errno.h>

#define USE_PROCESS

#define NUM_PAGES 16384
#define NUM_THREADS 32
#define PAGE_SIZE 4096
enum {
	NORMAL,
	DIRECT,
	MMAP,
};

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
struct thread_private
{
#ifdef USE_PROCESS
	pid_t id;
#else
	pthread_t id;
#endif
	int idx;
	char *file_name;
	int fd;
	void *addr;
	/* where the data read from the disk is stored */
	char *buf;
	/* shows the locations in the array where data has be to stored.*/
	off_t *buf_offset;

	/* the range in the file where we need to read data. */
	int start_i;
	int end_i;

	int (*thread_init) (struct thread_private *);
	ssize_t (*access) (struct thread_private *, char *, off_t, ssize_t);
};

int single_file_thread_init(struct thread_private *private)
{
	int ret;

	private->fd = open(private->file_name, flags);
	if (private->fd < 0) {
		perror("open");
		exit (1);
	}
	ret = posix_fadvise(private->fd, private->start_i,
			private->end_i, POSIX_FADV_RANDOM);
	if (ret < 0) {
		perror("posix_fadvise");
		exit(1);
	}
	return 0;
}

ssize_t single_file_access(struct thread_private *private, char *buf,
		off_t offset, ssize_t size)
{
	off_t ret = lseek(private->fd, offset, SEEK_SET);
	if (ret < 0) {
		perror("lseek");
		printf("%ld\n", offset);
		exit(1);
	}
	ret = read(private->fd, buf, size);
//	ret = read(fd, static_buf + idx * PAGE_SIZE, PAGE_SIZE);
	return ret;
}

ssize_t mmap_access(struct thread_private *private, char *buf,
		off_t offset, ssize_t size)
{
	int i;
	char *addr = private->addr + offset;
	/* I try to avoid gcc optimization eliminating the code below,
	 * and it works. */
	first[private->idx] = addr[0];
	__asm__ __volatile__("" : : "g"(&first[private->idx]));
	return size;
}

void rand_read(void *arg)
{
	ssize_t ret;
	int i, j, start_i, end_i;
	ssize_t read_bytes = 0;
	struct timeval start_time, end_time;
	struct thread_private *private = arg;
	char *buf;
	off_t *buf_offset;
	int fd;
	int idx;

	if (private->thread_init)
		private->thread_init(private);
	start_i = private->start_i;
	end_i = private->end_i;
	buf = private->buf;
	buf_offset = private->buf_offset;
	idx = private->idx;

	gettimeofday(&start_time, NULL);
	for (i = start_i, j = 0; i < end_i; i++, j++) {
		if (j == NUM_PAGES / nthreads)
			j = 0;
		ret = private->access(private, buf + buf_offset[j] * PAGE_SIZE,
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

int process_create(pid_t *pid, void (*func)(void *), void *private)
{
	pid_t id = fork();

	if (id < 0)
		return -1;

	if (id == 0) {	// child
		func(private);
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

int main(int argc, char *argv[])
{
	int ret;
	int i, j;
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;
	struct thread_private threads[NUM_THREADS];
	int is_mmap = 0;
	int num_files = 0;
	char *file_names[NUM_THREADS];

	if (argc < 5) {
		fprintf(stderr, "read files option num_pages num_threads\n");
		exit(1);
	}

	switch (atoi(argv[argc - 3])) {
		case NORMAL:
			break;
		case DIRECT:
			flags |= O_DIRECT;
			break;
		case MMAP:
			is_mmap = 1;
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
	offset = malloc(sizeof(*offset) * npages);
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
		char *buf;
		off_t *buf_offset;
		buf = valloc(PAGE_SIZE * (NUM_PAGES / nthreads));
		buf_offset = malloc(sizeof (*buf_offset) * NUM_PAGES / nthreads);

		if (buf == NULL){
			fprintf(stderr, "can't allocate buffer\n");
			exit(1);
		}
		/* trigger page faults and bring pages to memory. */
		for (i = 0; i < NUM_PAGES / nthreads; i++)
			buf[i * PAGE_SIZE] = 0;

		for (i = 0; i < NUM_PAGES / nthreads; i++)
			buf_offset[i] = i;
		permute_offset(buf_offset, NUM_PAGES / nthreads);
		
		threads[j].idx = j;
		threads[j].buf = buf;
		threads[j].buf_offset = buf_offset;
		if (num_files > 1) {
			threads[j].file_name = file_names[j];
			threads[j].start_i = 0;
			threads[j].end_i = npages;
		}
		else {
			threads[j].file_name = file_names[0];
			threads[j].start_i = npages / nthreads * j;
			threads[j].end_i = threads[j].start_i + npages / nthreads;
		}
		if (is_mmap) {
			static void *addr = NULL;
			static char *file_name = NULL;
			/* if we are mapping to a different file, do the real mapping. */
			if (file_name == NULL || strcmp(file_name, threads[j].file_name)) {
				int fd = open(threads[j].file_name, flags);
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
				file_name == threads[j].file_name;
			}
			threads[j].thread_init = NULL;
			threads[j].access = mmap_access;
			threads[j].addr = addr;
		}
		else {
			threads[j].thread_init = single_file_thread_init;
			threads[j].access = single_file_access;
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
		ret = process_create(&threads[i].id, rand_read, &threads[i]);
#else
		ret = pthread_create(&threads[i].id, NULL, rand_read, &threads[i]);
#endif
		if (ret) {
			perror("pthread_create");
			exit(1);
		}
	}

	for (i = 0; i < nthreads; i++) {
		ssize_t size;
#ifdef USE_PROCESS
		ret = process_join(threads[i].id);
#else
		ret = pthread_join(threads[i].id, (void **) &size);
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
}
