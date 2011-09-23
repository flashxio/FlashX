#define _GNU_SOURCE
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <stdlib.h>
#include <sys/resource.h>

#define NUM_THREADS 8
#define PAGE_SIZE 4096

off_t *offset;
int flags = O_RDONLY;
int npages;
char *file_name;
struct timeval global_start;
char static_buf[PAGE_SIZE * 8] __attribute__((aligned(PAGE_SIZE)));

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
	pthread_t id;
	int idx;
	char *file_name;
	/* where the data read from the disk is stored */
	char *buf;
	/* shows the locations in the array where data has be to stored.*/
	off_t *buf_offset;

	/* the range in the file where we need to read data. */
	int start_i;
	int end_i;
};

void rand_read(void *arg)
{
#define NUM_PAGES 16384
	ssize_t ret;
	int i, j, start_i, end_i;
	ssize_t read_bytes = 0;
	struct timeval start_time, end_time;
	struct thread_private *private = arg;
	char *buf;
	off_t *buf_offset;
	int fd;
	int idx;
	
	fd = open(private->file_name, flags);
	if (fd < 0) {
		perror("open");
		exit (1);
	}

	start_i = private->start_i;
	end_i = private->end_i;
	buf = private->buf;
	buf_offset = private->buf_offset;
	idx = private->idx;

	gettimeofday(&start_time, NULL);
	for (i = start_i, j = 0; i < end_i; i++, j++) {
		ret = lseek(fd, offset[i], SEEK_SET);
		if (ret < 0) {
			perror("lseek");
			printf("%ld\n", offset[i]);
			exit(1);
		}
		if (j == NUM_PAGES)
			j = 0;
		ret = read(fd, buf + buf_offset[j] * PAGE_SIZE, PAGE_SIZE);
//		ret = read(fd, static_buf + idx * PAGE_SIZE, PAGE_SIZE);
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
	
	pthread_exit((void *) read_bytes);
}

int main(int argc, char *argv[])
{
	int ret;
	int i, j;
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;
	struct thread_private threads[NUM_THREADS];
	int nfiles;

	if (argc < 4) {
		fprintf(stderr, "read files is_direct num_pages\n");
		exit(1);
	}

	nfiles = argc - 3;

	if (atoi(argv[argc - 2]) != 0)
		flags |= O_DIRECT;

	npages = atoi(argv[argc - 1]);
	offset = malloc(sizeof(*offset) * npages);
	for(i = 0; i < npages; i++) {
		offset[i] = ((off_t) i) * 4096L;
		if (offset[i] < 0) {
			printf("offset[%d]: %ld\n", i, offset[i]);
			exit(1);
		}
	}
	permute_offset(offset, npages);

	/* initialize the threads' private data.
	 * a thread is created to read one file. */
	for (j = 0; j < nfiles; j++) {
		char *buf;
		off_t *buf_offset;
		buf = valloc(PAGE_SIZE * (NUM_PAGES));
		buf_offset = malloc(sizeof (*buf_offset) * NUM_PAGES);

		if (buf == NULL){
			fprintf(stderr, "can't allocate buffer\n");
			exit(1);
		}
		/* trigger page faults and bring pages to memory. */
		for (i = 0; i < NUM_PAGES; i++)
			buf[i * PAGE_SIZE] = 0;

		for (i = 0; i < NUM_PAGES; i++)
			buf_offset[i] = i;
		permute_offset(buf_offset, NUM_PAGES);
		
		threads[j].file_name = argv[j + 1];
		threads[j].idx = j;
		threads[j].buf = buf;
		threads[j].buf_offset = buf_offset;
		threads[j].start_i = 0;
		threads[j].end_i = npages;
	}

	ret = setpriority(PRIO_PROCESS, getpid(), -20);
	if (ret < 0) {
		perror("setpriority");
		exit(1);
	}

	gettimeofday(&start_time, NULL);
	global_start = start_time;
	for (i = 0; i < nfiles; i++) {
		ret = pthread_create(&threads[i].id, NULL, rand_read, &threads[i]);
		if (ret) {
			perror("pthread_create");
			exit(1);
		}
	}

	for (i = 0; i < nfiles; i++) {
		ssize_t size;
		ret = pthread_join(threads[i].id, (void **) &size);
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
