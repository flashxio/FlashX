#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <numaif.h>
#include <numa.h>
#include <sys/time.h>

#include "common.h"

/**
 * generate a file filled with a sequence starting from 0.
 * The size of an element in the sequence is the size of long,
 * which is 8 bytes in 64-bit architecture.
 */

char buf[4096]__attribute__((aligned(4096)));
#define MAX_NODES 32

int main(int argc, char *argv[])
{
	int fd;
	char *file_name;
	long size;
	char last;

	if (argc < 3) {
		fprintf(stderr, "verify_file file_name size\n");
		exit(1);
	}

	file_name = argv[1];
	last = argv[2][strlen(argv[2]) - 1];
	argv[2][strlen(argv[2]) - 1] = 0;
	size = atol(argv[2]);
	switch (last) {
		case 'G':
		case 'g':
			size *= 1024 * 1024 * 1024;
			break;
		case 'M':
		case 'm':
			size *= 1024 * 1024;
			break;
		case 'K':
		case 'k':
			size *= 1024;
			break;
	}
	printf("verify a file of %ld bytes\n", size);

	/*
	 * because of my change in the kernel, creating a read-only file
	 * can crash the kernel because ext4 doesn't grab the semaphore
	 * when mapping blocks.
	 * So I have to make the file writable first, and then later
	 * make it read-only
	 */
	fd = open(file_name, O_DIRECT | O_RDONLY);
	if (fd < 0) {
		perror("open");
		exit(1);
	}

	/* write data of `size' bytes to the file. */
	struct timeval start_time, end_time;
	gettimeofday(&start_time, NULL);
	size_t tot_size = 0;
	for (off_t off = 0; off < size; off += sizeof(buf)) {
		ssize_t ret = pread(fd, buf, sizeof(buf), off);
		assert(ret == sizeof(buf));
		tot_size += ret;
		for (int i = 0; i < sizeof(buf); i += sizeof(long))
			printf("%ld\n", *(long *) (buf + i));
	}
	gettimeofday(&end_time, NULL);
	printf("read %ld bytes, takes %f seconds\n",
			tot_size, end_time.tv_sec - start_time.tv_sec
			+ ((float)(end_time.tv_usec - start_time.tv_usec))/1000000);

	return 0;
}
