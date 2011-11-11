#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>

/**
 * generate a file filled with a sequence starting from 0.
 * The size of an element in the sequence is the size of long,
 * which is 8 bytes in 64-bit architecture.
 */

char buf[1024 * 1024];

int main(int argc, char *argv[])
{
	int fd;
	char *file_name;
	long size;
	char last;
	long num = 0;
	ssize_t ret;

	if (argc < 2) {
		fprintf(stderr, "create_file file_name size\n");
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
	printf("create a file of %ld bytes\n", size);

	/*
	 * because of my change in the kernel, creating a read-only file
	 * can crash the kernel because ext4 doesn't grab the semaphore
	 * when mapping blocks.
	 * So I have to make the file writable first, and then later
	 * make it read-only
	 */
	fd = open(file_name, O_WRONLY | O_CREAT, 00644);
	if (fd < 0) {
		perror("open");
		exit(1);
	}

	/* write data of `size' bytes to the file. */
	while (size > 0) {
		int i;
		int write_size = sizeof(buf);

		if (write_size > size) {
			write_size = size;
		}
		size -= write_size;
		/* generate data */
		for (i = 0; i < sizeof(buf); i += sizeof(long)) {
			*(long *) &buf[i] = num++;
		}

		ret = write(fd, buf, write_size);
		if (ret < 0) {
			perror("write");
			exit(1);
		}
	}

	return 0;
}
