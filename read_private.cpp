#include "read_private.h"

int read_private::thread_init() {
	int ret;

	buf = new rand_buf(NUM_PAGES / (nthreads / NUM_NODES) * PAGE_SIZE, get_entry_size());
	for (int i = 0; i < num; i++) {
		fds[i] = open(file_names[i], flags);
		if (fds[i] < 0) {
			perror("open");
			exit (1);
		}
		ret = posix_fadvise(fds[i], 0, 0, POSIX_FADV_RANDOM);
		if (ret < 0) {
			perror("posix_fadvise");
			exit(1);
		}
	}
	return 0;
}

ssize_t read_private::access(char *buf, off_t offset, ssize_t size, int access_method) {
	int fd_idx = offset / (this->size / num);
	if (fd_idx >= num) {
		printf("offset: %ld, fd_idx: %d, size: %ld, num: %d\n", offset, fd_idx, this->size, num);
	}
#if NUM_NODES > 1
	int node_num = idx / (nthreads / NUM_NODES);
	if (node_num != fd_idx) {
		remote_reads++;
	}
#endif
	assert (fd_idx < num);
	int fd = fds[fd_idx];
#ifdef STATISTICS
	if (access_method == READ)
		num_reads++;
	struct timeval start, end;
	gettimeofday(&start, NULL);
#endif
	ssize_t ret;
	if (access_method == WRITE)
		ret = pwrite(fd, buf, size, offset);
	else
		ret = pread(fd, buf, size, offset);
#ifdef STATISTICS
	gettimeofday(&end, NULL);
	read_time += ((long) end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec;
#endif
	return ret;
}
