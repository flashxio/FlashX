#include "read_private.h"

int buffered_io::init() {
	int ret;

	for (int i = 0; i < partition.get_num_files(); i++) {
		fds[i] = open(partition.get_file_name(i).c_str(), flags);
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

ssize_t buffered_io::access(char *buf, off_t offset, ssize_t size, int access_method) {
	int fd = get_fd(offset);
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
	if (ret < 0) {
		perror("pread");
		abort();
	}
#ifdef STATISTICS
	gettimeofday(&end, NULL);
	read_time += ((long) end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec;
#endif
	return ret;
}
