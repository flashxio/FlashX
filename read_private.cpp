#include "read_private.h"

buffered_io::buffered_io(const logical_file_partition &partition_, long size,
		int node_id, int flags): io_interface(node_id), partition(
			partition_), fds(partition.get_num_files())
{
	this->flags = flags;
#ifdef STATISTICS
	read_time = 0;
	num_reads = 0;
#endif
	remote_reads = 0;
	this->size = size;

	for (int i = 0; i < partition.get_num_files(); i++)
		fds[i] = partition.get_fd(i);
}

int buffered_io::init() {
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
