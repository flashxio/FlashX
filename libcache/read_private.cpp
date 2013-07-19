#include "read_private.h"
#include "file_mapper.h"

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

	for (int i = 0; i < partition.get_num_files(); i++) {
		int ret;
		fds[i] = open(partition.get_file_name(i).c_str(), flags);
		if (fds[i] < 0) {
			char err_msg[128];
			snprintf(err_msg, sizeof(err_msg),
					"open %s\n", partition.get_file_name(i).c_str());
			perror(err_msg);
			exit (1);
		}
		ret = posix_fadvise(fds[i], 0, 0, POSIX_FADV_RANDOM);
		if (ret < 0) {
			perror("posix_fadvise");
			exit(1);
		}
	}
}

int buffered_io::init() {
	return 0;
}

io_status buffered_io::access(char *buf, off_t offset, ssize_t size, int access_method) {
	int fd;
	if (fds.size() == 1)
		fd = fds[0];
	else {
		struct block_identifier bid;
		partition.map(offset / PAGE_SIZE, bid);
		fd = fds[bid.idx];
		offset = bid.off * PAGE_SIZE;
	}
#ifdef STATISTICS
	if (access_method == READ)
		num_reads++;
	struct timeval start, end;
	gettimeofday(&start, NULL);
#endif
	ssize_t ret;
	// TODO I need to make sure all data is read or written to the file.
	if (access_method == WRITE)
		ret = pwrite(fd, buf, size, offset);
	else
		ret = pread(fd, buf, size, offset);
#ifdef STATISTICS
	gettimeofday(&end, NULL);
	read_time += ((long) end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec;
#endif
	if (ret < 0)
		return IO_FAIL;
	else
		return IO_OK;
}
