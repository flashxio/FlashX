#ifndef __READ_PRIVATE_H__
#define __READ_PRIVATE_H__

#include <sys/time.h>

#include "thread_private.h"
#include "file_partition.h"
#include "parameters.h"

class buffered_io: public io_interface
{
	logical_file_partition partition;
	/* the array of files that it's going to access */
	std::vector<int> fds;
	/* the size of data it's going to access, and it'll be divided for each file */
	long size;

	int flags;
	long remote_reads;
#ifdef STATISTICS
	long read_time; // in us
	long num_reads;
#endif
public:
	buffered_io(const logical_file_partition &partition_, long size,
			int node_id, int flags = O_RDWR);

	virtual ~buffered_io() {
	}

	long get_size() const {
		return size;
	}

	long get_local_size() const {
		// TODO I assume all files are connected to the local node.
		return size;
	}

	/* get the file descriptor corresponding to the offset. */
	int get_fd(long offset) {
		if (fds.size() == 1)
			return fds[0];

		int idx = partition.map2file(offset / PAGE_SIZE);
		return fds[idx];
	}

	int num_open_files() {
		return partition.get_num_files();
	}

	const logical_file_partition &get_partition() const {
		return partition;
	}

	int init();

	void cleanup() {
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method);

#ifdef STATISTICS
	virtual void print_stat() {
		static int seen_threads = 0;
		static long tot_nreads;
		static long tot_read_time;
		static long tot_remote_reads;
		tot_remote_reads += remote_reads;
		tot_nreads += num_reads;
		tot_read_time += read_time;
		seen_threads++;
		if (seen_threads == nthreads) {
			printf("there are %ld reads and takes %ldus\n", tot_nreads, tot_read_time);
#if NUM_NODES > 1
			printf("total remote reads: %ld\n", tot_remote_reads);
#endif
		}
	}
#endif
};

#endif
