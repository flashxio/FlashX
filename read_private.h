#ifndef __READ_PRIVATE_H__
#define __READ_PRIVATE_H__

#include <sys/time.h>

#include "thread_private.h"

#define MIN_BLOCK_SIZE 512

class buffered_io: public io_interface
{
	/* the array of files that it's going to access */
	const char **file_names;
	int *fds;
	/* the number of files */
	int num;
	/* the size of data it's going to access, and it'll be divided for each file */
	long size;

	int flags;
	long remote_reads;
#ifdef STATISTICS
	long read_time; // in us
	long num_reads;
#endif
public:
	buffered_io(const char *names[], int num,
			long size, int flags = O_RDONLY) {
		this->flags = flags;
#ifdef STATISTICS
		read_time = 0;
		num_reads = 0;
#endif
		remote_reads = 0;
		file_names = new const char *[num];
		for (int i = 0; i < num; i++)
			file_names[i] = names[i];
		fds = new int[num];
		this->num = num;
		this->size = size;
	}

	~buffered_io() {
		delete [] file_names;
		delete [] fds;
	}

	long get_size() {
		return size;
	}

	int get_fd_idx(long offset) {
		int fd_idx = (offset / PAGE_SIZE) % num;
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
		return fd_idx;
	}

	/* get the file descriptor corresponding to the offset. */
	int get_fd(long offset) {
		return fds[get_fd_idx(offset)];
	}

	int num_open_files() {
		return num;
	}

	int init();

	void cleanup() {
		for (int i = 0; i < num; i++)
			close(fds[i]);
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
