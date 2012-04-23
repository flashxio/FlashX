#ifndef __READ_PRIVATE_H__
#define __READ_PRIVATE_H__

#include <sys/time.h>

#include "thread_private.h"

class read_private: public thread_private
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
	read_private(const char *names[], int num, long size, int idx, int entry_size,
			int flags = O_RDWR): thread_private(idx, entry_size) {
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

	~read_private() {
		delete [] file_names;
		delete [] fds;
		delete buf;
	}

	long get_size() {
		return size;
	}

	/* get the file descriptor corresponding to the offset. */
	int get_fd(long offset) {
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
		return fds[fd_idx];
	}

	int thread_init();

	int thread_end() {
		for (int i = 0; i < num; i++)
			close(fds[i]);
		return 0;
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method);

#ifdef STATISTICS
	virtual void print_stat() {
		thread_private::print_stat();
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
