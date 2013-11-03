#ifndef __READ_PRIVATE_H__
#define __READ_PRIVATE_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "io_interface.h"
#include "file_partition.h"
#include "parameters.h"

class buffered_io: public io_interface
{
	logical_file_partition partition;
	/* the array of files that it's going to access */
	std::vector<int> fds;

	int flags;
	long remote_reads;
#ifdef STATISTICS
	long read_time; // in us
	long num_reads;
#endif
public:
	buffered_io(const logical_file_partition &partition_,
			thread *t, int flags = O_RDWR);

	virtual ~buffered_io() {
	}

	/* get the file descriptor corresponding to the offset. */
	int get_fd(long offset) {
		if (fds.size() == 1)
			return fds[0];

		int idx = partition.map2file(offset / PAGE_SIZE);
		return fds[idx];
	}

	const std::vector<int> &get_fds() const {
		return fds;
	}

	int num_open_files() {
		return partition.get_num_files();
	}

	const logical_file_partition &get_partition() const {
		return partition;
	}

	int get_file_id() const {
		return partition.get_file_id();
	}

	int init();

	void cleanup() {
		for (size_t i = 0; i < fds.size(); i++)
			fsync(fds[i]);
	}

	io_status access(char *buf, off_t offset, ssize_t size, int access_method);

#ifdef STATISTICS
	virtual void print_stat(int nthreads) {
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
		}
	}
#endif
};

#endif
