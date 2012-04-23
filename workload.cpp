#include "workload.h"

off_t *file_workload::offsets;

off_t stride_workload::next_offset() {
	off_t ret = curr;
	num++; 

	/*
	 * we stride with PAGE_SIZE.
	 * When we reach the end of the range,
	 * we start over but move one ahead from the last startover.
	 */
	curr += stride;
	if (curr >= last) {
		curr = first + (curr & (stride - 1));
		curr++;
	}
	ret *= entry_size;
	return ret;
}

file_workload::file_workload(const std::string &file, int nthreads) {
	static off_t file_size;
	static int remainings;
	static int shift = 0;
	static long start;
	static long end = 0;

	if (offsets == NULL) {
		int fd = open(file.c_str(), O_RDONLY);
		if (fd < 0) {
			perror("open");
			exit(1);
		}
		printf("%s's fd is %d\n", file.c_str(), fd);

		/* get the file size */
		struct stat stats;
		if (fstat(fd, &stats) < 0) {
			perror("fstat");
			exit(1);
		}
		file_size = stats.st_size;
		remainings = file_size / sizeof(off_t) % nthreads;
		assert(file_size % sizeof(off_t) == 0);

		offsets = (off_t *) malloc(file_size);
		/* read data of the file to a buffer */
		char *buf = (char *) offsets;
		long size = file_size;
		while (size > 0) {
			ssize_t ret = read(fd, buf, size);
			if (ret < 0) {
				perror("read");
				exit(1);
			}
			buf += ret;
			size -= ret;
		}
		close(fd);
	}

	/* the range in `offsets' */
	start = end;
	end = start + file_size / sizeof(off_t) / nthreads + (shift < remainings);
	this->curr = start;
	this->end = end;
	printf("start at %ld end at %ld\n", curr, end);
}

bool stride_workload_chunk::get_workload(off_t *offsets, int num) {
	long start;
	long end;

	pthread_spin_lock(&_lock);
	start = curr;
	curr += stride * num;
	end = curr;
	/*
	 * if the chunk we try to get is in the range,
	 * get the chunk. 
	 */
	if (end < last + stride)
		goto unlock;

	/*
	 * the chunk is out of the range,
	 * let's start over but move the first entry forward.
	 */
	curr = first + (curr & (stride - 1));
	curr++;
	/*
	 * if the first entry is in the second page,
	 * it means we have accessed all pages, so no more work to do.
	 */
	if (curr == first + stride) {
		pthread_spin_unlock(&_lock);
		curr = end;
		return false;
	}
	start = curr;
	curr += stride * num;
	end = curr;
unlock:
	pthread_spin_unlock(&_lock);

	for (long i = 0; start < end; i++, start += stride)
		offsets[i] = start * entry_size;
	return true;
}

workload_chunk *balanced_workload::chunks;
local_rand_permute_workload *RAID0_rand_permute_workload::gen;
