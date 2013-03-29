#include "workload.h"
#include "common.h"
#include "container.cpp"

int workload_gen::default_entry_size = PAGE_SIZE;
int workload_gen::default_access_method = READ;

off_t cache_hit_defined_workload::next_offset()
{
	if (seq < cache_hit_ratio * 100 && cached_pages.size() > 0) {
		// cache hit
		seq = (seq + 1) % 100;
		off_t ret = cached_pages[cache_hit_seq];
		cache_hit_seq = (cache_hit_seq + 1) % cached_pages.size();
		return ret;
	}
	else {
		// cache miss
		seq = (seq + 1) % 100;
		off_t ret = global_rand_permute_workload::next_offset();
		if (cached_pages.size() >= (size_t) num_pages)
			cached_pages.pop_front();
		cached_pages.push_back(ret);
		return ret;
	}
}

off_t *load_java_dump(const std::string &file, long &num_offsets)
{
	int fd = open(file.c_str(), O_RDONLY);
	if (fd < 0) {
		perror("open");
		exit(1);
	}
	printf("%s's fd is %d\n", file.c_str(), fd);

	/* get the file size */
	long file_size = get_file_size(file.c_str());
	assert(file_size % sizeof(off_t) == 0);

	off_t *offsets = (off_t *) malloc(file_size);
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
	num_offsets = file_size / sizeof(off_t);
	return offsets;
}

workload_t *load_file_workload(const std::string &file, long &num)
{
	int fd = open(file.c_str(), O_RDONLY);
	if (fd < 0) {
		perror("open");
		exit(1);
	}
	printf("%s's fd is %d\n", file.c_str(), fd);

	/* get the file size */
	long file_size = get_file_size(file.c_str());
	assert(file_size % sizeof(workload_t) == 0);

	workload_t *workloads = (workload_t *) malloc(file_size);
	/* read data of the file to a buffer */
	char *buf = (char *) workloads;
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
	num = file_size / sizeof(workload_t);
	return workloads;
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

template class thread_safe_FIFO_queue<workload_t>;

workload_chunk *balanced_workload::chunks;
thread_safe_FIFO_queue<off_t> *global_rand_permute_workload::permuted_offsets;
thread_safe_FIFO_queue<workload_t> *file_workload::workload_queue;
