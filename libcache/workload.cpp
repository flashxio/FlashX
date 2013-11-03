#include "workload.h"
#include "common.h"

int workload_gen::default_entry_size = PAGE_SIZE;
int workload_gen::default_access_method = -1;

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

template class thread_safe_FIFO_queue<workload_t>;

thread_safe_FIFO_queue<off_t> *global_rand_permute_workload::permuted_offsets;
thread_safe_FIFO_queue<workload_t> *file_workload::workload_queue;
