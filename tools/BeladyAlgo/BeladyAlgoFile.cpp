#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <google/profiler.h>

#include <vector>

#include "BeladyAlgo.h"
#include "common.h"
#include "workload.h"

static std::string prof_file = "BeladyAlgo.prof";

void int_handler(int sig_num)
{
	if (!prof_file.empty())
		ProfilerStop();
	exit(0);
}

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "%s file num_pages\n", argv[0]);
		return -1;
	}

	int cache_size = atoi(argv[2]);		// the number of pages
	int fd = open(argv[1], O_RDONLY);
	if (fd < 0) {
		perror("open");
		return -1;
	}

	// Get the file size
	struct stat stat_buf;
	fstat(fd, &stat_buf);
	size_t file_size = stat_buf.st_size;
	assert(file_size % sizeof(workload_t) == 0);

	// Load accesses
	long start = get_curr_ms();
	int tot_num_accesses = file_size / sizeof(workload_t);
	int num_accesses = 0;
	std::vector<int> offs;
	workload_t workload;
	// Find the location of the access sequence where we start
	// to count cache hits.
	// Basically, we use the first half to warm up the cache
	// and the second half to measure the cache hit rate.
	int half_access_pointer = 0;
	while (read(fd, &workload, sizeof(workload)) > 0) {
		num_accesses++;
		off_t off = workload.off;
		off_t end = off + workload.size;
		while (off < end) {
			if (num_accesses < tot_num_accesses / 2) {
				half_access_pointer++;
			}
			int tmp = (int) (off / 4096);
			off = ((long) tmp) * 4096 + 4096;
			offs.push_back(tmp);
		}
	}
	long end = get_curr_ms();
	printf("loading all access data takes %ld ms\n", end - start);
	close(fd);

	signal(SIGINT, int_handler);

	belady_algo algo(cache_size);
	indexed_offset_scanner scanner(offs.data(), (int) offs.size());
	if (!prof_file.empty())
		ProfilerStart(prof_file.c_str());
	int nhits = algo.access(scanner, half_access_pointer);
	if (!prof_file.empty())
		ProfilerStop();
	printf("There are %d hits among %ld accesses\n", nhits, offs.size() - half_access_pointer);
}
