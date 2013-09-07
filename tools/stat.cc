#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <limits.h>

#include <vector>
#include <algorithm>
#include <tr1/unordered_map>

#define PAGE_SIZE 4096
#define SCALE_FACTOR 100

#include "workload.h"

bool histogram = true;

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "stat file\n");
		exit(1);
	}

	char *name = argv[1];
	int fd = open(name, O_RDONLY);
	if (fd < 0) {
		perror("open");
		exit(1);
	}
	long max_pos = 0;

	/* get the file size */
	struct stat stats;
	if (fstat(fd, &stats) < 0) {
		perror("fstat");
		exit(1);
	}
	long file_size = stats.st_size;

	/* the numbers of accesses of each page */
	int num_accesses = (int) (file_size / sizeof(workload_t));
	workload_t *workloads = new workload_t[num_accesses];
	ssize_t ret = read(fd, (void *) workloads, file_size);

	std::tr1::unordered_map<off_t, int> page_map;
	int num_reads = 0;
	int num_writes = 0;
	int max_size = 0;
	int min_size = INT_MAX;
	long tot_size = 0;
	for (int i = 0; i < num_accesses; i++) {
		off_t page_off = ROUND_PAGE(workloads[i].off);
		if (max_size < workloads[i].size)
			max_size = workloads[i].size;
		if (min_size > workloads[i].size)
			min_size = workloads[i].size;
		tot_size += workloads[i].size;
		off_t last_page = ROUNDUP_PAGE(workloads[i].off + workloads[i].size);
		int num_pages = (last_page - page_off) / PAGE_SIZE;
		for (; page_off < last_page; page_off += PAGE_SIZE) {
			std::tr1::unordered_map<off_t, int>::iterator it = page_map.find(page_off);
			if (it == page_map.end()) {
				page_map.insert(std::pair<off_t, int>(page_off, 1));
			}
			else
				it->second++;
		}
		if (workloads[i].read)
			num_reads += num_pages;
		else
			num_writes += num_pages;
	}

	if (histogram) {
		std::map<int, int> count_nums;
		for (std::tr1::unordered_map<off_t, int>::iterator it = page_map.begin();
				it != page_map.end(); it++) {
			int num_accesses = it->second;
			std::map<int, int>::iterator it1 = count_nums.find(num_accesses);
			if (it1 != count_nums.end())
				it1->second++;
			else
				count_nums.insert(std::pair<int, int>(num_accesses, 1));
		}
		for (std::map<int, int>::const_iterator it = count_nums.begin();
				it != count_nums.end(); it++) {
			printf("%d pages get %d hits\n", it->second, it->first);
		}
	}
	printf("min request size: %d, max req size: %d, avg size: %ld\n",
			min_size, max_size, tot_size / num_accesses);
	printf("there are %d accesses\n", num_accesses);
	printf("there are %d reads\n", num_reads);
	printf("there are %d writes\n", num_writes);
	printf("there are %ld accessed pages\n", page_map.size());
}
