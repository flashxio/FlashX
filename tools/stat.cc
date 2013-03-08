#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <vector>
#include <algorithm>
#include <tr1/unordered_set>

#define PAGE_SIZE 4096
#define SCALE_FACTOR 100

#include "../workload.h"

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

	std::tr1::unordered_set<off_t> pages;
	int num_reads = 0;
	int num_writes = 0;
	for (int i = 0; i < num_accesses; i++) {
		off_t page_off = ROUND_PAGE(workloads[i].off);
		off_t last_page = ROUNDUP_PAGE(workloads[i].off + workloads[i].size);
		int num_pages = (last_page - page_off) / PAGE_SIZE;
		for (; page_off < last_page; page_off += PAGE_SIZE)
			pages.insert(page_off);
		if (workloads[i].read)
			num_reads += num_pages;
		else
			num_writes += num_pages;
	}

	printf("there are %d reads\n", num_reads);
	printf("there are %d writes\n", num_writes);
	printf("there are %ld accessed pages\n", pages.size());
}
