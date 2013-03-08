#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <vector>
#include <algorithm>

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
	std::vector<int> nums(file_size / 8, 0);

	char buf[8];
	long *src = (long *) buf;
	int count = 0;
	while(true) {
		ssize_t ret = read(fd, buf, sizeof(buf));
		if (ret != 8)
			break;
		count++;
		long n = java_dump_workload::swap_bytesl(*src);
		nums[n / 4096]++;
		if (n > max_pos)
			max_pos = n;
	}

	std::sort(nums.begin(), nums.end());

	long naccesses = 0;
	long accessed_pages = 0;
	for (int i = nums.size() - 1; i >= 0; i--) {
		if (nums[i] == 0)
			break;
		naccesses += nums[i];
		if (i % 100 == 0)
			printf("%ld\t%ld\n", nums.size() - i, naccesses);
		accessed_pages++;
	}
	printf("there are %ld numbers accessed\n", naccesses);
	printf("there are %ld accessed pages\n", accessed_pages);
	printf("the max location is %ld\n", max_pos);
}
