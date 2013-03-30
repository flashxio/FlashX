#include <stdlib.h>
#include <stdio.h>

#include <string>

long str2size(std::string str)
{
	int len = str.length();
	long multiply = 1;
	if (str[len - 1] == 'M' || str[len - 1] == 'm') {
		multiply *= 1024 * 1024;
		str[len - 1] = 0;
	}
	else if (str[len - 1] == 'K' || str[len - 1] == 'k') {
		multiply *= 1024;
		str[len - 1] = 0;
	}
	else if (str[len - 1] == 'G' || str[len - 1] == 'g') {
		multiply *= 1024 * 1024 * 1024;
		str[len - 1] = 0;
	}
	return atol(str.c_str()) * multiply;
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "memory-fill size\n");
		return -1;
	}

	std::string size_str = argv[1];
	long size = str2size(size_str);
	long one_alloc = 1024 * 1024 * 1024;
	long remaining_size = size;
	long allocated_size = 0;
	for (; allocated_size < size;) {
		long alloc = one_alloc;
		if (remaining_size < alloc)
			alloc = remaining_size;
		char *buf = (char *) malloc(alloc);
		for (int off = 0; off < alloc; off += 4096)
			buf[off] = 0;
		remaining_size -= alloc;
		allocated_size += alloc;
	}
	unsigned int sleep_time = 1 << 31;
	printf("allocate %ld bytes\n", allocated_size);
	printf("sleep %d seconds\n", sleep_time);
	sleep(sleep_time);
}
