/**
 * This file is to read block IO trace format and convert it
 * to the workload structure.
 */

#include <stdio.h>

#include <map>

#include "../workload.h"

bool histogram = true;

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "workload-gen input_file output_file\n");
		return 1;
	}

	off_t min_off = (1L << 62) - 1;
	off_t max_off = 0;
	size_t min_size = (1L << 62) - 1;
	size_t max_size = 0;
	int num_reads = 0;
	int num_writes = 0;
	int num_breakdown_reads = 0;
	int num_breakdown_writes = 0;
	size_t num_read_bytes = 0;
	size_t num_write_bytes = 0;
	// calculate which pages are accessed and how many times a page is accessed.
	std::map<off_t, int> page_map;

	FILE *in = fopen(argv[1], "r");
	FILE *out = fopen(argv[2], "w");
	char buf[1024];
	while (fgets(buf, sizeof(buf), in)) {
		workload_t workload;
		char *str = strchr(buf, ',');
		if (str == NULL)
			break;
		char *off_str = str + 1;
		str = strchr(off_str, ',');
		if (str == NULL)
			break;
		*str = 0;
		char *size_str = str + 1;
		str = strchr(size_str, ',');
		if (str == NULL)
			break;
		*str = 0;
		char *flag_str = str + 1;
		workload.off = strtol(off_str, NULL, 10);
		if (min_off > workload.off)
			min_off = workload.off;
		if (max_off < workload.off)
			max_off = workload.off;
		workload.size = atoi(size_str);
		if (min_size > workload.size)
			min_size = workload.size;
		if (max_size < workload.size)
			max_size = workload.size;
		workload.read = flag_str[0] == 'r' || flag_str[0] == 'R';
		off_t begin_page_off = ROUND_PAGE(workload.off);
		off_t end_page_off = ROUNDUP_PAGE(workload.off + workload.size);
		for (off_t off = begin_page_off; off < end_page_off; off += PAGE_SIZE) {
			std::map<off_t, int>::iterator it = page_map.find(off);
			if (it == page_map.end()) {
				page_map.insert(std::pair<off_t, int>(off, 1));
			}
			else
				it->second++;
		}
		if (workload.read) {
			num_reads++;
			num_breakdown_reads += (end_page_off - begin_page_off) / PAGE_SIZE;
			num_read_bytes += workload.size;
		}
		else {
			num_writes++;
			num_breakdown_writes += (end_page_off - begin_page_off) / PAGE_SIZE;
			num_write_bytes += workload.size;
		}
		size_t ret = fwrite(&workload, sizeof(workload), 1, out);
		assert(ret == 1);
	}
	fclose(in);
	fclose(out);

	if (histogram) {
		for (std::map<off_t, int>::iterator it = page_map.begin(); it != page_map.end(); it++) {
			printf("page %ld: %d\n", it->first, it->second);
		}
	}
	printf("There are %ld pages accessed\n", page_map.size());
	printf("The min offset is %ld, the max offset is %ld\n", min_off, max_off);
	printf("The min request size is %ld, the max request size is %ld\n", min_size, max_size);
	printf("There are %d reads and %d writes\n", num_reads, num_writes);
	printf("There are %d breakdown reads (in pages) and %d breakdown writes (in pages)\n",
			num_breakdown_reads, num_breakdown_writes);
	printf("It reads %ld bytes and writes %ld bytes\n", num_read_bytes, num_write_bytes);
}
