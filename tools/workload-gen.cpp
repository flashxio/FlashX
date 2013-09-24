/**
 * This file is to read block IO trace format and convert it
 * to the workload structure.
 */

#include <stdio.h>

#include <tr1/unordered_map>

#include "../workload.h"

bool histogram = true;

/**
 * This class can convert a sequence of numbers into another sequence,
 * but the frequency of numbers should be preserved.
 */
class rand_rehasher
{
	long *rehash_map;
	long rehash_zero;
	long st;
	long range;
public:
	rand_rehasher(long st, long n) {
		srandom(time(NULL));
		this->range = n - st;
		this->st = st;
		rehash_map = new long[range];
		memset(rehash_map, 0, sizeof(rehash_map[0]) * range);
		rehash_zero = random() % range + st;
	}

	~rand_rehasher() {
		delete rehash_map;
	}

	long rehash(long v) {
		assert(v >= st && v < st + range);
		if (v == 0)
			return rehash_zero;
		else if (rehash_map[v - st] == 0) {
			rehash_map[v - st] = random() % range + st;
		}
		return rehash_map[v - st];
	}
};

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
	std::tr1::unordered_map<off_t, int> page_map;

	rand_rehasher rehasher(0, 100 * 1024 * 1024);

	FILE *in = fopen(argv[1], "r");
	if (in == NULL) {
		perror("fopen");
		assert(in);
	}
	FILE *out = fopen(argv[2], "w");
	if (out == NULL) {
		perror("fopen");
		assert(out);
	}
	assert(out);
	char buf[1024];
	char orig[1024];
	int num_lines = 0;
	while (fgets(buf, sizeof(buf), in)) {
		memcpy(orig, buf, sizeof(buf));
		num_lines++;
		workload_t workload;
		char *str = strchr(buf, ',');
		if (str == NULL) {
			printf("skip %s\n", orig);
			continue;
		}
		char *off_str = str + 1;
		str = strchr(off_str, ',');
		if (str == NULL) {
			printf("skip %s\n", orig);
			continue;
		}
		*str = 0;
		char *size_str = str + 1;
		str = strchr(size_str, ',');
		if (str == NULL) {
			printf("skip %s\n", orig);
			continue;
		}
		*str = 0;
		char *flag_str = str + 1;
		workload.off = strtol(off_str, NULL, 10);
		assert(workload.off >= 0);
//		workload.off = rehasher.rehash(workload.off / 4096) * 4096;
		if (min_off > workload.off)
			min_off = workload.off;
		if (max_off < workload.off)
			max_off = workload.off;
		workload.size = atoi(size_str);
		if (workload.size <= 0) {
			printf("skip %s\n", orig);
			continue;
		}
		if (min_size > workload.size)
			min_size = workload.size;
		if (max_size < workload.size)
			max_size = workload.size;
		workload.read = flag_str[0] == 'r' || flag_str[0] == 'R';
		off_t begin_page_off = ROUND_PAGE(workload.off);
		off_t end_page_off = ROUNDUP_PAGE(workload.off + workload.size);
		for (off_t off = begin_page_off; off < end_page_off; off += PAGE_SIZE) {
			std::tr1::unordered_map<off_t, int>::iterator it = page_map.find(off);
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
	printf("read %d lines\n", num_lines);
	printf("There are %ld pages accessed\n", page_map.size());
	printf("The min offset is %ld, the max offset is %ld\n", min_off, max_off);
	printf("The min request size is %ld, the max request size is %ld\n", min_size, max_size);
	printf("There are %d reads and %d writes\n", num_reads, num_writes);
	printf("There are %d breakdown reads (in pages) and %d breakdown writes (in pages)\n",
			num_breakdown_reads, num_breakdown_writes);
	printf("It reads %ld bytes and writes %ld bytes\n", num_read_bytes, num_write_bytes);
}
