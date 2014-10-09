/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sched.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <fcntl.h>
#include <sys/mman.h>

#include <sstream>

#include "common.h"
#include "messaging.h"

extern "C" {

int isnumeric(char *str)
{
	int len = strlen(str);
	for (int i = 0; i < len; i++) {
		if (!isdigit(str[i]))
			return 0;
	}
	return 1;
}

static const int huge_page_size = 2 * 1024 * 1024;
static bool use_huge_page = false;

void set_use_huge_page(bool v)
{
	printf("use huge page: %d\n", v);
	use_huge_page = v;
}

void *malloc_large(size_t size)
{
#define PROTECTION (PROT_READ | PROT_WRITE)
	/* Only ia64 requires this */
#ifdef __ia64__
#define ADDR (void *)(0x8000000000000000UL)
#define FLAGS (MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | MAP_FIXED)
#else
#define ADDR (void *)(0x0UL)
#define FLAGS (MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB)
#endif
	if (use_huge_page) {
		size = ROUNDUP(size, huge_page_size);
		void *addr = mmap(ADDR, size, PROTECTION, FLAGS, 0, 0);
		if (addr == MAP_FAILED) {
			perror("mmap");
			return NULL;
		}
		return addr;
	}
	else {
		return numa_alloc_local(size);
	}
}

void free_large(void *addr, size_t size)
{
	if (use_huge_page) {
		size = ROUNDUP(size, huge_page_size);
		munmap(addr, size);
	}
	else {
		numa_free(addr, size);
	}
}

}

bool check_read_content(char *buf, int size, off_t off, int file_id)
{
	// I assume the space in the buffer is larger than 8 bytes.
	off_t aligned_off = off & (~(sizeof(off_t) - 1));
	long data[2];
	data[0] = aligned_off / sizeof(off_t) + file_id;
	data[1] = aligned_off / sizeof(off_t) + 1 + file_id;
	long expected = 0;
	int copy_size = size < (int) sizeof(off_t) ? size : (int) sizeof(off_t);
	memcpy(&expected, ((char *) data) + (off - aligned_off), copy_size);
	long read_value = 0;
	memcpy(&read_value, buf, copy_size);
	if(read_value != expected)
		printf("off: %ld, size: %d, read: %ld, expect: %ld\n",
				off, size, read_value, expected);
	return read_value == expected;
}

void create_write_data(char *buf, int size, off_t off, int file_id)
{
	off_t aligned_start = off & (~(sizeof(off_t) - 1));
	off_t aligned_end = (off + size) & (~(sizeof(off_t) - 1));
	long start_data = aligned_start / sizeof(off_t) + file_id;
	long end_data = aligned_end / sizeof(off_t) + file_id;

	/* If all data is in one 8-byte word. */
	if (aligned_start == aligned_end) {
		memcpy(buf, ((char *) &start_data) + (off - aligned_start), size);
		return;
	}

	int first_size =  (int)(sizeof(off_t) - (off - aligned_start));
	int last_size = (int) (off + size - aligned_end);

	if (first_size == sizeof(off_t))
		first_size = 0;
	if (first_size)
		memcpy(buf, ((char *) &start_data) + (off - aligned_start),
				first_size);
	for (int i = first_size; i < aligned_end - off; i += sizeof(off_t)) {
		*((long *) (buf + i)) = (off + i) / sizeof(off_t) + file_id;
	}
	if (aligned_end > aligned_start
			|| (aligned_end == aligned_start && first_size == 0)) {
		if (last_size)
			memcpy(buf + (aligned_end - off), (char *) &end_data, last_size);
	}

	check_read_content(buf, size, off, file_id);
}

bool align_check(size_t alignment)
{
	assert(alignment >= 0);
	if (alignment == 0 || alignment == 1) {
		return false;
	}
	bool aligned = true;
	while (alignment > 1) {
		// If it's not a power of 2.
		if (alignment % 2) {
			aligned = false;
			break;
		}
		alignment /= 2;
	}
	return aligned;
}

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

int split_string(const std::string &str, char delim,
		std::vector<std::string> &strs)
{
	std::stringstream ss(str);
	std::string item;
	while (std::getline(ss, item, delim)) {
		strs.push_back(item);
	}
	return strs.size();
}
