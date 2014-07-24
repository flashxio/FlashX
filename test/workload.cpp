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

#include "workload.h"
#include "common.h"
#include "native_file.h"

int workload_gen::default_entry_size = PAGE_SIZE;
int workload_gen::default_access_method = -1;

void permute_offsets(int num, int repeats, int stride, off_t start,
		off_t offsets[])
{
	int idx = 0;
	for (int k = 0; k < repeats; k++) {
		for (int i = 0; i < num; i++) {
			offsets[idx++] = ((off_t) i) * stride + start * stride;
		}
	}
	int tot_length = idx;

	for (int i = tot_length - 1; i >= 1; i--) {
		int j = random() % tot_length;
		off_t tmp = offsets[j];
		offsets[j] = offsets[i];
		offsets[i] = tmp;
	}
}

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
	native_file f(file);
	long file_size = f.get_size();
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
	native_file f(file);
	long file_size = f.get_size();
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
thread_safe_FIFO_queue<fifo_queue<workload_t> *> *file_workload::workload_queue;

thread_safe_FIFO_queue<workload_pack> *dynamic_rand_workload::workload_queue;
int dynamic_rand_workload::pack_size = 1000;
