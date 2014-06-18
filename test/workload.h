#ifndef __WORKLOAD_H__
#define __WORKLOAD_H__

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

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <malloc.h>
#include <stdlib.h>

#include <string>
#include <deque>

#include "container.h"
#include "cache.h"

#define CHUNK_SLOTS 1024

const int WORKLOAD_BUF_SIZE = 1024 * 32;

typedef struct workload_type
{
	off_t off;
	int size: 31;
	int read: 1;
} workload_t;

class workload_pack
{
	workload_t *queue_ptr;
	int curr;
	int size;
public:
	workload_pack() {
		queue_ptr = NULL;
		curr = 0;
		size = 0;
	}

	workload_pack(workload_t *workloads, int num) {
		queue_ptr = new workload_t[num];
		memcpy(queue_ptr, workloads, sizeof(*queue_ptr) * num);
		curr = 0;
		size = num;
	}

	bool is_empty() const {
		return curr == size;
	}

	workload_t pop_front() {
		return queue_ptr[curr++];
	}

	int get_num_remaining() const {
		return size - curr;
	}
};

class workload_gen
{
	workload_t access;
	static int default_entry_size;
	static int default_access_method;
public:
	static void set_default_entry_size(int entry_size) {
		default_entry_size = entry_size;
	}
	static int get_default_entry_size() {
		return default_entry_size;
	}
	static void set_default_access_method(int access_method) {
		default_access_method = access_method;
	}
	static int get_default_access_method() {
		return default_access_method;
	}

	workload_gen() {
		memset(&access, 0, sizeof(access));
	}

	/**
	 * This is a wrapper to the original interface `next_offset'
	 * so more info of an access is provided.
	 */
	virtual const workload_t &next() {
		assert(default_access_method >= 0);
		access.off = next_offset();
		access.size = default_entry_size;
		access.read = default_access_method == READ;
		return access;
	}

	/**
	 * The enxt offset in bytes.
	 */
	virtual off_t next_offset() = 0;
	virtual bool has_next() = 0;
	virtual ~workload_gen() { }
	virtual void print_state() {
	}
};

class seq_workload: public workload_gen
{
	long start;
	long end;
	long cur;
	int entry_size;

public:
	seq_workload(long start, long end, int entry_size) {
		this->start = start;
		this->end = end;
		this->cur = start;
		this->entry_size = entry_size;
	}

	off_t next_offset() {
		off_t next = cur;
		cur++;
		return next * entry_size;
	}

	bool has_next() {
		return cur < end;
	}
};

class global_rand_permute_workload: public workload_gen
{
	static thread_safe_FIFO_queue<off_t> *permuted_offsets;
	int num_reads_in_100;
	int num_accesses;
	workload_t access;
	fifo_queue<off_t> local_buf;

public:
	/**
	 * @start: the index of the first entry.
	 * @end: the index of the last entry.
	 */
	global_rand_permute_workload(int stride, int length, int repeats,
			double read_ratio): local_buf(-1, WORKLOAD_BUF_SIZE) {
		if (permuted_offsets == NULL) {
			int tot_length = length * repeats;
			off_t *offsets = new off_t[tot_length];
			permute_offsets(length, repeats, stride, 0, offsets);
			permuted_offsets = thread_safe_FIFO_queue<off_t>::create(
					std::string("permuted_offsets"), -1, tot_length);
			int ret = permuted_offsets->add(offsets, tot_length);
			assert(ret == tot_length);
			delete offsets;
		}
		this->num_reads_in_100 = (int) (read_ratio * 100);
		this->num_accesses = 0;
	}

	virtual ~global_rand_permute_workload() {
		if (permuted_offsets) {
			delete permuted_offsets;
			permuted_offsets = NULL;
		}
	}

	bool is_initialized() const {
		return permuted_offsets != NULL;
	}

	off_t next_offset() {
		assert(!local_buf.is_empty());
		return local_buf.pop_front();
	}

	bool has_next() {
		if (local_buf.is_empty()) {
			off_t offs[WORKLOAD_BUF_SIZE];
			int num = permuted_offsets->fetch(offs, WORKLOAD_BUF_SIZE);
			local_buf.add(offs, num);
		}
		return !local_buf.is_empty();
	}

	virtual const workload_t &next() {
		off_t off = next_offset();
		access.off = off;
		access.size = workload_gen::get_default_entry_size();
		if (num_accesses < num_reads_in_100)
			access.read = 1;
		else
			access.read = 0;
		num_accesses++;
		if (num_accesses >= 100)
			num_accesses = 0;
		return access;
	}
};

/**
 * In this workload generator, the expected cache hit ratio
 * can be defined by the user.
 */
class cache_hit_defined_workload: public global_rand_permute_workload
{
	long num_pages;
	double cache_hit_ratio;
	std::deque<off_t> cached_pages;
	long seq;				// the sequence number of accesses
	long cache_hit_seq;		// the sequence number of cache hits
public:
	cache_hit_defined_workload(int stride, int length, long cache_size,
			double hit_ratio, double read_ratio): global_rand_permute_workload(
				stride, length, 1, read_ratio) {
		// only to access the most recent pages.
		this->num_pages = cache_size / PAGE_SIZE / 100;
		cache_hit_ratio = hit_ratio;
		seq = 0;
		cache_hit_seq = 0;
	}

	off_t next_offset();
};

off_t *load_java_dump(const std::string &file, long &num_offsets);
workload_t *load_file_workload(const std::string &file, long &num);

/**
 * This workload generator statistically divides workloads between threads.
 */
class divided_file_workload: public workload_gen
{
	fifo_queue<workload_t> local_buf;
	workload_t curr;
public:
	divided_file_workload(workload_t workloads[], long length, int thread_id,
			int nthreads): local_buf(-1, length / nthreads) {
		long part_size = length / nthreads;
		local_buf.add(workloads + thread_id * part_size, part_size);
	}

	const workload_t &next() {
		assert(!local_buf.is_empty());
		curr = local_buf.pop_front();
		if (get_default_access_method() >= 0)
			curr.read = get_default_access_method() == READ;
		return curr;
	}

	off_t next_offset() {
		assert(!local_buf.is_empty());
		workload_t access = local_buf.pop_front();
		return access.off;
	}

	bool has_next() {
		return !local_buf.is_empty();
	}
};

/**
 * This workload generator dynamically assigns workloads to threads.
 */
class file_workload: public workload_gen
{
	static thread_safe_FIFO_queue<fifo_queue<workload_t> *> *workload_queue;
	fifo_queue<workload_t> *local_buf;
	workload_t curr;
public:
	file_workload(workload_t workloads[], long length, int thread_id,
			int nthreads, int read_percent = -1) {
		if (workload_queue == NULL) {
			workload_queue = thread_safe_FIFO_queue<fifo_queue<workload_t> *>::create(
					"file_workload_queue", -1, length / WORKLOAD_BUF_SIZE + 1);
			if (read_percent >= 0) {
				for (long i = 0; i < length; i++) {
					if (random() % 100 < read_percent)
						workloads[i].read = 1;
					else
						workloads[i].read = 0;
				}
			}
			for (long i = 0; i < length; i += WORKLOAD_BUF_SIZE) {
				fifo_queue<workload_t> *q = new fifo_queue<workload_t>(-1,
						WORKLOAD_BUF_SIZE);
				int local_len = min(WORKLOAD_BUF_SIZE, length - i);
				q->add(workloads + i, local_len);
				workload_queue->push_back(q);
			}
			printf("There are %ld I/O accesses, and there are %d in the queue\n",
					length, workload_queue->get_num_entries() * WORKLOAD_BUF_SIZE);
		}
		local_buf = NULL;
	}

	virtual ~file_workload() {
		if (workload_queue) {
			while (!workload_queue->is_empty()) {
				fifo_queue<workload_t> *q = workload_queue->pop_front();
				delete q;
			}
			thread_safe_FIFO_queue<fifo_queue<workload_t> *>::destroy(workload_queue);
			workload_queue = NULL;
		}
	}

	const workload_t &next() {
		assert(local_buf && !local_buf->is_empty());
		curr = local_buf->pop_front();
		return curr;
	}

	off_t next_offset() {
		assert(local_buf && !local_buf->is_empty());
		workload_t access = local_buf->pop_front();
		return access.off;
	}

	bool has_next() {
		if (local_buf && local_buf->is_empty()) {
			delete local_buf;
			local_buf = NULL;
		}

		if (local_buf == NULL) {
			int num = workload_queue->fetch(&local_buf, 1);
			if (num == 0)
				local_buf = NULL;
		}
		return local_buf != NULL;
	}

	virtual void print_state() {
		printf("file workload has %d global works\n",
				workload_queue->get_num_entries());
	}
};

class dynamic_rand_workload: public workload_gen
{
	static thread_safe_FIFO_queue<workload_pack> *workload_queue;
	static int pack_size;
	workload_pack local_pack;
	workload_t curr;
public:
	static void create_workload(long length, int stride, int read_percent) {
		if (workload_queue == NULL) {
			workload_t *workloads = new workload_t[length];
			for (int i = 0; i < length; i++) {
				workloads[i].off = (random() % length) * stride;
				workloads[i].size = get_default_entry_size();
				workloads[i].read = random() % 100 < read_percent;
			}
			workload_queue = thread_safe_FIFO_queue<workload_pack>::create(
					"rand_workload_queue", -1,
					(length + pack_size - 1) / pack_size);

			long curr = 0;
			while (curr < length) {
				workload_pack pack(workloads + curr, min(pack_size,
							length - curr));
				curr += pack.get_num_remaining();
				int ret = workload_queue->add(&pack, 1);
				assert(ret == 1);
			}
			delete [] workloads;
		}
	}

	dynamic_rand_workload() {
		assert(workload_queue);
	}

	virtual ~dynamic_rand_workload() {
		if (workload_queue) {
			thread_safe_FIFO_queue<workload_pack>::destroy(workload_queue);
			workload_queue = NULL;
		}
	}

	const workload_t &next() {
		assert(!local_pack.is_empty());
		curr = local_pack.pop_front();
		return curr;
	}

	off_t next_offset() {
		assert(!local_pack.is_empty());
		workload_t acc = local_pack.pop_front();
		return acc.off;
	}

	bool has_next() {
		if (local_pack.is_empty()) {
			int ret = workload_queue->fetch(&local_pack, 1);
			if (ret == 1)
				assert(!local_pack.is_empty());
			return ret == 1;
		}
		else
			return true;
	}
};

class rand_workload: public workload_gen
{
	workload_t access;
	long start;
	long range;
	long num;
	long tot_accesses;
	off_t *offsets;
	bool *access_methods;
public:
	rand_workload(long start, long end, int stride, long tot_accesses,
			int read_percent) {
		this->start = start;
		this->range = end - start;
		this->tot_accesses = tot_accesses;
		num = 0;
		offsets = (off_t *) valloc(sizeof(*offsets) * tot_accesses);
		for (int i = 0; i < tot_accesses; i++) {
			offsets[i] = (start + random() % range) * stride;
		}
		access_methods = (bool *) valloc(sizeof(bool) * tot_accesses);
		for (int i = 0; i < tot_accesses; i++) {
			if (random() % 100 < read_percent)
				access_methods[i] = READ;
			else
				access_methods[i] = WRITE;
		}
	}

	virtual ~rand_workload() {
		free(offsets);
		free(access_methods);
	}

	off_t next_offset() {
		return offsets[num++];
	}

	bool has_next() {
		return num < tot_accesses;
	}

	virtual const workload_t &next() {
		access.off = offsets[num];
		access.size = get_default_entry_size();
		access.read = access_methods[num] == READ;
		num++;
		return access;
	}

	virtual void print_state() {
		printf("rand workload has %ld works left\n", tot_accesses - num);
	}
};

#endif
