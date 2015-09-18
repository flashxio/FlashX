#ifndef __SHADOW_CELL_H__
#define __SHADOW_CELL_H__

/*
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

#include "cache.h"

/* 36 shadow pages makes exactly 4 cache lines. */
#define NUM_SHADOW_PAGES 36

class shadow_page
{
	int offset;
	unsigned char hits;
	char flags;
public:
	shadow_page() {
		offset = -1;
		hits = 0;
		flags = 0;
	}
	shadow_page(page &pg) {
		offset = pg.get_offset() >> LOG_PAGE_SIZE;
		hits = pg.get_hits();
		flags = 0;
	}

	void set_referenced(bool referenced) {
		if (referenced)
			flags |= 0x1 << REFERENCED_BIT;
		else
			flags &= ~(0x1 << REFERENCED_BIT);
	}
	bool referenced() {
		return flags & (0x1 << REFERENCED_BIT);
	}

	off_t get_offset() const {
		return ((off_t) offset) << LOG_PAGE_SIZE;
	}

	int get_hits() {
		return hits;
	}

	void set_hits(int hits) {
		assert (hits <= 0xff);
		this->hits = hits;
	}

	bool is_valid() {
		return offset != -1;
	}
};

#ifdef USE_SHADOW_PAGE

class shadow_cell
{
public:
	virtual void add(shadow_page pg) = 0;
	virtual shadow_page search(off_t off) = 0;
	virtual void scale_down_hits() = 0;
};

/**
 * The elements in the queue stored in the same piece of memory
 * as the queue metadata. The size of the queue is defined 
 * during the compile time.
 */
template<class T, int SIZE>
class embedded_queue
{
	unsigned short start;
	unsigned short num;
	/* the size of the buffer is specified by SIZE. */
	T buf[SIZE];
public:
	embedded_queue() {
		assert(SIZE < 0xffff);
		start = 0;
		num = 0;
	}

	void push_back(T v) {
		assert(num < SIZE);
		buf[(start + num) % SIZE] = v;
		num++;
	}

	void pop_front() {
		assert(num > 0);
		start = (start + 1) % SIZE;
		num--;
	}

	void remove(int idx);

	bool is_empty() {
		return num == 0;
	}

	bool is_full() {
		return num == SIZE;
	}

	int size() {
		return num;
	}

	T &back() {
		assert(num > 0);
		return buf[(start + num - 1) % SIZE];
	}

	T &front() {
		assert(num > 0);
		return buf[start];
	}

	T &get(int idx) {
		assert(num > 0);
		return buf[(start + idx) % SIZE];
	}

	void set(T &v, int idx) {
		buf[(start + idx) % SIZE] = v;
	}

	void print_state() {
		printf("start: %d, num: %d\n", start, num);
		for (int i = 0; i < this->size(); i++)
			printf("%ld\t", this->get(i));
		printf("\n");
	}
};

class clock_shadow_cell: public shadow_cell
{
	int last_idx;
	// TODO adjust to make it fit in cache lines.
	embedded_queue<shadow_page, NUM_SHADOW_PAGES> queue;
public:
	clock_shadow_cell() {
		last_idx = 0;
	}

	void add(shadow_page pg);

	shadow_page search(off_t off);

	void scale_down_hits();
};

class LRU_shadow_cell: public shadow_cell
{
	embedded_queue<shadow_page, NUM_SHADOW_PAGES> queue;
public:
	LRU_shadow_cell() {
	}

	/*
	 * add a page to the cell, which is evicted from hash_cell.
	 * the only thing we need to record is the number of hits
	 * of this page.
	 */
	void add(shadow_page pg) {
		if (queue.is_full())
			queue.pop_front();
		queue.push_back(pg);
	}

	shadow_page search(off_t off);

	void scale_down_hits();
};

#endif

#endif
