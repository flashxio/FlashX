#ifndef __SHADOW_CELL_H__
#define __SHADOW_CELL_H__

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
