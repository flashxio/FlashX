#ifndef __ASSOCIATIVE_CACHE_H__
#define __ASSOCIATIVE_CACHE_H__

#include <pthread.h>
#include <math.h>

#include <vector>

#include "memory_manager.h"
#include "cache.h"
#include "concurrency.h"
#include "container.h"

#define CELL_SIZE 8

/* 36 shadow pages makes exactly 4 cache lines. */
#define NUM_SHADOW_PAGES 36

#ifdef STATISTICS
volatile extern int avail_cells;
volatile extern int num_wait_unused;
volatile extern int lock_contentions;
#endif

const int CACHE_LINE = 128;

/**
 * This data structure is to implement LRU.
 */
template<class T, int BUF_SIZE>
class page_cell
{
	unsigned int idx;		// to the point where we can evict a page in the buffer
	T buf[BUF_SIZE];			// a circular buffer to keep pages.

public:
	/*
	 * @size: the size of the page buffer
	 * @page_buf: the offset of the page array in the global page cache.
	 */
	page_cell() {
		idx = 0;
	}

	void set_pages(char *pages[]) {
		for (int i = 0; i < BUF_SIZE; i++) {
			buf[i] = T(-1, pages[i]);
		}
		idx = 0;
	}

	/**
	 * return an empty page.
	 * I expected the page will be filled with data,
	 * so I change the begin and end index of the circular buffer.
	 */
	T *get_empty_page() {
		/* TODO I ignore the case of integer overflow */
		T *ret = &buf[idx % BUF_SIZE];
		idx++;
		return ret;
	}

	T *get_page(int i) {
		if (i >= BUF_SIZE)
			return NULL;
		return &buf[i];
	}

	int get_idx(T *page) {
		int idx = page - buf;
		assert (idx >= 0 && idx < BUF_SIZE);
		return idx;
	}

	void scale_down_hits() {
		for (int i = 0; i < BUF_SIZE; i++) {
			buf[i].set_hits(buf[i].get_hits() / 2);
		}
	}
};

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

class associative_cache;

class hash_cell
{
	// for testing
	long hash;

	pthread_spinlock_t _lock;
	page_cell<thread_safe_page, CELL_SIZE> buf;
	bool overflow;
	associative_cache *table;
#ifdef USE_LRU
	std::vector<int> pos_vec;
#endif
#ifdef USE_SHADOW_PAGE
	clock_shadow_cell shadow;
#endif

	thread_safe_page *get_empty_page();

public:
	hash_cell() {
		overflow = false;
		table = NULL;
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

	hash_cell(associative_cache *cache, long hash);

	~hash_cell() {
		pthread_spin_destroy(&_lock);
	}

	void *operator new[](size_t size) {
		printf("allocate %ld bytes\n", size);
		void *addr = memalign(CACHE_LINE, size + CACHE_LINE);
		// TODO 8 might be architecture specific. It's 8 for 64-bit machines.
		return (void *) ((long) addr + CACHE_LINE - 8);
	}

	void operator delete[](void *p) {
		free((void *) ((long) p - (CACHE_LINE - 8)));
	}

	bool is_overflow() {
		return overflow;
	}

	page *search(off_t off, off_t &old_off);
	/* 
	 * this is to rehash the pages in the current cell
	 * to the cell in the parameter.
	 */
	void rehash(hash_cell *cell);

	void print_cell() {
		for (int i = 0; i < CELL_SIZE; i++)
			printf("%lx\t", buf.get_page(i)->get_offset());
		printf("\n");
	}
};

const long init_cache_size = 128 * 1024 * 1024;

class associative_cache: public page_cache
{
	enum {
		TABLE_EXPANDING,
	};
	/* 
	 * this table contains cell arrays.
	 * each array contains N cells;
	 * TODO we might need to use map to improve performance.
	 */
	std::vector<hash_cell*> cells_table;
	atomic_integer ncells;
	
	seq_lock table_lock;
	atomic_flags flags;
	/* the initial number of cells in the table. */
	int init_ncells;

	memory_manager *manager;

	/* used for linear hashing */
	int level;
	int split;

public:
	associative_cache(memory_manager *manager);

	~associative_cache() {
		for (unsigned int i = 0; i < cells_table.size(); i++)
			if (cells_table[i])
				delete [] cells_table[i];
		manager->unregister_cache(this);
	}

	/* the hash function used for the current level. */
	int hash(off_t offset) {
		return offset / PAGE_SIZE % (init_ncells * (long) pow(2, level));
	}

	/* the hash function used for the next level. */
	int hash1(off_t offset) {
		return offset / PAGE_SIZE % (init_ncells * (long) pow(2, level + 1));
	}

	int hash1_locked(off_t offset) {
		unsigned long count;
		int ret;
		do {
			table_lock.read_lock(count);
			ret = offset / PAGE_SIZE % (init_ncells * (long) pow(2, level + 1));
		} while (!table_lock.read_unlock(count));
		return ret;
	}

	page *search(off_t offset, off_t &old_off);

	bool expand(hash_cell *cell);

	memory_manager *get_manager() {
		return manager;
	}

	void print_cell(off_t off) {
		get_cell(off)->print_cell();
	}

	long size() {
		return ((long) ncells.get())
			* init_ncells * CELL_SIZE * PAGE_SIZE;
	}

	bool shrink(int npages, char *pages[]) {
		// TODO shrink the cache
		return false;
	}

	hash_cell *get_cell(unsigned int global_idx) {
		unsigned int cells_idx = global_idx / init_ncells;
		int idx = global_idx % init_ncells;
		assert(cells_idx < cells_table.size());
		hash_cell *cells = cells_table[cells_idx];
		if (cells)
			return &cells[idx];
		else
			return NULL;
	}

	hash_cell *get_cell_offset(off_t offset) {
		int global_idx;
		unsigned long count;
		hash_cell *cell = NULL;
		do {
			table_lock.read_lock(count);
			global_idx = hash(offset);
			if (global_idx < split)
				global_idx = hash1(offset);
			cell = get_cell(global_idx);
		} while (!table_lock.read_unlock(count));
		assert(cell);
		return cell;
	}
};

#endif
