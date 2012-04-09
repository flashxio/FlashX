#ifndef __ASSOCIATIVE_CACHE_H__
#define __ASSOCIATIVE_CACHE_H__

#include <pthread.h>
#include <math.h>

#include <vector>

#include "memory_manager.h"
#include "cache.h"

#define CELL_SIZE 8

#define USE_SHADOW_PAGE

/* 36 shadow pages makes exactly 4 cache lines. */
#define NUM_SHADOW_PAGES 36

#ifdef STATISTICS
volatile extern int avail_cells;
volatile extern int num_wait_unused;
volatile extern int lock_contentions;
#endif
extern int end_evicts;
extern int middle_evicts;

const int CACHE_LINE = 128;

class expand_exception
{
};

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

template<class T, int SIZE>
class generic_queue
{
	unsigned short start;
	unsigned short num;
	/* the size of the buffer is specified by SIZE. */
	T buf[SIZE];
public:
	generic_queue() {
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
		end_evicts++;
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
	generic_queue<shadow_page, NUM_SHADOW_PAGES> queue;
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
	generic_queue<shadow_page, NUM_SHADOW_PAGES> queue;
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

class atomic_flags
{
	volatile int flags;
public:
	atomic_flags() {
		flags = 0;
	}

	void set_flags(int flag) {
		__sync_fetch_and_or(&flags, 0x1 << flag);
	}

	void clear_flags(int flag) {
		__sync_fetch_and_and(&flags, ~(0x1 << flag));
	}

	bool test_flags(int flag) {
		return flags & (0x1 << flag);
	}
};

class seq_lock
{
	volatile unsigned long count;
	volatile int lock;
public:
	seq_lock() {
		count = 0;
		lock = 0;
	}

	void read_lock(unsigned long &count) {
		count = this->count;
	}

	bool read_unlock(unsigned long count) {
		return this->count == count;
	}

	void write_lock() {
		while(__sync_fetch_and_or(&lock, 1)) {}
		__sync_fetch_and_add(&count, 1);
	}

	void write_unlock() {
		__sync_fetch_and_and(&lock, 0);
	}
};

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
	volatile std::vector<hash_cell*> cells_table;
	// TODO it might be better to use seq_lock
	// because the table doesn't change much.
	seq_lock table_lock;
	atomic_flags flags;
	/* the initial number of cells in the table. */
	int init_ncells;

	memory_manager *manager;

	/* used for linear hashing */
	volatile int level;
	int split;

public:
	associative_cache(memory_manager *manager) {
		printf("associative cache is used\n");
		level = 0;
		split = 0;
		this->manager = manager;
		manager->register_cache(this);
		int npages = init_cache_size / PAGE_SIZE;
		assert(init_cache_size >= CELL_SIZE * PAGE_SIZE);
		init_ncells = npages / CELL_SIZE;
		hash_cell *cells = new hash_cell[init_ncells];
		printf("%d cells: %p\n", init_ncells, cells);
		for (int i = 0; i < init_ncells; i++)
			cells[i] = hash_cell(this, i);
		cells_table.push_back(cells);
	}

	~associative_cache() {
		for (unsigned int i = 0; i < cells_table.size(); i++)
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
		return ((long) cells_table.size())
			* init_ncells * CELL_SIZE * PAGE_SIZE;
	}

	bool shrink(int npages, char *pages[]) {
		// TODO shrink the cache
		return false;
	}

	hash_cell *get_cell(unsigned int global_idx) {
		unsigned int cells_idx = global_idx / init_ncells;
		int idx = global_idx % init_ncells;
		// TODO this check isn't enough after I use seq_lock.
		assert(cells_idx < cells_table.size());
		hash_cell *cells = cells_table[cells_idx];
		return &cells[idx];
	}

	hash_cell *get_cell_offset(off_t offset) {
		int global_idx;
		unsigned long count;
		hash_cell *cell;
		do {
			table_lock.read_lock(count);
			global_idx = hash(offset);
			if (global_idx < split)
				global_idx = hash1(offset);
			cell = get_cell(global_idx);
		} while (!table_lock.read_unlock(count));
		return cell;
	}
};

#endif
