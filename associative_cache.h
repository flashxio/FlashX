#ifndef __ASSOCIATIVE_CACHE_H__
#define __ASSOCIATIVE_CACHE_H__

#include <pthread.h>
#include <math.h>

#include <vector>

#include "memory_manager.h"
#include "cache.h"
#include "concurrency.h"
#include "container.h"
#include "parameters.h"
#include "shadow_cell.h"

#ifdef STATISTICS
volatile extern int avail_cells;
volatile extern int num_wait_unused;
volatile extern int lock_contentions;
#endif

const int CACHE_LINE = 128;

/**
 * This data structure is to contain page data structures in the hash cell.
 * It has space large enough for maximal CELL_SIZE, but only some of them
 * are used. The actual number of pages in the data structure varies,
 * and has a minimal limit.
 */
template<class T>
class page_cell
{
	// to the point where we can evict a page in the buffer
	unsigned char idx;
	unsigned char min_num_pages;
	unsigned char num_pages;
	// There are gaps in the `buf' array but we expose a virtual array without
	// gaps, so this mapping is to help remove the gaps in the physical array.
	// The number of elements in `maps' is `num_pages'.
	unsigned char maps[CELL_SIZE];

	T buf[CELL_SIZE];			// a circular buffer to keep pages.

public:
	/*
	 * @size: the size of the page buffer
	 * @page_buf: the offset of the page array in the global page cache.
	 */
	page_cell() {
		idx = 0;
		min_num_pages = (CELL_SIZE + 1) / 2;
		num_pages = 0;
		memset(maps, 0, sizeof(maps));
	}

	void set_pages(char *pages[], int num) {
		assert(num <= CELL_SIZE);
		for (int i = 0; i < num; i++) {
			buf[i] = T(-1, pages[i]);
		}
		idx = 0;
		num_pages = num;
		for (int i = 0; i < num; i++) {
			maps[i] = i;
		}
	}

	/**
	 * The two methods might be expensive because it reconstruct
	 * the entire mapping array.
	 */

	void add_pages(char *pages[], int num) {
		assert(num + num_pages <= CELL_SIZE);
		for (int i = 0; i < CELL_SIZE; i++) {
			if (buf[i].get_data() == NULL)
				buf[i] = T(-1, pages[i]);
		}
		num_pages += num;
		int j = 0;
		for (int i = 0; i < CELL_SIZE; i++) {
			if (buf[i].get_data())
				maps[j++] = i;
		}
		assert(j == num_pages);
	}

	/**
	 * Steal the specified number of pages and the pages are
	 * stored in the array.
	 */
	void steal_page(T *pg) {
		assert(pg->get_ref() == 0);
		num_pages--;
		int j = 0;
		for (int i = 0; i < CELL_SIZE; i++) {
			if (buf[i].get_data())
				maps[j++] = i;
		}
		assert(j == num_pages);
	}

	unsigned int get_num_pages() const {
		return num_pages;
	}

	unsigned int get_min_num_pages() const {
		return min_num_pages;
	}

	/**
	 * return an empty page.
	 * I expected the page will be filled with data,
	 * so I change the begin and end index of the circular buffer.
	 */
	T *get_empty_page() {
		T *ret;
		// the condition check in the while loop should be false
		// most of time, and it only iterates once.
		while (idx >= num_pages)
			idx -= num_pages;
		ret = get_page(idx);
		idx++;
		return ret;
	}

	/**
	 * return a page pointed by the iterator or after the iterator.
	 */
	T *get_page(int i) {
		int real_idx = maps[i];
		T *ret = &buf[real_idx];
		assert(ret->get_data());
		return ret;
	}

	int get_idx(T *page) const {
		int idx = page - buf;
		assert (idx >= 0 && idx < num_pages);
		return idx;
	}

	void scale_down_hits() {
		for (int i = 0; i < num_pages; i++) {
			T *pg = get_page(i);
			pg->set_hits(pg->get_hits() / 2);
		}
	}
};

class eviction_policy
{
public:
	thread_safe_page *evict_page(page_cell<thread_safe_page> &buf);
	void access_page(thread_safe_page *pg,
			page_cell<thread_safe_page> &buf) {
		// We don't need to do anything if a page is accessed for many policies.
	}
	bool expand_buffer(const thread_safe_page &pg) {
		return false;
	}
};

class LRU_eviction_policy: public eviction_policy
{
	std::vector<int> pos_vec;
public:
	LRU_eviction_policy() {
		static bool has_print = false;
		if (!has_print)
			printf("use LRU eviction policy\n");
		has_print = true;
	}
	thread_safe_page *evict_page(page_cell<thread_safe_page> &buf);
	void access_page(thread_safe_page *pg,
			page_cell<thread_safe_page> &buf);
	bool expand_buffer(const thread_safe_page &pg) {
		return pg.get_hits() > 0;
	}
};

class clock_eviction_policy: public eviction_policy
{
	unsigned int clock_head;
public:
	clock_eviction_policy() {
		clock_head = 0;
		static bool has_print = false;
		if (!has_print)
			printf("use clock eviction policy\n");
		has_print = true;
	}

	thread_safe_page *evict_page(page_cell<thread_safe_page> &buf);
};

class gclock_eviction_policy: public eviction_policy
{
	unsigned int clock_head;
public:
	gclock_eviction_policy() {
		clock_head = 0;
		static bool has_print = false;
		if (!has_print)
			printf("use gclock eviction policy\n");
		has_print = true;
	}

	thread_safe_page *evict_page(page_cell<thread_safe_page> &buf);
};

class LFU_eviction_policy: public eviction_policy
{
public:
	LFU_eviction_policy() {
		static bool has_print = false;
		if (!has_print)
			printf("use LFU eviction policy\n");
		has_print = true;
	}
	thread_safe_page *evict_page(page_cell<thread_safe_page> &buf);
};

class FIFO_eviction_policy: public eviction_policy
{
public:
	FIFO_eviction_policy() {
		static bool has_print = false;
		if (!has_print)
			printf("use FIFO eviction policy\n");
		has_print = true;
	}
	thread_safe_page *evict_page(page_cell<thread_safe_page> &buf);
};

class associative_cache;

class hash_cell
{
	enum {
		CELL_OVERFLOW,
		IN_QUEUE,
	};
	// It's actually a virtual index of the cell on the hash table.
	int hash;
	atomic_flags<int> flags;

	pthread_spinlock_t _lock;
	page_cell<thread_safe_page> buf;
	associative_cache *table;
#ifdef USE_LRU
	LRU_eviction_policy policy;
#elif defined USE_LFU
	LFU_eviction_policy policy;
#elif defined USE_FIFO
	FIFO_eviction_policy policy;
#elif defined USE_CLOCK
	clock_eviction_policy policy;
#elif defined USE_GCLOCK
	gclock_eviction_policy policy;
#endif
#ifdef USE_SHADOW_PAGE
	clock_shadow_cell shadow;
#endif

	thread_safe_page *get_empty_page();

public:
	hash_cell() {
		table = NULL;
		hash = -1;
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
		return flags.test_flag(OVERFLOW);
	}

	page *search(off_t off, off_t &old_off);
	page *search(off_t offset);
	/* 
	 * this is to rehash the pages in the current cell
	 * to the cell in the parameter.
	 */
	void rehash(hash_cell *cell);

	void print_cell() {
		for (unsigned int i = 0; i < buf.get_num_pages(); i++)
			printf("%lx\t", buf.get_page(i)->get_offset());
		printf("\n");
	}

	void get_dirty_pages(std::map<off_t, thread_safe_page *> &pages);

	long get_hash() const {
		return hash;
	}

	bool is_in_queue() const {
		return flags.test_flag(IN_QUEUE);
	}

	bool set_in_queue(bool v) {
		if (v)
			return flags.set_flag(IN_QUEUE);
		else
			return flags.clear_flag(IN_QUEUE);
	}

	int num_pages(char set_flags, char clear_flags);
};

class flush_thread;
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
	// The number of pages in the cache.
	// Cells may have different numbers of pages, so we can't
	// deduct the size of cache from `ncells'.
	atomic_integer npages;
	
	seq_lock table_lock;
	atomic_flags<int> flags;
	/* the initial number of cells in the table. */
	int init_ncells;

	memory_manager *manager;

	bool expandable;
	/* used for linear hashing */
	int level;
	int split;

	flush_thread *_flush_thread;

public:
	associative_cache(long cache_size, bool expandable = false);

	~associative_cache() {
		for (unsigned int i = 0; i < cells_table.size(); i++)
			if (cells_table[i])
				delete [] cells_table[i];
		manager->unregister_cache(this);
	}

	/* the hash function used for the current level. */
	int hash(off_t offset) {
		return offset / PAGE_SIZE % (init_ncells * (long) (1 << level));
	}

	/* the hash function used for the next level. */
	int hash1(off_t offset) {
		return offset / PAGE_SIZE % (init_ncells * (long) (1 << (level + 1)));
	}

	int hash1_locked(off_t offset) {
		unsigned long count;
		int ret;
		do {
			table_lock.read_lock(count);
			ret = offset / PAGE_SIZE % (init_ncells * (long) (1 << (level + 1)));
		} while (!table_lock.read_unlock(count));
		return ret;
	}

	/**
	 * This method searches for the specified page. If the page doesn't exist,
	 * it tries to evict a page. Therefore, it also triggers some code of
	 * maintaining eviction policy.
	 */
	page *search(off_t offset, off_t &old_off);
	/**
	 * This method just searches for the specified page, nothing more.
	 * So if the request isn't issued from the workload, and we don't need
	 * evict a page if the specified page doesn't exist, we should use
	 * this method.
	 */
	page *search(off_t offset);

	bool expand(hash_cell *cell);

	bool shrink(int npages, char *pages[]) {
		// TODO shrink the cache
		return false;
	}

	memory_manager *get_manager() {
		return manager;
	}

	void print_cell(off_t off) {
		get_cell(off)->print_cell();
	}

	/**
	 * The size of allocated pages in the cache.
	 */
	long size() {
		return ((long) npages.get()) * PAGE_SIZE;
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

	bool is_expandable() const {
		return expandable;
	}

	/* Methods for flushing dirty pages. */

	flush_thread *create_flush_thread(io_interface *io);
	flush_thread *get_flush_thread() const {
		return _flush_thread;
	}

	hash_cell *get_prev_cell(hash_cell *cell);
	hash_cell *get_next_cell(hash_cell *cell);
};

#endif
