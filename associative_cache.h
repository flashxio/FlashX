#ifndef __ASSOCIATIVE_CACHE_H__
#define __ASSOCIATIVE_CACHE_H__

#include "cache.h"

#define CELL_SIZE 8

#ifdef STATISTICS
static volatile int avail_cells;
static volatile int num_wait_unused;
static volatile int lock_contentions;
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

	void set_pages(long page_buf) {
		for (int i = 0; i < BUF_SIZE; i++) {
			buf[i] = T(-1, page_buf + i * PAGE_SIZE);
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
};

// TODO the entire cell should be put in the same cache line
// so each access to the hash table has only one cache miss.
class hash_cell
{
	pthread_spinlock_t _lock;
	page_cell<thread_safe_page, CELL_SIZE> buf;
	char stuffing[CACHE_LINE - sizeof(_lock) - sizeof(buf)];

	/* this function has to be called with lock held */
	thread_safe_page *get_empty_page() {
		thread_safe_page *ret = NULL;

		do {
			int min_hits = 0x7fffffff;
			for (int i = 0; i < CELL_SIZE; i++) {
				thread_safe_page *pg = buf.get_page(i);
				if (pg->get_ref())
					continue;

				/* 
				 * refcnt only increases within the lock of the cell,
				 * so if the page's refcnt is 0 above,
				 * it'll be always 0 within the lock.
				 */

				if (min_hits > pg->get_hits()) {
					min_hits = pg->get_hits();
					ret = pg;
				}
			}
			/* it happens when all pages in the cell is used currently. */
		} while (ret == NULL);

		ret->reset_hits();
		ret->set_data_ready(false);
		return ret;
	}

public:
	hash_cell() {
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

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

	void set_pages(long page_buf) {
		buf.set_pages(page_buf);
	}

	/**
	 * search for a page with the offset.
	 * If the page doesn't exist, return an empty page.
	 */
	page *search(off_t off, off_t &old_off) {
		thread_safe_page *ret = NULL;
#ifndef STATISTICS
		pthread_spin_lock(&_lock);
#else
		if (pthread_spin_trylock(&_lock) == EBUSY) {
			__sync_fetch_and_add(&lock_contentions, 1);
			pthread_spin_lock(&_lock);
		}
#endif

		for (int i = 0; i < CELL_SIZE; i++) {
			if (buf.get_page(i)->get_offset() == off) {
				ret = buf.get_page(i);
				break;
			}
		}
		if (ret == NULL) {
			ret = get_empty_page();
			old_off = ret->get_offset();
			/*
			 * I have to change the offset in the spinlock,
			 * to make sure when the spinlock is unlocked, 
			 * the page can be seen by others even though
			 * it might not have data ready.
			 */
			ret->set_offset(off);
		}
		/* it's possible that the data in the page isn't ready */
		ret->inc_ref();
		pthread_spin_unlock(&_lock);
		return ret;
	}

	void print_cell() {
		for (int i = 0; i < CELL_SIZE; i++)
			printf("%lx\t", buf.get_page(i)->get_offset());
		printf("\n");
	}
};

class associative_cache: public page_cache
{
	hash_cell *cells;
	int ncells;
	// TODO maybe it's not a good hash function
	int hash(off_t offset) {
		return offset / PAGE_SIZE % ncells;
	}

public:
	associative_cache(long cache_size) {
		printf("associative cache is used\n");
		int npages = cache_size / PAGE_SIZE;
		ncells = npages / CELL_SIZE;
		cells = new hash_cell[ncells];
		printf("%d cells: %p\n", ncells, cells);
		for (int i = 0; i < ncells; i++)
			cells[i].set_pages(i * PAGE_SIZE * CELL_SIZE);
	}

	~associative_cache() {
		delete [] cells;
	}

	page *search(off_t offset, off_t &old_off) {
		hash_cell *cell = &cells[hash(offset)];
		return cell->search(offset, old_off);
	}

	void print_cell(off_t off) {
		cells[hash(off)].print_cell();
	}
};

#endif
