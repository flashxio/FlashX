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

	void scale_down_hits() {
		for (int i = 0; i < BUF_SIZE; i++) {
			buf[i].set_hits(buf[i].get_hits() / 2);
		}
	}
};

#define SHADOW_FACTOR 4
class shadow_cell
{
	unsigned int num;
	page buf[CELL_SIZE * SHADOW_FACTOR];
public:
	shadow_cell() {
		num = 0;
	}

	/*
	 * add a page to the cell, which is evicted from hash_cell.
	 * the only thing we need to record is the number of hits
	 * of this page.
	 */
	void add(page &pg) {
		/* if the cell isn't full */
		if (num < sizeof (buf) / sizeof (page)) {
			buf[num] = pg;
			num++;
		}
		else {
			int idx = -1;
			int min_hits = 0x7fffffff;
			for (unsigned int i = 0; i < num; i++)
				if (buf[i].get_hits() < min_hits) {
					min_hits = buf[i].get_hits();
					idx = i;
				}
			assert(idx != -1);
			/*
			 * If the cell is full and the number of hits of the new page
			 * is higher than a page in the cell, we need to
			 * evicted the shadown page from the cell.
			 */
			if (min_hits < pg.get_hits())
				buf[idx] = pg;
			/*
			 * If the cell is full and the number of hits of the new page
			 * is lower than any pages in the cell. It's ignored.
			 */
		}
	}

	page *search(off_t off) {
		for (unsigned int i = 0; i < sizeof(buf) / sizeof(page); i++) {
			if (buf[i].get_offset() == off)
				return &buf[i];
		}
		return NULL;
	}

	void scale_down_hits() {
		for (int i = 0; i < CELL_SIZE * SHADOW_FACTOR; i++) {
			buf[i].set_hits(buf[i].get_hits() / 2);
		}
	}
};

// TODO the entire cell should be put in the same cache line
// so each access to the hash table has only one cache miss.
class hash_cell
{
	pthread_spinlock_t _lock;
	page_cell<thread_safe_page, CELL_SIZE> buf;
	shadow_cell shadow;
//	char stuffing[CACHE_LINE - sizeof(_lock) - sizeof(buf)];

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

		/* we record the hit info of the page in the shadow cell. */
		if (ret->get_hits() > 0)
			shadow.add(*ret);

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
			page *shadow_pg = shadow.search(off);
			/*
			 * if the page has been seen before,
			 * we should set the hits info.
			 */
			if (shadow_pg)
				ret->set_hits(shadow_pg->get_hits());
		}
		/* it's possible that the data in the page isn't ready */
		ret->inc_ref();
		if (ret->get_hits() == 0xff) {
			buf.scale_down_hits();
			shadow.scale_down_hits();
		}
		ret->hit();
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
