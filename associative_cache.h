#ifndef __ASSOCIATIVE_CACHE_H__
#define __ASSOCIATIVE_CACHE_H__

#include "cache.h"

#define CELL_SIZE 4

#ifdef STATISTICS
static volatile int avail_cells;
static int num_wait_unused;
#endif

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

	/* this function has to be called with lock held */
	thread_safe_page *get_empty_page() {
		thread_safe_page *ret = buf.get_empty_page();
		/*
		 * each time we select a page to evict,
		 * it's possible that it's still used by some
		 * other threads. We need to wait for other threads
		 * to finish using it.
		 */
		// TODO I have to improve here.
		// by removing this part of code, I can get 2s faster.
		while (ret->get_ref()) {
#ifdef STATISTICS
			num_wait_unused++;
#endif
			pthread_spin_unlock(&_lock);
			ret->wait_unused();
			pthread_spin_lock(&_lock);
		}
		ret->set_data_ready(false);
		ret->inc_ref();
		return ret;
	}

public:
	hash_cell() {
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

	~hash_cell() {
		pthread_spin_destroy(&_lock);
	}

	void set_pages(long page_buf) {
		buf.set_pages(page_buf);
	}

	/**
	 * search for a page with the offset.
	 * If the page doesn't exist, return an empty page.
	 */
	page *search(off_t off) {
		thread_safe_page *ret = NULL;
		pthread_spin_lock(&_lock);

		for (int i = 0; i < CELL_SIZE; i++) {
			if (buf.get_page(i)->get_offset() == off) {
				ret = buf.get_page(i);
				break;
			}
		}
		if (ret == NULL) {
			ret = get_empty_page();
			ret->set_offset(off);
			pthread_spin_unlock(&_lock);
			return ret;
		}
		/* it's possible that the data in the page isn't ready */
		ret->inc_ref();
		pthread_spin_unlock(&_lock);
		return ret;
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
		for (int i = 0; i < ncells; i++)
			cells[i].set_pages(i * PAGE_SIZE * CELL_SIZE);
	}

	~associative_cache() {
		delete [] cells;
	}

	page *search(off_t offset) {
		hash_cell *cell = &cells[hash(offset)];
		return cell->search(offset);
	}
};

#endif
