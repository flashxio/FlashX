#ifndef __CACHE_H__
#define __CACHE_H__

#include <unistd.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>

#include <map>

#define PAGE_SIZE 4096
#define ROUND_PAGE(off) (((long) off) & (~(PAGE_SIZE - 1)))

class page
{
	off_t offset;
	void *data;
	char flags;

public:
	page():offset(-1), data(NULL) { }

	page(off_t off, void *d): offset(off), data(d) { }

	void set_offset(off_t off) {
		offset = off;
	}

	off_t get_offset() const { return offset; }
	void *get_data() const { return data; }

	bool data_ready() const { return flags & 0x1; }
	void set_data_ready(bool ready) {
		if (ready)
			flags |= 0x1;
		else
			flags &= ~0x1;
	}
};

class thread_safe_page: public page
{
	volatile char flags;
	volatile short refcnt;
	volatile char io_pending;
public:
	thread_safe_page(): page(), flags(0), refcnt(0), io_pending(0) {
	}

	thread_safe_page(off_t off, void *d): page(off, d), flags(0), refcnt(0), io_pending(0) {
	}

	/* this is enough for x86 architecture */
	bool data_ready() const { return flags & 0x1; }
	void wait_ready() {
		while (!data_ready()) {}
	}
	void set_data_ready(bool ready) {
		if (ready)
			__sync_fetch_and_or(&flags, 0x1);
		else
			__sync_fetch_and_and(&flags, ~0x1);
	}

	void set_io_pending(bool pending) {
		io_pending = pending;
	}
	/* we set the status to io pending,
	 * and return the original status */
	bool test_and_set_io_pending() {
		return !__sync_bool_compare_and_swap(&io_pending, false, true);
	}

	void inc_ref() {
		__sync_fetch_and_add(&refcnt, 1);
	}
	void dec_ref() {
		__sync_fetch_and_sub(&refcnt, 1);
	}
	short get_ref() {
		return refcnt;
	}
	void wait_unused() {
		while(get_ref()) {}
	}
};

class page_cache
{
	pthread_spinlock_t _lock;
public:
	page_cache() {
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

	int lock() {
		return pthread_spin_lock(&_lock);
	}

	int unlock() {
		return pthread_spin_unlock(&_lock);
	}

	~page_cache() {
		pthread_spin_destroy(&_lock);
	}
	virtual page *search(off_t offset) = 0;
};

/**
 * This data structure is to implement LRU.
 */
template<class T>
class page_buffer
{
	long size;			// the number of pages that can be buffered
	T *buf;			// a circular buffer to keep pages.
	char *pages;		// the pointer to the space for all pages.
	int beg_idx;		// the index of the beginning of the buffer
	int end_idx;		// the index of the end of the buffer

public:
	page_buffer(long size) {
		this->size = size;
		buf = new T[size];
		pages = (char *) valloc(size * PAGE_SIZE);
		for (int i = 0; i < size; i++) {
			buf[i] = T(-1, &pages[i * PAGE_SIZE]);
		}
		beg_idx = 0;
		end_idx = 0;
	}

	~page_buffer() {
		delete [] buf;
		free(pages);
	}

	bool is_full() const {
		return end_idx - beg_idx == -1
			|| end_idx - beg_idx == size - 1;
	}

	/**
	 * return an empty page.
	 * I expected the page will be filled with data,
	 * so I change the begin and end index of the circular buffer.
	 */
	T *get_empty_page() {
		if (is_full()) {
			beg_idx = (beg_idx + 1) % size;
		}
		T *ret = &buf[end_idx];
		end_idx = (end_idx + 1) % size;
		return ret;
	}

	T *search(off_t off) {
		T *ret = NULL;
		for (int i = beg_idx; i != end_idx;
				i = (i + 1) % size) {
			if (buf[i].get_offset() == off) {
				ret = &buf[i];
				break;
			}
		}
		return ret;
	}

#if 0
	/* push a page to the buffer.
	 * if the buffer is full, the first page is removed,
	 * and return the user */
	page *push_back(const page &page, page &beg) {
		page *ret;
		if (is_full()) {
			beg = buf[beg_idx];
			beg_idx = (beg_idx + 1) % size;
		}
		buf[end_idx] = page;
		ret = &buf[end_idx];
		end_idx = (end_idx + 1) % size;
		return ret;
	}
#endif
};

#endif
