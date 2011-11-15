#ifndef __CACHE_H__
#define __CACHE_H__

#include <unistd.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>

#include <map>

#define PAGE_SIZE 4096
#define LOG_PAGE_SIZE 12
#define ROUND_PAGE(off) (((long) off) & (~(PAGE_SIZE - 1)))

enum {
	DATA_READY_BIT,
	IO_PENDING_BIT,
};

class page
{
	/*
	 * the offset in the file in pages.
	 * it can cover a file of 8 TB.
	 */
	int offset;

	/* 
	 * all data in a page is in a buffer,
	 * so can we use the start of the buffer
	 * and the offset in the buffer to calculate
	 * the address of the page.
	 */
	static void *data_start;
	/*
	 * in pages.
	 */
	int buf_offset;

protected:
	volatile char flags;

public:
	page():offset(-1), buf_offset(0), flags(0) { }

	page(off_t off, long data) {
		set_offset(off);
		flags = 0;
		buf_offset = data >> LOG_PAGE_SIZE;
	}

	/* offset in the file in bytes */
	void set_offset(off_t off) {
		offset = off >> LOG_PAGE_SIZE;
	}

	bool initialized() const {
		return offset != -1;
	}

	// TODO is off_t unsigned?
	off_t get_offset() const { return ((off_t) offset) << LOG_PAGE_SIZE; }
	void *get_data() const { return (void *) ((long) data_start
			+ (((long) buf_offset) << LOG_PAGE_SIZE)); }

	bool data_ready() const { return flags & (0x1 << DATA_READY_BIT); }
	void set_data_ready(bool ready) {
		if (ready)
			flags |= 0x1 << DATA_READY_BIT;
		else
			flags &= ~(0x1 << DATA_READY_BIT);
	}

	static void allocate_cache(long size) {
		data_start = valloc(size);
	}
};
void *page::data_start;

class thread_safe_page: public page
{
	volatile unsigned char refcnt;

	void set_flags_bit(int i, bool v) {
		if (v)
			__sync_fetch_and_or(&flags, 0x1 << i);
		else
			__sync_fetch_and_and(&flags, ~(0x1 << i));
	}

	bool get_flags_bit(int i) const {
		return flags & (0x1 << i);
	}

public:
	thread_safe_page(): page(), refcnt(0) {
	}

	thread_safe_page(off_t off, long d): page(off, d), refcnt(0) {
	}

	/* this is enough for x86 architecture */
	bool data_ready() const { return get_flags_bit(DATA_READY_BIT); }
	void wait_ready() {
		while (!data_ready()) {}
	}
	void set_data_ready(bool ready) {
		set_flags_bit(DATA_READY_BIT, ready);
	}

	void set_io_pending(bool pending) {
		set_flags_bit(IO_PENDING_BIT, pending);
	}
	/* we set the status to io pending,
	 * and return the original status */
	bool test_and_set_io_pending() {
		char old = __sync_fetch_and_or(&flags, 0x1 << IO_PENDING_BIT);
		return old & (0x1 << IO_PENDING_BIT);
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
	int beg_idx;		// the index of the beginning of the buffer
	int end_idx;		// the index of the end of the buffer

public:
	/*
	 * @size: the size of the page buffer
	 * @page_buf: the offset of the page array in the global page cache.
	 */
	page_buffer(long size, long page_buf) {
		this->size = size;
		buf = new T[size];
		for (int i = 0; i < size; i++) {
			buf[i] = T(-1, page_buf + i * PAGE_SIZE);
		}
		beg_idx = 0;
		end_idx = 0;
	}

	~page_buffer() {
		delete [] buf;
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
