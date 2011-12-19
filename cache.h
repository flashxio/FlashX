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
	FLAGS_NBITS = 6,
	BUF_OFFSET_NBITS = 24,
	DATA_READY_BIT = 30,
	IO_PENDING_BIT = 31,
};

const int BUF_OFFSET_BITS = (1 << BUF_OFFSET_NBITS) - 1;

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
protected:
	volatile int buf_offset;

	int get_buf_offset() const {
		return buf_offset & BUF_OFFSET_BITS;
	}
	void set_buf_offset(int off) {
		int flags = buf_offset & (~BUF_OFFSET_BITS);
		buf_offset = flags | (off & BUF_OFFSET_BITS);
	}

public:
	page():offset(-1), buf_offset(0) { }

	page(off_t off, long data) {
		set_offset(off);
		set_buf_offset(data >> LOG_PAGE_SIZE);
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
			+ (((long) get_buf_offset()) << LOG_PAGE_SIZE)); }

	bool data_ready() const { return buf_offset & (0x1 << DATA_READY_BIT); }
	void set_data_ready(bool ready) {
		if (ready)
			buf_offset |= 0x1 << DATA_READY_BIT;
		else
			buf_offset &= ~(0x1 << DATA_READY_BIT);
	}

	static void allocate_cache(long size) {
		data_start = valloc(size);
	}
};
void *page::data_start;

class thread_safe_page: public page
{
	void set_flags_bit(int i, bool v) {
		if (v)
			__sync_fetch_and_or(&buf_offset, 0x1 << i);
		else
			__sync_fetch_and_and(&buf_offset, ~(0x1 << i));
	}

	bool get_flags_bit(int i) const {
		return buf_offset & (0x1 << i);
	}

public:
	thread_safe_page(): page() {
	}

	thread_safe_page(off_t off, long d): page(off, d) {
	}

	/* this is enough for x86 architecture */
	bool data_ready() const { return get_flags_bit(DATA_READY_BIT); }
	void wait_ready() {
		while (!data_ready()) {
#ifdef DEBUG
			printf("thread %ld wait for data ready\n", pthread_self());
#endif
		}
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
		int old = __sync_fetch_and_or(&buf_offset, 0x1 << IO_PENDING_BIT);
		return old & (0x1 << IO_PENDING_BIT);
	}

	void inc_ref() {
		char *refcnt = ((char *) &buf_offset) + 3;
		__sync_fetch_and_add(refcnt, 1);
	}
	void dec_ref() {
		char *refcnt = ((char *) &buf_offset) + 3;
		__sync_fetch_and_sub(refcnt, 1);
	}
	short get_ref() {
		char *refcnt = ((char *) &buf_offset) + 3;
		return (*refcnt) & ((1 << FLAGS_NBITS) - 1);
	}
	void wait_unused() {
		while(get_ref()) {
#ifdef DEBUG
			printf("thread %ld wait for used\n", pthread_self());
#endif
		}
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
 * this is a first-in-first-out queue.
 * However, the location of an entry in the queue never changes.
 */
template<class T>
class queue
{
	long size;			// the number of pages that can be buffered
	volatile unsigned long idx;		// to the point where we can evict a page in the buffer
protected:
	T *buf;			// a circular buffer to keep pages.
public:
	queue(long size) {
		idx = 0;
		this->size = size;
		buf = new T[size];
	}

	~queue() {
		delete [] buf;
	}

	T *get_empty_entry() {
		/* TODO I ignore the case of integer overflow */
		long orig = __sync_fetch_and_add(&idx, 1);
		T *ret = &buf[orig % size];
		return ret;
	}

	T *get_entry(int i) {
		if (i >= size)
			return NULL;
		return &buf[i];
	}

	int get_idx(T *p) {
		return p - buf;
	}
};

/**
 * This data structure is to implement LRU.
 */
template<class T>
class page_buffer: public queue<T>
{
public:
	/*
	 * @size: the size of the page buffer
	 * @page_buf: the offset of the page array in the global page cache.
	 */
	page_buffer(long size, long page_buf): queue<T>(size) {
		for (int i = 0; i < size; i++) {
			queue<T>::buf[i] = T(-1, page_buf + i * PAGE_SIZE);
		}
	}

	/**
	 * return an empty page.
	 * I expected the page will be filled with data,
	 * so I change the begin and end index of the circular buffer.
	 */
	T *get_empty_page() {
		return queue<T>::get_empty_entry();
	}

	T *get_page(int i) {
		return queue<T>::get_entry(i);
	}
};

#endif
