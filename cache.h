#ifndef __CACHE_H__
#define __CACHE_H__

#include <unistd.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>
#include <assert.h>
#include <numa.h>

#include <map>

#include "container.h"
#include "concurrency.h"

#define PTHREAD_WAIT

#define PAGE_SIZE 4096
#define LOG_PAGE_SIZE 12
#define ROUND_PAGE(off) (((long) off) & (~(PAGE_SIZE - 1)))

enum {
	DATA_READY_BIT = 0,
	IO_PENDING_BIT,
	DIRTY_BIT,
	LOCK_BIT,
	ACTIVE_BIT,
	REFERENCED_BIT,
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
protected:
	volatile void *data;
	volatile short refcnt;
	volatile char flags;
	volatile unsigned char hits;

public:
	page(): data(NULL) {
		offset = -1;
		refcnt = 0;
		flags = 0;
		hits = 0;
	}

	page(off_t off, char *data) {
		set_offset(off);
		this->data = data;
		refcnt = 0;
		flags = 0;
		hits = 0;
	}

	page(off_t off, long data_off) {
		set_offset(off);
		data = (void *) ((long) data_start + data_off);
		refcnt = 0;
		flags = 0;
		hits = 0;
	}

	/* offset in the file in bytes */
	void set_offset(off_t off) {
		offset = off >> LOG_PAGE_SIZE;
	}

	bool initialized() const {
		return offset != -1;
	}

	off_t get_offset() const { return ((off_t) offset) << LOG_PAGE_SIZE; }
	void *get_data() const { return (void *) data; }

	bool data_ready() const { return flags & (0x1 << DATA_READY_BIT); }
	void set_data_ready(bool ready) {
		if (ready)
			flags |= 0x1 << DATA_READY_BIT;
		else
			flags &= ~(0x1 << DATA_READY_BIT);
	}
	bool is_dirty() const { return flags & (0x1 << DIRTY_BIT); }
	void set_dirty(bool dirty) {
		if (dirty)
			flags |= 0x1 << DIRTY_BIT;
		else
			flags &= ~(0x1 << DIRTY_BIT);
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

	void set_active(bool active) {
		if (active)
			flags |= 0x1 << ACTIVE_BIT;
		else
			flags &= ~(0x1 << ACTIVE_BIT);
	}
	bool active() {
		return flags & (0x1 << ACTIVE_BIT);
	}

	void reset_hits() {
		hits = 0;
	}
	int get_hits() {
		return hits;
	}
	void set_hits(int hits) {
		assert (hits <= 0xff);
		this->hits = hits;
	}
	/* the page is accessed */
	void hit() {
		hits++;
	}

	static void allocate_cache(long size) {
		data_start = numa_alloc_local(size);
	}

	virtual void inc_ref() {
		refcnt++;
	}
	virtual void dec_ref() {
		refcnt--;
	}
	virtual short get_ref() {
		return refcnt;
	}
};

class thread_safe_page: public page
{
#ifdef PTHREAD_WAIT
	pthread_cond_t ready_cond;
	pthread_cond_t dirty_cond;
	pthread_mutex_t mutex;
#endif

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
	thread_safe_page(): page() {
#ifdef PTHREAD_WAIT
		pthread_cond_init(&ready_cond, NULL);
		pthread_cond_init(&dirty_cond, NULL);
		pthread_mutex_init(&mutex, NULL);
#endif
	}

	thread_safe_page(off_t off, long d): page(off, d) {
#ifdef PTHREAD_WAIT
		pthread_cond_init(&ready_cond, NULL);
		pthread_cond_init(&dirty_cond, NULL);
		pthread_mutex_init(&mutex, NULL);
#endif
	}

	thread_safe_page(off_t off, char *data): page(off, data) {
#ifdef PTHREAD_WAIT
		pthread_cond_init(&ready_cond, NULL);
		pthread_cond_init(&dirty_cond, NULL);
		pthread_mutex_init(&mutex, NULL);
#endif
	}

	~thread_safe_page() {
#ifdef PTHREAD_WAIT
		pthread_mutex_destroy(&mutex);
		pthread_cond_destroy(&ready_cond);
		pthread_cond_destroy(&dirty_cond);
#endif
	}

	/* this is enough for x86 architecture */
	bool data_ready() const { return get_flags_bit(DATA_READY_BIT); }
	void wait_ready() {
#ifdef PTHREAD_WAIT
		printf("wait data to be ready\n");
		pthread_mutex_lock(&mutex);
#endif
		while (!data_ready()) {
#ifdef PTHREAD_WAIT
			pthread_cond_wait(&ready_cond, &mutex);
			printf("thread %ld wait for data ready\n", pthread_self());
#endif
		}
#ifdef PTHREAD_WAIT
		pthread_mutex_unlock(&mutex);
#endif
	}
	void set_data_ready(bool ready) {
#ifdef PTHREAD_WAIT
		pthread_mutex_lock(&mutex);
#endif
		set_flags_bit(DATA_READY_BIT, ready);
#ifdef PTHREAD_WAIT
		pthread_cond_signal(&ready_cond);
		pthread_mutex_unlock(&mutex);
#endif
	}

	void set_dirty(bool dirty) {
#ifdef PTHREAD_WAIT
		pthread_mutex_lock(&mutex);
#endif
		set_flags_bit(DIRTY_BIT, dirty);
#ifdef PTHREAD_WAIT
		pthread_cond_signal(&dirty_cond);
		pthread_mutex_unlock(&mutex);
#endif
	}
	void wait_cleaned() {
#ifdef PTHREAD_WAIT
		pthread_mutex_lock(&mutex);
#endif
		while (is_dirty()) {
#ifdef PTHREAD_WAIT
			pthread_cond_wait(&dirty_cond, &mutex);
			printf("thread %ld wait for the page to be written\n",
					pthread_self());
#endif
		}
#ifdef PTHREAD_WAIT
		pthread_mutex_unlock(&mutex);
#endif
	}

	bool is_io_pending() {
		return get_flags_bit(IO_PENDING_BIT);
	}
	void set_io_pending(bool pending) {
		set_flags_bit(IO_PENDING_BIT, pending);
	}
	/* we set the status to io pending,
	 * and return the original status */
	bool test_and_set_io_pending() {
		int old = __sync_fetch_and_or(&flags, 0x1 << IO_PENDING_BIT);
		return old & (0x1 << IO_PENDING_BIT);
	}

	void lock() {
		int old;
		do {
			old = __sync_fetch_and_or(&flags, 0x1 << LOCK_BIT);
		} while (old & (0x1 << LOCK_BIT));
	}
	void unlock() {
		set_flags_bit(LOCK_BIT, false);
	}

	void inc_ref() {
		__sync_fetch_and_add(&refcnt, 1);
	}
	void dec_ref() {
		__sync_fetch_and_sub(&refcnt, 1);
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
	virtual page *search(off_t offset, off_t &old_off) = 0;
	virtual long size() {
		return 0;
	}
	virtual bool shrink(int npages, char *pages[]) {
		return false;
	}
};

/**
 * This data structure is to implement LRU.
 */
template<class T>
class page_buffer: public fifo_queue<T>
{
public:
	/*
	 * @size: the size of the page buffer
	 * @page_buf: the offset of the page array in the global page cache.
	 */
	page_buffer(long size, long page_buf): fifo_queue<T>(size) {
		for (int i = 0; i < size; i++) {
			fifo_queue<T>::buf[i] = T(-1, page_buf + i * PAGE_SIZE);
		}
	}

	/**
	 * return an empty page.
	 * I expected the page will be filled with data,
	 * so I change the begin and end index of the circular buffer.
	 */
	T *get_empty_page() {
		return fifo_queue<T>::get_empty_entry();
	}

	T *get_page(int i) {
		return fifo_queue<T>::get_entry(i);
	}
};

/**
 * This is used to implement hashtable-based indexing.
 */
class frame: public thread_safe_page
{
	/* equivalent to the number of hits */
	atomic_integer wcount;
	/* equivalent to the number of references */
	atomic_integer pinning;
	frame *prev, *next;

public:
	frame(): thread_safe_page() {
		prev = next = this;
	}

	frame(off_t offset, char *data): thread_safe_page(offset, data) {
		prev = next = this;
	}

	frame(off_t off, long d): thread_safe_page(off, d) {
		prev = next = this;
	}

	void *volatileGetValue() {
		/*
		 * maybe a memory fence here,
		 * but it's not needed in the x86 machine.
		 */
		return get_data();
	}

	bool CASValue(void *expect,void *update) {
		return __sync_bool_compare_and_swap(&data, expect, update);
	}

	int getWC() {
		return wcount.get();
	}

	void incrWC(int num = 1) {
		wcount.inc(num);
	}

	int decrWC(int num = 1) {
		return wcount.dec(num);
	}

	bool tryEvict() {
		return pinning.CAS(0, -1);
	}

	void evictUnshared() {
		pinning.CAS(1, -1);
	}

	int pinCount() {
		return pinning.get();
	}

	bool pin() {
		int x;
		do {
			x = pinning.get();
			if (x <= -1)
				return false;
		} while (!pinning.CAS(x, x + 1));
		return true;
	}

	void dec_ref() {
		unpin();
	}

	void unpin() {
		pinning.dec(1);
	}

	/* for linked pages */

private:
	/* Add a frame behind the frame in the list. */
	void add_front(frame *pg) {
		frame *next = this->next;
		pg->next = next;
		pg->prev = this;
		this->next = pg;
		next->prev = pg;
	}

	/* Add a frame before the frame in the list. */
	void add_back(frame *pg) {
		frame *prev = this->prev;
		pg->next = this;
		pg->prev = prev;
		this->prev = pg;
		prev->next = pg;
	}

	void remove_from_list() {
		frame *prev = this->prev;
		frame *next = this->next;
		prev->next = next;
		next->prev = prev;
		this->next = this;
		this->prev = this;
	}

	bool is_empty() {
		return this->next == this;
	}

public:
	frame *front() {
		return next;
	}

	frame *back() {
		return prev;
	}

	friend class linked_page_queue;
};

class linked_page_queue {
	frame head;
	int _size;
public:
	linked_page_queue() {
		_size = 0;
	}

	bool is_head(frame *pg) {
		return pg == &head;
	}

	void push_back(frame *pg) {
		head.add_back(pg);
		_size++;
	}

	void pop_front() {
		frame *pg = front();
		remove(pg);
	}

	void remove(frame *pg) {
		if (pg->is_empty())
			return;
		// TODO do I need to check whether the page is in the queue.
		pg->remove_from_list();
		_size--;
	}

	void remove(int idx) {
		frame *f = front();
		for (int i = 0; i < idx; i++)
			f = f->front();
		f->remove_from_list();
		_size--;
	}

	bool empty() {
		return head.is_empty();
	}

	frame *front() {
		return head.front();
	}

	frame *back() {
		return head.back();
	}

	int size() {
		return _size;
	}

	void replace(frame *old_frame, frame *new_frame) {
		frame *prev = old_frame->back();
		old_frame->remove_from_list();
		prev->add_front(new_frame);
	}

	void merge(linked_page_queue *list) {
		if (list->empty())
			return;

		frame *list_begin = list->front();
		frame *list_end = list->back();
		frame *this_end = this->back();
		
		/* Remove all pages in the list. */
		list->head.next = &list->head;
		list->head.prev = &list->head;
		this->_size += list->_size;
		list->_size = 0;

		/* Add `list' to the end of this list. */
		this_end->next = list_begin;
		list_begin->prev = this_end;

		/* Link the end of `list' to the end of this list. */
		list_end->next = &this->head;
		this->head.prev = list_end;
	}

	void print() {
		frame *f = front();
		printf("queue size: %d\n", _size);
		while (f != &head) {
			printf("%ld\t", f->get_offset());
			f = f->front();
		}
		printf("\n");
	}
};

#endif
