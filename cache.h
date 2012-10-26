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

#include "common.h"
#include "container.h"
#include "concurrency.h"
#include "slab_allocator.h"

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
	/* This method should be called within each thread. */
	virtual void init() {
	}
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

public:
	frame(): thread_safe_page() {
	}

	frame(off_t offset, char *data): thread_safe_page(offset, data) {
	}

	frame(off_t off, long d): thread_safe_page(off, d) {
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
};

class linked_page_queue;

class linked_obj
{
	linked_obj *prev, *next;
	void *payload;
public:
	linked_obj() {
		prev = next = this;
		payload = NULL;
	}

	linked_obj(void *payload) {
		prev = next = this;
		this->payload = payload;
	}

	void set_payload(void *payload) {
		this->payload = payload;
	}

	void *get_payload() const {
		return payload;
	}

	/* Add an obj behind the obj in the list. */
	void add_front(linked_obj *obj) {
		linked_obj *next = this->next;
		obj->next = next;
		obj->prev = this;
		this->next = obj;
		next->prev = obj;
	}

	/* Add an obj before the obj in the list. */
	void add_back(linked_obj *obj) {
		linked_obj *prev = this->prev;
		obj->next = this;
		obj->prev = prev;
		this->prev = obj;
		prev->next = obj;
	}

	void remove_from_list() {
		linked_obj *prev = this->prev;
		linked_obj *next = this->next;
		prev->next = next;
		next->prev = prev;
		this->next = this;
		this->prev = this;
	}

	bool is_empty() const {
		return this->next == this;
	}

	linked_obj *front() const {
		return next;
	}

	linked_obj *back() const {
		return prev;
	}

	friend class linked_page_queue;
};

/**
 * The queue is formed as a linked list, so it only supports sequential access.
 * The queue is designed to reduce writes to the shared memory when we need to
 * reorganize the page list. The idea is to split the page metadata from the 
 * data structure that forms the linked list. As long as the linked page list
 * isn't shared by multiple CPUs, any modification on the page list occurs
 * on the local memory.
 */
class linked_page_queue {
	linked_obj head;
	int _size;

	bool local_allocator;
	obj_allocator<linked_obj> *allocator;

protected:
	void remove(linked_obj *obj) {
		obj->remove_from_list();
		_size--;
		allocator->free(obj);
	}
public:
	class iterator {
		linked_obj *curr_loc;
		linked_page_queue *queue;
		int num_iter;	// number of pages that have been accessed.

		iterator(linked_obj *head, linked_page_queue *queue) {
			this->curr_loc = head;
			this->queue = queue;
			num_iter = 0;
		}
	public:
		iterator() {
			curr_loc = NULL;
			queue = NULL;
			num_iter = 0;
		}

		bool has_next() const {
			if (curr_loc == NULL || queue == NULL)
				return false;
			return num_iter < queue->size();
		}

		/* move to the next object and return the next object. */
		frame *next() {
			assert(curr_loc != NULL && queue != NULL);
			curr_loc = curr_loc->front();
			num_iter++;
			return (frame *) curr_loc->get_payload();
		}

		/* 
		 * These methods are extensions of the basic iterator.
		 * They can only be called after next() is called at least once.
		 */

		/*
		 * return the current object.
		 * return NULL if it's pointing to the head of the queue.
		 */
		frame *curr() {
			assert(curr_loc != NULL && queue != NULL);
			if (curr_loc == &queue->head)
				return NULL;
			return (frame *) curr_loc->get_payload();
		}

		void set(frame *payload) {
			assert(curr_loc != NULL && queue != NULL);
			if (curr_loc != &queue->head)
				curr_loc->set_payload(payload);
		}

		/* 
		 * remove the current frame in the queue
		 * and move to the next page automatically.
		 */
		void remove() {
			assert(curr_loc != NULL && queue != NULL);
			if (queue->size() <= 0 && curr_loc != &queue->head)
				return;
			linked_obj *tmp = curr_loc;
			curr_loc = curr_loc->back();
			num_iter--;
			queue->remove(tmp);
		}

		// for test
		linked_page_queue *owner() const {
			return queue;
		}

		friend class linked_page_queue;
	};

	virtual ~linked_page_queue() {
		if (local_allocator)
			delete allocator;
	}

	linked_page_queue(obj_allocator<linked_obj> *allocator) {
		_size = 0;
		this->allocator = allocator;
		local_allocator = false;
	}

	linked_page_queue() {
		_size = 0;
		allocator = new obj_allocator<linked_obj>(PAGE_SIZE);
		local_allocator = true;
	}

	virtual iterator begin() {
		return iterator(&head, this);
	}

	virtual linked_obj *const push_back(frame *pg) {
		linked_obj *obj = allocator->alloc_obj();
		assert(obj);
		obj->set_payload(pg);
		head.add_back(obj);
		_size++;
		return obj;
	}

	virtual void pop_front() {
		if (size() == 0)
			return;
		assert(size() > 0);

		linked_obj *obj = head.front();
		obj->remove_from_list();
		_size--;
		allocator->free(obj);
	}

	bool empty() const {
		return head.is_empty();
	}

	frame *front() {
		if (size() <= 0)
			return NULL;
		return (frame *) head.front()->get_payload();
	}

	frame *back() {
		if (size() <= 0)
			return NULL;
		return (frame *) head.back()->get_payload();
	}

	int size() const {
		return _size;
	}

	void merge(linked_page_queue *list) {
		if (list->empty())
			return;

		linked_obj *list_begin = list->head.front();
		linked_obj *list_end = list->head.back();
		linked_obj *this_end = this->head.back();
		
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

	// for test
	void print() const {
		linked_obj *f = head.front();
		printf("queue size: %d\n", _size);
		while (f != &head) {
			printf("%ld\t", ((frame *) f->get_payload())->get_offset());
			f = f->front();
		}
		printf("\n");
	}

	void remove(int idx) {
		int i = 0;
		for (iterator it = begin(); it.has_next(); i++) {
			it.next();
			if (idx == i)
				it.remove();
		}
	}
};

#endif
