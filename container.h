#ifndef __CONTAINER_H__
#define __CONTAINER_H__

#include <stdio.h>
#include <numa.h>
#include <assert.h>
#include <pthread.h>

#include <string>

#include "common.h"

/**
 * The elements in the queue stored in the same piece of memory
 * as the queue metadata. The size of the queue is defined 
 * during the compile time.
 */
template<class T, int SIZE>
class embedded_queue
{
	unsigned short start;
	unsigned short num;
	/* the size of the buffer is specified by SIZE. */
	T buf[SIZE];
public:
	embedded_queue() {
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

/**
 * this is a first-in-first-out queue.
 * However, the location of an entry in the queue never changes.
 */
template<class T>
class fifo_queue
{
	long size;			// the number of pages that can be buffered
	volatile unsigned long idx;		// to the point where we can evict a page in the buffer
protected:
	T *buf;			// a circular buffer to keep pages.
public:
	fifo_queue(long size) {
		idx = 0;
		this->size = size;
		buf = new T[size];
	}

	~fifo_queue() {
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
 * this is a thread-safe FIFO queue.
 * It supports bulk operations.
 */
template<class T>
class thread_safe_FIFO_queue
{
	volatile T *buf;
	const int capacity;			// capacity of the buffer

	/**
	 * The buffer virtually has infinite space. The three offsets
	 * shows the location in the buffer of infinite size. 
	 * We can easy convert the three offsets to the real location 
	 * in the buffer with modulo operation.
	 */

	/* The location, up to where the space in the buffer has been allocated. */
	volatile long alloc_offset;
	/* The location where new entries have been added to the buffer. */
	volatile long add_offset;
	/* The location where we can fetch entries in the buffer. */
	volatile long fetch_offset;

	/* 
	 * lock is still needed because we need to check whether the buffer
	 * has entries or has space.
	 */
	pthread_spinlock_t _lock;

	int get_virtual_num_entries() const {
		return (int) (alloc_offset - fetch_offset);
	}

	int get_actual_num_entries() const {
		return (int) (add_offset - fetch_offset);
	}

	int get_remaining_space() const {
		return (int) (capacity - get_virtual_num_entries());
	}

public:
	thread_safe_FIFO_queue(int size): capacity(size) {
		buf = new T[size];
		alloc_offset = 0;
		add_offset = 0;
		fetch_offset = 0;
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

	virtual ~thread_safe_FIFO_queue() {
		pthread_spin_destroy(&_lock);
		delete [] buf;
	}

	virtual int fetch(T *entries, int num);

	virtual int add(T *entries, int num);

	T pop_front() {
		T entry;
		int num = fetch(&entry, 1);
		assert(num == 1);
		return entry;
	}

	void push_back(T &entry) {
		while (add(&entry, 1) == 0) {
		}
	}

	/**
	 * Get the existing entries from the queue.
	 * It locks the buffer, so don't over use it.
	 */
	int get_num_entries() {
		pthread_spin_lock(&_lock);
		int num_entries = (int) (add_offset - fetch_offset);
		pthread_spin_unlock(&_lock);
		return num_entries;
	}

	bool is_full() {
		return get_virtual_num_entries() == capacity;
	}

	bool is_empty() {
		return get_actual_num_entries() == 0;
	}
};

/**
 * This FIFO queue can block the thread if
 * a thread wants to add more entries when the queue is full;
 * or
 * a thread wants to fetch more entries when the queue is empty.
 */
template<class T>
class blocking_FIFO_queue: public thread_safe_FIFO_queue<T>
{
	/* when the queue becomes empty */
	pthread_cond_t empty_cond;
	pthread_mutex_t empty_mutex;

	/* when the queue becomes full */
	pthread_cond_t full_cond;
	pthread_mutex_t full_mutex;

	std::string name;
public:
	blocking_FIFO_queue(const std::string name, int size): thread_safe_FIFO_queue<T>(size) {
		pthread_mutex_init(&empty_mutex, NULL);
		pthread_cond_init(&empty_cond, NULL);
		pthread_mutex_init(&full_mutex, NULL);
		pthread_cond_init(&full_cond, NULL);
		this->name = name;
	}

	virtual int fetch(T *entries, int num);

	virtual int add(T *entries, int num);
};

#endif
