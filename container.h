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
	int size;
	T *buf;			// a circular buffer to keep pages.
	long start;
	long end;
	bool resizable;

public:
	fifo_queue(int size, bool resizable = false) {
		this->size = size;
		buf = new T[size];
		start = 0;
		end = 0;
		this->resizable = resizable;
	}

	~fifo_queue() {
		delete [] buf;
	}

	bool expand_queue(int new_size);

	virtual T pop_front() {
		assert(start < end);
		T ret = buf[start % size];
		start++;
		return ret;
	}

	virtual void push_back(T &v) {
		assert(end - start < size);
		buf[end % size] = v;
		end++;
	}

	virtual int fetch(T *entries, int num) {
		int num_fetches = 0;
		while (!is_empty() && num_fetches < num) {
			entries[num_fetches++] = buf[start % size];
			start++;
		}
		return num_fetches;
	}

	virtual int add(fifo_queue<T> *queue) {
		int idx = (int) (end % size);
		int length = min(size - idx, get_num_remaining());
		int num_added = 0;
		int num = queue->fetch(buf + idx, length);
		end += num;
		num_added += num;
		// If we fetch fewer entries than we ask for, it means we have fetched
		// all entries in the queue.
		// or the current queue is full.
		if (num < length || get_num_remaining() == 0)
			return num_added;
		assert(end % size == 0);
		length = get_num_remaining();
		num = queue->fetch(buf, length);
		end += num;
		num_added += num;
		return num_added;
	}

	virtual int add(T *entries, int num) {
		int num_pushes = 0;
		while (!is_full() && num_pushes < num) {
			buf[end % size] = entries[num_pushes++];
			end++;
		}
		return num_pushes;
	}

	int get_num_remaining() {
		return size - fifo_queue<T>::get_num_entries();
	}

	virtual int get_num_entries() {
		return (int) (end - start);
	}

	int get_size() const {
		return size;
	}

	virtual bool is_full() {
		return end - start >= size;
	}

	virtual bool is_empty() {
		return start >= end;
	}
};

/**
 * this is a thread-safe FIFO queue.
 * It supports bulk operations.
 */
template<class T>
class thread_safe_FIFO_queue: public fifo_queue<T>
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
	/* The location where we have completed fetching. */
	volatile long fetched_offset;

	/* 
	 * lock is still needed because we need to check whether the buffer
	 * has entries or has space.
	 */
	pthread_spinlock_t _lock;

	int get_virtual_num_entries() const {
		return (int) (alloc_offset - fetched_offset);
	}

	int get_actual_num_entries() const {
		return (int) (add_offset - fetch_offset);
	}

	int get_remaining_space() const {
		return (int) (capacity - get_virtual_num_entries());
	}

public:
	thread_safe_FIFO_queue(int size): fifo_queue<T>(size), capacity(size) {
		// TODO I don't need to allocate this buffer.
		buf = new T[size];
		alloc_offset = 0;
		add_offset = 0;
		fetch_offset = 0;
		fetched_offset = 0;
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

	virtual ~thread_safe_FIFO_queue() {
		pthread_spin_destroy(&_lock);
		delete [] buf;
	}

	virtual int fetch(T *entries, int num);

	virtual int add(T *entries, int num);
	virtual int add(fifo_queue<T> *queue) {
		assert(0);
		return 0;
	}

	/**
	 * It guarantees to be able to add n entries to the queue.
	 * If there isn't enough space left, it will increase the capacity
	 * of the queue.
	 */
	virtual void addByForce(T *entries, int num);

	// TODO I should return reference.
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
		int num_entries = get_actual_num_entries();
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
class blocking_FIFO_queue: public fifo_queue<T>
{
	pthread_cond_t cond;
	pthread_mutex_t mutex;
	int num_empty;
	int num_full;
	int max_size;

	std::string name;
public:
	blocking_FIFO_queue(const std::string name, int init_size,
			int max_size): fifo_queue<T>(init_size, max_size > init_size) {
		assert(init_size <= max_size);
		pthread_mutex_init(&mutex, NULL);
		pthread_cond_init(&cond, NULL);
		this->name = name;
		num_empty = 0;
		num_full = 0;
		this->max_size = max_size;
	}

	virtual int fetch(T *entries, int num);
	int non_blocking_fetch(T *entries, int num);

	virtual int add(T *entries, int num);

	virtual int add(fifo_queue<T> *queue);
	// Add at least `min_added' elements or all elements in `queue' are added.
	int add_partial(fifo_queue<T> *queue, int min_added = 1);
	int non_blocking_add(fifo_queue<T> *queue);

	int get_num_empty() const {
		return num_empty;
	}

	int get_num_full() const {
		return num_full;
	}

	int get_num_entries() {
		pthread_mutex_lock(&mutex);
		int ret = fifo_queue<T>::get_num_entries();
		pthread_mutex_unlock(&mutex);
		return ret;
	}
};

#endif
