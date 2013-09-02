#ifndef __CONTAINER_H__
#define __CONTAINER_H__

#include <stdio.h>
#include <math.h>
#include <numa.h>
#include <assert.h>
#include <pthread.h>
#include <limits.h>

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
	int size_mask;
	T *buf;			// a circular buffer to keep pages.
	bool allocated;	// indicates whether the buffer is allocated by the queue.
	long start;
	long end;
	bool resizable;
	int node_id;

	long loc_in_queue(long idx) {
		return idx & size_mask;
	}

	T *alloc_buf(int size) {
		void *addr;
		if (node_id < 0)
			addr = numa_alloc_local(sizeof(T) * size);
		else
			addr = numa_alloc_onnode(sizeof(T) * size, node_id);
		T *buf = (T *) addr;
		for (int i = 0; i < size; i++)
			new(&buf[i]) T();
		return buf;
	}

	void free_buf(T *buf) {
		assert(allocated);
		int size = size_mask + 1;
		for (int i = 0; i < size; i++)
			buf[i].~T();
		numa_free(buf, sizeof(T) * size);
	}

public:
	fifo_queue(T *entries, int num) {
		size_mask = INT_MAX;
		buf = entries;
		allocated = false;
		start = 0;
		end = num;
		resizable = false;
		node_id = -1;
	}

	// the queue has to be 2^n. If it's not, the smallest number of 2^n
	// is used.
	fifo_queue(int node_id, int size, bool resizable = false) {
		int log_size = (int) ceil(log2(size));
		size = 1 << log_size;
		this->size_mask = size - 1;
		this->node_id = node_id;
		buf = alloc_buf(size);
		allocated = true;
		start = 0;
		end = 0;
		this->resizable = resizable;
	}

	virtual ~fifo_queue() {
		if (allocated)
			free_buf(buf);
	}

	static fifo_queue<T> *create(int node_id, int size,
			bool resizable = false) {
		void *addr;
		if (node_id < 0)
			addr = numa_alloc_local(sizeof(fifo_queue<T>));
		else
			addr = numa_alloc_onnode(sizeof(fifo_queue<T>), node_id);
		return new(addr) fifo_queue<T>(node_id, size, resizable);
	}

	static void destroy(fifo_queue<T> *q) {
		q->~fifo_queue();
		numa_free(q, sizeof(*q));
	}

	bool expand_queue(int new_size);

	/**
	 * Get the first element in the queue.
	 */
	T &front() {
		assert(!is_empty());
		return buf[loc_in_queue(start)];
	}

	/**
	 * Get the last element in the queue.
	 */
	T &back() {
		assert(!is_empty());
		return buf[loc_in_queue(end - 1)];
	}

	virtual T pop_front() {
		assert(start < end);
		T ret = buf[loc_in_queue(start)];
		start++;
		return ret;
	}

	virtual void push_back(T &v) {
		assert(end - start < get_size());
		buf[loc_in_queue(end)] = v;
		end++;
	}

	virtual int fetch(T *entries, int num) {
		int num_fetches = 0;
		while (!is_empty() && num_fetches < num) {
			entries[num_fetches++] = buf[loc_in_queue(start)];
			start++;
		}
		return num_fetches;
	}

	virtual int add(fifo_queue<T> *queue) {
		int idx = (int) (loc_in_queue(end));
		int length = min(get_size() - idx, get_num_remaining());
		int num_added = 0;
		int num = queue->fetch(buf + idx, length);
		end += num;
		num_added += num;
		// If we fetch fewer entries than we ask for, it means we have fetched
		// all entries in the queue.
		// or the current queue is full.
		if (num < length || get_num_remaining() == 0)
			return num_added;
		assert(loc_in_queue(end) == 0);
		length = get_num_remaining();
		num = queue->fetch(buf, length);
		end += num;
		num_added += num;
		return num_added;
	}

	virtual int add(T *entries, int num) {
		int num_pushes = 0;
		while (!is_full() && num_pushes < num) {
			buf[loc_in_queue(end)] = entries[num_pushes++];
			end++;
		}
		return num_pushes;
	}

	int get_num_remaining() {
		return get_size() - fifo_queue<T>::get_num_entries();
	}

	virtual int get_num_entries() {
		return (int) (end - start);
	}

	int get_size() const {
		assert(allocated);
		return size_mask + 1;
	}

	virtual bool is_full() {
		return end - start >= get_size();
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
	/* 
	 * lock is still needed because we need to check whether the buffer
	 * has entries or has space.
	 */
	pthread_spinlock_t _lock;

public:
	thread_safe_FIFO_queue(int node_id, int size): fifo_queue<T>(node_id, size) {
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

	virtual ~thread_safe_FIFO_queue() {
		pthread_spin_destroy(&_lock);
	}

	static thread_safe_FIFO_queue<T> *create(int node_id, int size) {
		void *addr;
		if (node_id < 0)
			addr = numa_alloc_local(sizeof(thread_safe_FIFO_queue<T>));
		else
			addr = numa_alloc_onnode(sizeof(thread_safe_FIFO_queue<T>), node_id);
		return new(addr) thread_safe_FIFO_queue<T>(node_id, size);
	}

	static void destroy(thread_safe_FIFO_queue<T> *q) {
		q->~thread_safe_FIFO_queue();
		numa_free(q, sizeof(*q));
	}

	virtual int fetch(T *entries, int num) {
		pthread_spin_lock(&_lock);
		int ret = fifo_queue<T>::fetch(entries, num);
		pthread_spin_unlock(&_lock);
		return ret;
	}

	virtual int add(T *entries, int num) {
		pthread_spin_lock(&_lock);
		int ret = fifo_queue<T>::add(entries, num);
		pthread_spin_unlock(&_lock);
		return ret;
	}
	virtual int add(fifo_queue<T> *queue) {
		assert(0);
		return 0;
	}

	/**
	 * It guarantees to be able to add n entries to the queue.
	 * If there isn't enough space left, it will increase the capacity
	 * of the queue.
	 */
	virtual void addByForce(T *entries, int num) {
		int added = add(entries, num);
		assert(added == num);
		// TODO I should make the queue extensible.
	}

	// TODO I should return reference.
	T pop_front() {
		T entry;
		int num = fetch(&entry, 1);
		assert(num == 1);
		return entry;
	}

	void push_back(T &entry) {
		while (add(&entry, 1) == 0);
	}

	/**
	 * Get the existing entries from the queue.
	 * It locks the buffer, so don't over use it.
	 */
	int get_num_entries() {
		pthread_spin_lock(&_lock);
		int num_entries = fifo_queue<T>::get_num_entries();
		pthread_spin_unlock(&_lock);
		return num_entries;
	}

	// TODO these are bugs. They should be protected by locks.
	bool is_full() {
		bool ret = fifo_queue<T>::is_full();
		return ret;
	}

	bool is_empty() {
		bool ret = fifo_queue<T>::is_empty();
		return ret;
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
	bool interrupted;

	std::string name;
public:
	blocking_FIFO_queue(int node_id, const std::string name, int init_size,
			int max_size): fifo_queue<T>(node_id, init_size, max_size > init_size) {
		assert(init_size <= max_size);
		pthread_mutex_init(&mutex, NULL);
		pthread_cond_init(&cond, NULL);
		this->name = name;
		num_empty = 0;
		num_full = 0;
		this->max_size = max_size;
		interrupted = false;
	}

	static blocking_FIFO_queue<T> *create(int node_id, const std::string name,
			int init_size, int max_size) {
		void *addr;
		if (node_id < 0)
			addr = numa_alloc_local(sizeof(blocking_FIFO_queue<T>));
		else
			addr = numa_alloc_onnode(sizeof(blocking_FIFO_queue<T>), node_id);
		return new(addr) blocking_FIFO_queue<T>(node_id, name, init_size, max_size);
	}

	static void destroy(blocking_FIFO_queue<T> *q) {
		q->~blocking_FIFO_queue();
		numa_free(q, sizeof(*q));
	}

	virtual int fetch(T *entries, int num) {
		return fetch(entries, num, true, false);
	}
	virtual int add(T *entries, int num) {
		fifo_queue<T> tmp(entries, num);
		return add(&tmp);
	}
	virtual int add(fifo_queue<T> *queue) {
		return add_partial(queue, INT_MAX);
	}

	// Add at least `min_added' elements or all elements in `queue' are added.
	int add_partial(fifo_queue<T> *queue, int min_added = 1);
	int non_blocking_add(fifo_queue<T> *queue) {
		return add_partial(queue, 0);
	}
	int non_blocking_add(T *entries, int num) {
		fifo_queue<T> tmp(entries, num);
		return non_blocking_add(&tmp);
	}
	int non_blocking_fetch(T *entries, int num) {
		return fetch(entries, num, false, false);
	}

	int fetch(T *entries, int num, bool blocking, bool interruptible);
	int add(T *entries, int num, bool blocking, bool interruptible);

	/**
	 * This method wakes up the thread that is waiting on the queue
	 * and can be interrupted.
	 * This can only wake up one thread.
	 */
	void wakeup() {
		pthread_mutex_lock(&mutex);
		// We only try to wake up the thread when it's waiting for requests.
		if (this->is_empty()) {
			interrupted = true;
			pthread_cond_signal(&cond);
		}
		pthread_mutex_unlock(&mutex);
	}

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

template class thread_safe_FIFO_queue<long>;

#endif
