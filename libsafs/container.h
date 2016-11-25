#ifndef __CONTAINER_H__
#define __CONTAINER_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <limits.h>

#include <string>
#include <boost/assert.hpp>

#include "common.h"
#include "concurrency.h"

template<class T>
class queue_interface
{
public:
	virtual ~queue_interface() {
	}

	virtual T pop_front() = 0;
	virtual void push_back(T &v) = 0;
	virtual int get_num_entries() = 0;
	virtual bool is_empty() = 0;
};

/*
 * this is a first-in-first-out queue.
 * However, the location of an entry in the queue never changes.
 */
template<class T>
class fifo_queue: public queue_interface<T>
{
	int size_mask;
	T *buf;			// a circular buffer to keep pages.
	bool allocated;	// indicates whether the buffer is allocated by the queue.
	long start;
	long end;
	bool resizable;
	int node_id;

	long loc_in_queue(long idx) const {
		return idx & size_mask;
	}

	T *alloc_buf(int size) {
		return new T[size];
	}

	void free_buf(T *buf) {
		assert(allocated);
		delete [] buf;
	}

public:
	class const_iterator {
		const fifo_queue<T> *q;
		long idx;
	public:
		const_iterator(const fifo_queue<T> *q) {
			this->q = q;
			this->idx = q->start;
		}

		const T &operator*() const {
			long real_idx = q->loc_in_queue(idx);
			return q->buf[real_idx];
		}

		const_iterator &operator++() {
			idx++;
			return *this;
		}

		bool operator==(const const_iterator &it) const {
			return this->idx == it.idx;
		}

		bool operator!=(const const_iterator &it) const {
			return this->idx != it.idx;
		}

		const_iterator &operator+=(int num) {
			assert(num >= 0);
			idx += num;
			assert(idx <= q->end);
			return *this;
		}
	};

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
		return new fifo_queue<T>(node_id, size, resizable);
	}

	static void destroy(fifo_queue<T> *q) {
		delete q;
	}

	const_iterator get_begin() const {
		return const_iterator(this);
	}

	const_iterator get_end() const {
		const_iterator it(this);
		assert(end >= start);
		it += (end - start);
		return it;
	}

	bool expand_queue(int new_size);

	/*
	 * Get the first element in the queue.
	 */
	T &front() {
		assert(!is_empty());
		return buf[loc_in_queue(start)];
	}

	/*
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
		while (!fifo_queue<T>::is_empty() && num_fetches < num) {
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
		while (!fifo_queue<T>::is_full() && num_pushes < num) {
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

template<class T>
class thread_safe_FIFO_queue;

/*
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
	spin_lock _lock;
	int max_size;
	std::string name;

public:
	thread_safe_FIFO_queue(const std::string &name, int node_id,
			int size): fifo_queue<T>(node_id, size, false) {
		this->name = name;
		this->max_size = size;
	}

	thread_safe_FIFO_queue(const std::string &name, int node_id, int init_size,
			int max_size): fifo_queue<T>(node_id, init_size,
				max_size > init_size) {
		this->name = name;
		this->max_size = max_size;
	}

	virtual ~thread_safe_FIFO_queue() {
	}

	static thread_safe_FIFO_queue<T> *create(const std::string &name,
			int node_id, int size) {
		return new thread_safe_FIFO_queue<T>(name, node_id, size);
	}

	static void destroy(thread_safe_FIFO_queue<T> *q) {
		delete q;
	}

	virtual int fetch(T *entries, int num) {
		_lock.lock();
		int ret = fifo_queue<T>::fetch(entries, num);
		_lock.unlock();
		return ret;
	}

	virtual int add(T *entries, int num) {
		_lock.lock();
		int ret = fifo_queue<T>::add(entries, num);
		int orig_size = fifo_queue<T>::get_size();
		if (ret < num && orig_size < max_size) {
			int new_size = orig_size;
			int min_required_size = orig_size + num - ret;
			while (new_size < min_required_size && new_size < max_size)
				new_size *= 2;
			fifo_queue<T>::expand_queue(new_size);
			int ret2 = fifo_queue<T>::add(entries + ret, num - ret);
			ret += ret2;
			assert(ret == num);
		}
		_lock.unlock();
		return ret;
	}
	virtual int add(fifo_queue<T> *queue) {
		_lock.lock();
		int ret = fifo_queue<T>::add(queue);
		int orig_size = fifo_queue<T>::get_size();
		if (!queue->is_empty() && orig_size < max_size) {
			int new_size = orig_size;
			int min_required_size = orig_size + queue->get_num_entries();
			while (new_size < min_required_size && new_size < max_size)
				new_size *= 2;
			fifo_queue<T>::expand_queue(new_size);
			int ret2 = fifo_queue<T>::add(queue);
			ret += ret2;
			assert(queue->is_empty());
		}
		_lock.unlock();
		return 0;
	}

	/*
	 * It guarantees to be able to add n entries to the queue.
	 * If there isn't enough space left, it will increase the capacity
	 * of the queue.
	 */
	virtual void addByForce(T *entries, int num) {
		BOOST_VERIFY(add(entries, num) == num);
	}

	// TODO I should return reference.
	T pop_front() {
		T entry;
		BOOST_VERIFY(fetch(&entry, 1) == 1);
		return entry;
	}

	void push_back(T &entry) {
		while (add(&entry, 1) == 0);
	}

	/*
	 * Get the existing entries from the queue.
	 * It locks the buffer, so don't over use it.
	 */
	int get_num_entries() {
		_lock.lock();
		int num_entries = fifo_queue<T>::get_num_entries();
		_lock.unlock();
		return num_entries;
	}

	// TODO these are bugs. They should be protected by locks.
	bool is_full() {
		_lock.lock();
		bool ret = fifo_queue<T>::is_full();
		_lock.unlock();
		return ret;
	}

	bool is_empty() {
		_lock.lock();
		bool ret = fifo_queue<T>::is_empty();
		_lock.unlock();
		return ret;
	}

	const std::string &get_name() const {
		return name;
	}
};

/*
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
		return new blocking_FIFO_queue<T>(node_id, name, init_size, max_size);
	}

	static void destroy(blocking_FIFO_queue<T> *q) {
		delete q;
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
	int add(T *entries, int num, bool blocking, bool interruptible) {
		// TODO
		return -1;
	}

	/*
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

/*
 * This is used to allocate an array on the stack.
 * Some code needs to allocate a large array sometimes, but most of time,
 * it only needs a small array.
 * This wrapper class is to avoid allocating a large array on the stack,
 * while being efficient most of times.
 */
template<class T, int size = 32>
class stack_array
{
	T buf[size];
	T *real_buf;
	int capacity;
public:
	stack_array(int capacity) {
		if (capacity <= size) {
			this->capacity = size;
			real_buf = buf;
		}
		else {
			this->capacity = capacity;
			real_buf = new T[capacity];
		}
	}

	~stack_array() {
		if (real_buf != buf)
			delete [] real_buf;
	}

	T &operator[](int idx) {
		assert(idx < capacity);
		return real_buf[idx];
	}

	T *data() {
		return real_buf;
	}

	const T *data() const {
		return real_buf;
	}
};

template<class T, int size = 32>
class embedded_array
{
	T buf[size];
	T *real_buf;
	int capacity;

	void assign(embedded_array<T, size> &arr) {
		this->clear();

		if (arr.real_buf == arr.buf) {
			for (int i = 0; i < capacity; i++) {
				this->buf[i] = arr.buf[i];
				arr.buf[i] = T();
			}
			this->real_buf = this->buf;
			this->capacity = size;
		}
		else {
			this->real_buf = arr.real_buf;
			this->capacity = arr.capacity;
			arr.real_buf = arr.buf;
			arr.capacity = size;
		}
	}
public:
	embedded_array() {
		real_buf = buf;
		capacity = size;
	}

	embedded_array(embedded_array<T, size> &arr) {
		assign(arr);
	}

	embedded_array &operator=(embedded_array<T, size> &arr) {
		assign(arr);
		return *this;
	}

	~embedded_array() {
		if (real_buf != buf)
			delete [] real_buf;
	}

	T &operator[](int idx) {
		assert(idx < capacity);
		return real_buf[idx];
	}

	const T &operator[](int idx) const {
		assert(idx < capacity);
		return real_buf[idx];
	}

	// This can only increase the capacity of the array.
	void resize(int new_size) {
		if (new_size <= capacity)
			return;

		if (real_buf == buf) {
			real_buf = new T[new_size];
			memcpy(real_buf, buf, sizeof(buf));
		}
		else {
			T *tmp = new T[new_size];
			memcpy(tmp, real_buf, capacity * sizeof(T));
			delete [] real_buf;
			real_buf = tmp;
		}
		capacity = new_size;
	}

	void clear() {
		if (this->real_buf != this->buf) {
			delete [] this->real_buf;
			this->real_buf = this->buf;
			this->capacity = size;
		}
	}

	int get_capacity() const {
		return capacity;
	}

	T *data() {
		return real_buf;
	}

	const T *data() const {
		return real_buf;
	}
};

template<class T>
bool fifo_queue<T>::expand_queue(int new_size)
{
	int log_size = (int) ceil(log2(new_size));
	new_size = 1 << log_size;
	assert(resizable && get_size() < new_size);
	assert(allocated);

	// Allocate new memory for the array and initialize it.
	T *tmp = alloc_buf(new_size);

	// Copy the old array to the new one.
	int num = fifo_queue<T>::get_num_entries();
	for (int i = 0; i < num; i++) {
		tmp[i] = buf[loc_in_queue(start + i)];
	}

	// Destroy the old array.
	free_buf(buf);

	buf = tmp;
	size_mask = new_size - 1;
	start = 0;
	end = num;
	return true;
}

template<class T>
int blocking_FIFO_queue<T>::add_partial(fifo_queue<T> *queue, int min_added)
{
	int num_added = 0;
	while (!queue->is_empty()) {
		int num = queue->get_num_entries();
		pthread_mutex_lock(&mutex);
		bool empty = this->is_empty();
		if (this->get_size() - fifo_queue<T>::get_num_entries() < num
				&& this->get_size() < max_size) {
			int new_size = this->get_size() * 2;
			new_size = max(new_size, fifo_queue<T>::get_num_entries() + num);
			new_size = min(new_size, max_size);
#ifdef DEBUG
			printf("try to expand queue %s to %d\n", name.c_str(), new_size);
#endif
			BOOST_VERIFY(fifo_queue<T>::expand_queue(new_size));
		}
		int ret = fifo_queue<T>::add(queue);
		num_added += ret;
		/* signal the thread of reading disk to wake up. */
		if (empty)
			pthread_cond_broadcast(&cond);

		/* We only block the thread when it doesn't send enough data. */
		if (num_added < min_added) {
			while (this->is_full() && !queue->is_empty()) {
				num_full++;
				pthread_cond_wait(&cond, &mutex);
			}
			pthread_mutex_unlock(&mutex);
		}
		else {
			pthread_mutex_unlock(&mutex);
			break;
		}
	}
	return num_added;
}

template<class T>
int blocking_FIFO_queue<T>::fetch(T *entries, int num, bool blocking,
		bool interruptible)
{
	/* we have to wait for coming requests. */
	pthread_mutex_lock(&mutex);
	if (blocking) {
		while(this->is_empty()) {
			num_empty++;
			if (interruptible && interrupted) {
				// We need to reset the interrupt signal.
				interrupted = false;
				break;
			}
			pthread_cond_wait(&cond, &mutex);
		}
	}
	bool full = this->is_full();
	int ret = fifo_queue<T>::fetch(entries, num);
	pthread_mutex_unlock(&mutex);

	/* wake up all threads to send more requests */
	if (full)
		pthread_cond_broadcast(&cond);

	return ret;
}

#endif
