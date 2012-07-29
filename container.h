#ifndef __CONTAINER_H__
#define __CONTAINER_H__

#include <stdio.h>

#include "common.h"

template<class T>
class dynamic_queue
{
	int start;
	int num;
	T *buf;
	int size;
public:
	dynamic_queue(int size) {
		buf = (T *) numa_alloc_local(sizeof(T) * size);
		this->size = size;
		start = 0;
		num = 0;
	}

	void push_back(T v) {
		assert(num < size);
		buf[(start + num) % size] = v;
		num++;
	}

	void pop_front() {
		assert(num > 0);
		start = (start + 1) % size;
		num--;
	}

	bool is_empty() {
		return num == 0;
	}

	bool is_full() {
		return num == size;
	}

	T &back() {
		assert(num > 0);
		return buf[(start + num - 1) % size];
	}

	T &front() {
		assert(num > 0);
		return buf[start];
	}
};

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
class bulk_queue
{
	T *buf;
	volatile int size;
	int start;
	int num_entries;
	pthread_spinlock_t _lock;
public:
	bulk_queue(int size) {
		buf = new T[size];
		this->size = size;
		start = 0;
		num_entries = 0;
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

	virtual ~bulk_queue() {
		pthread_spin_destroy(&_lock);
		delete [] buf;
	}

	virtual int fetch(T *entries, int num);

	virtual int add(T *entries, int num);

	int get_num_entries() {
		return num_entries;
	}

	bool is_full() {
		return num_entries == size;
	}

	bool is_empty() {
		return num_entries == 0;
	}
};

template<class T>
int bulk_queue<T>::fetch(T *entries, int num) {
	pthread_spin_lock(&_lock);
	int n = min(num, num_entries);
	for (int i = 0; i < n; i++) {
		entries[i] = buf[(start + i) % this->size];
	}
	start = (start + n) % this->size;
	num_entries -= n;
	pthread_spin_unlock(&_lock);
	return n;
}

/**
 * this is non-blocking. 
 * It adds entries to the queue as much as possible,
 * and returns the number of entries that have been
 * added.
 */
template<class T>
int bulk_queue<T>::add(T *entries, int num) {
	pthread_spin_lock(&_lock);
	int n = min(num, this->size - num_entries);
	int end = (start + num_entries) % this->size;
	for (int i = 0; i < n; i++) {
		buf[(end + i) % this->size] = entries[i];
	}
	num_entries += n;
	pthread_spin_unlock(&_lock);
	return n;
}

#ifdef USE_SHADOW_PAGE

/*
 * remove the idx'th element in the queue.
 * idx is the logical position in the queue,
 * instead of the physical index in the buffer.
 */
template<class T, int SIZE>
void embedded_queue<T, SIZE>::remove(int idx) {
	assert(idx < num);
	/* the first element in the queue. */
	if (idx == 0) {
		pop_front();
	}
	/* the last element in the queue. */
	else if (idx == num - 1){
		num--;
	}
	/*
	 * in the middle.
	 * now we need to move data.
	 */
	else {
		T tmp[num];
		T *p = tmp;
		/* if the end of the queue is physically behind the start */
		if (start + num <= SIZE) {
			/* copy elements in front of the removed element. */
			memcpy(p, &buf[start], sizeof(T) * idx);
			p += idx;
			/* copy elements behind the removed element. */
			memcpy(p, &buf[start + idx + 1], sizeof(T) * (num - idx - 1));
		}
		/* 
		 * the removed element is between the first element
		 * and the end of the buffer.
		 */
		else if (idx + start < SIZE) {
			/* copy elements in front of the removed element. */
			memcpy(p, &buf[start], sizeof(T) * idx);
			p += idx;
			/*
			 * copy elements behind the removed element
			 * and before the end of the buffer.
			 */
			memcpy(p, &buf[start + idx + 1], sizeof(T) * (SIZE - start - idx - 1));
			p += (SIZE - start - idx - 1);
			/* copy the remaining elements in the beginning of the buffer. */
			memcpy(p, buf, sizeof(T) * (num - (SIZE - start)));
		}
		/*
		 * the removed element is between the beginning of the buffer
		 * and the last element.
		 */
		else {
			/* copy elements between the first element and the end of the buffer. */
			memcpy(p, &buf[start], sizeof(T) * (SIZE - start));
			p += (SIZE - start);
			/* copy elements between the beginning of the buffer and the removed element. */
			idx = (idx + start) % SIZE;
			memcpy(p, buf, sizeof(T) * idx);
			p += idx;
			/* copy elements after the removed element and before the last element */
			memcpy(p, &buf[idx + 1], sizeof(T) * ((start + num) % SIZE - idx - 1));
		}
		memcpy(buf, tmp, sizeof(T) * (num - 1));
		start = 0;
		num--;
	}
}    

#endif

#endif
