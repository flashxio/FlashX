#ifndef __CONCURRENCY_H__
#define __CONCURRENCY_H__

#include <pthread.h>
#include <numa.h>

class atomic_unsigned_integer
{
	volatile unsigned int v;
public:
	atomic_unsigned_integer() {
		v = 0;
	}

	atomic_unsigned_integer(unsigned int init) {
		v = init;
	}

	unsigned int inc(unsigned int by) {
		return __sync_add_and_fetch(&v, by);
	}

	unsigned int dec(unsigned int by) {
		return __sync_sub_and_fetch(&v, by);
	}

	unsigned int get() const {
		return v;
	}

	bool CAS(unsigned int expected, unsigned int value) {
		return __sync_bool_compare_and_swap(&v, expected, value);
	}
};

class atomic_integer
{
	volatile int v;
public:
	atomic_integer() {
		v = 0;
	}

	atomic_integer(int init) {
		v = init;
	}

	int inc(int by) {
		return __sync_add_and_fetch(&v, by);
	}

	int dec(int by) {
		return __sync_sub_and_fetch(&v, by);
	}

	int get() const {
		return v;
	}

	bool CAS(int expected, int value) {
		return __sync_bool_compare_and_swap(&v, expected, value);
	}
};

template<class T>
class atomic_array
{
	volatile T *arr;
	int size;
public:
	atomic_array(int size) {
		this->size = size;
		arr = (T *) numa_alloc_local(sizeof(T) * size);
		memset((void *) arr, 0, sizeof(T) * size);
	}

	~atomic_array() {
		numa_free((void *) arr, sizeof(T) * size);
	}

	T get(int idx) const {
		return arr[idx];
	}

	bool CAS(int idx, T expected, T value) {
		return __sync_bool_compare_and_swap(&arr[idx], expected, value);
	}
};

class seq_lock
{
	volatile unsigned long count;
	pthread_spinlock_t lock;
public:
	seq_lock() {
		count = 0;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	~seq_lock() {
		pthread_spin_destroy(&lock);
	}

	void read_lock(unsigned long &count) {
		/*
		 * If the count is odd, it means another thread is changing
		 * the data structure the lock tries to protect. 
		 */
		do {
			count = this->count;
		} while (count & 1);
	}

	bool read_unlock(unsigned long count) {
		return this->count == count;
	}

	void write_lock() {
		pthread_spin_lock(&lock);
		__sync_fetch_and_add(&count, 1);
	}

	void write_unlock() {
		__sync_fetch_and_add(&count, 1);
		pthread_spin_unlock(&lock);
	}
};

class atomic_flags
{
	volatile int flags;
public:
	atomic_flags() {
		flags = 0;
	}

	void set_flags(int flag) {
		__sync_fetch_and_or(&flags, 0x1 << flag);
	}

	void clear_flags(int flag) {
		__sync_fetch_and_and(&flags, ~(0x1 << flag));
	}

	bool test_flags(int flag) {
		return flags & (0x1 << flag);
	}

	bool test_and_set_flags(int flag) {
		int orig = __sync_fetch_and_or(&flags, 0x1 << flag);
		return orig & (0x1 << flag);
	}
};

#endif
