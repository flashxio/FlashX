#ifndef __CONCURRENCY_H__
#define __CONCURRENCY_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <pthread.h>
#include <numa.h>

template<class T>
class atomic_number
{
	volatile T v;
public:
	atomic_number() {
		v = 0;
	}

	atomic_number(T init) {
		v = init;
	}

	T inc(T by) {
		return __sync_add_and_fetch(&v, by);
	}

	T dec(T by) {
		return __sync_sub_and_fetch(&v, by);
	}

	T get() const {
		return v;
	}

	bool CAS(T expected, T value) {
		return __sync_bool_compare_and_swap(&v, expected, value);
	}
};

class atomic_unsigned_integer: public atomic_number<unsigned int>
{
public:
	atomic_unsigned_integer() {
	}

	atomic_unsigned_integer(unsigned int init): atomic_number<unsigned int>(
			init) {
	}
};

class atomic_integer: public atomic_number<int>
{
public:
	atomic_integer() {
	}

	atomic_integer(int init): atomic_number<int>(init) {
	}
};

class atomic_long: public atomic_number<long>
{
public:
	atomic_long() {
	}

	atomic_long(long init): atomic_number<long>(init) {
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

	void read_lock(unsigned long &count) const {
		/*
		 * If the count is odd, it means another thread is changing
		 * the data structure the lock tries to protect. 
		 */
		do {
			count = this->count;
		} while (count & 1);
	}

	bool read_unlock(unsigned long count) const {
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

template<class T>
class atomic_flags
{
	volatile T flags;
public:
	atomic_flags() {
		flags = 0;
	}

	/**
	 * Set the flag and return the orig flag.
	 */
	bool set_flag(int flag) {
		int orig = __sync_fetch_and_or(&flags, 0x1 << flag);
		return orig & (0x1 << flag);
	}

	/**
	 * Clear the flag and return the orig flag.
	 */
	bool clear_flag(int flag) {
		int orig = __sync_fetch_and_and(&flags, ~(0x1 << flag));
		return orig & (0x1 << flag);
	}

	bool test_flag(int flag) const {
		return flags & (0x1 << flag);
	}

	int get_num_tot_flags() const {
		return sizeof(T) * 8;
	}
};

#endif
