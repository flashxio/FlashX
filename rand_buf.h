#ifndef __RAND_BUF_H__
#define __RAND_BUF_H__

#include "workload.h"

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

class rand_buf
{
	/* where the data read from the disk is stored */
	char *buf;
	char *marks;
	int entry_size;
	int num_entries;
	dynamic_queue<int> free_refs;
	pthread_spinlock_t lock;

	int current;
public:
	rand_buf(int buf_size, int entry_size);

	~rand_buf() {
		free(buf);
	}

	bool is_full() {
		return free_refs.is_empty();
	}

	void free_entry(char *buf);

	char *next_entry();

	int get_entry_size() {
		return entry_size;
	}

	int get_num_entries() {
		return num_entries;
	}
};

#endif
