#ifndef __RAND_BUF_H__
#define __RAND_BUF_H__

#include "workload.h"
#include "container.h"
#include "aligned_allocator.h"

class rand_buf
{
	/* where the data read from the disk is stored */
	char *buf;
	char *marks;
	int entry_size;
	int num_entries;
	fifo_queue<int> free_refs;
	pthread_spinlock_t lock;

	int current;
#ifdef MEMCHECK
	aligned_allocator allocator;
#endif
public:
	rand_buf(int buf_size, int entry_size, int nodeid = -1);

	~rand_buf() {
		free(buf);
	}

	bool is_full() {
		return free_refs.is_empty();
	}

	void free_entry(char *buf);

	char *next_entry(int size);

	int get_entry_size() {
		return entry_size;
	}

	int get_num_entries() {
		return num_entries;
	}
};

#endif
