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
	thread_safe_FIFO_queue<off_t> free_refs;
	pthread_spinlock_t lock;

	// The buffers pre-allocated to serve allocation requests
	// from the local threads.
	pthread_key_t local_buf_key;
	// The buffers freed in the local threads, which hasn't been
	// added the main buffer.
	pthread_key_t local_free_key;

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
