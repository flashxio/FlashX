#ifndef __RAND_BUF_H__
#define __RAND_BUF_H__

#include "container.h"
#include "aligned_allocator.h"

/**
 * The class maintains a set of buffers of a mixed size for a single thread.
 * It works the same as slab_allocator.
 * To reduce overhead, the class isn't thread-safe.
 */
class rand_buf
{
	/* where the data read from the disk is stored */
	char *buf;
	char *marks;
	int entry_size;
	int num_entries;
	fifo_queue<off_t> free_refs;

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
