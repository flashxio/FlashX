#ifndef __RAND_BUF_H__
#define __RAND_BUF_H__

#include "workload.h"
#include "container.h"

class rand_buf
{
	/* where the data read from the disk is stored */
	char *buf;
	const int entry_size;
	const int num_entries;
	blocking_FIFO_queue<int> free_refs;

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
};

#endif
