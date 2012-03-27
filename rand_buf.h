#ifndef __RAND_BUF_H__
#define __RAND_BUF_H__

#include "workload.h"

class rand_buf
{
	/* where the data read from the disk is stored */
	char *buf;
	char *marks;
	/* shows the locations in the array where data has be to stored.*/
	rand_permute buf_offset;
	int entry_size;
	int num_entries;
	int used_entries;

	int current;
public:
	rand_buf(int buf_size, int entry_size);

	~rand_buf() {
		free(buf);
	}

	bool is_full() {
		return used_entries == num_entries;
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
