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
	rand_buf(int buf_size, int entry_size): buf_offset(buf_size / entry_size, entry_size) {
		this->entry_size = entry_size;
		num_entries = buf_size / entry_size;
		printf("there are %d entries in the rand buffer\n", num_entries);
		used_entries = 0;
		buf = (char *) numa_alloc_local(buf_size);
		marks = (char *) numa_alloc_local(num_entries);
		memset(marks, 0, num_entries);

		if (buf == NULL){
			fprintf(stderr, "can't allocate buffer\n");
			exit(1);
		}
		/* trigger page faults and bring pages to memory. */
		for (int i = 0; i < buf_size / PAGE_SIZE; i++)
			buf[i * PAGE_SIZE] = 0;

		current = 0;
	}

	~rand_buf() {
		free(buf);
	}

	bool is_full() {
		return used_entries == num_entries;
	}

	void free_entry(char *buf) {
		int off = (buf - this->buf) / entry_size;
		if (marks[off] == 0)
			printf("free %p error\n", buf);
		assert(marks[off]);
		marks[off] = 0;
		used_entries--;
//		printf("free %p\n", buf);
	}

	char *next_entry() {
		int off;
		int num = 0;
		while(1) {
			off = buf_offset.get_offset(current);
			current = (current + 1) % num_entries;;
			if (!marks[off / entry_size])
				break;
			num++;
		}
//		if (num > 1)
//			printf("takes %d iterations\n", num);
		used_entries++;
		marks[off / entry_size] = 1;
//		printf("allocate %p\n", &buf[off]);
		return &buf[off];
	}

	int get_entry_size() {
		return entry_size;
	}

	int get_num_entries() {
		return num_entries;
	}
};

#endif
