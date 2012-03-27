#include "rand_buf.h"

rand_buf::rand_buf(int buf_size, int entry_size): buf_offset(buf_size / entry_size, entry_size) {
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

char *rand_buf::next_entry() {
	int off;
	int num = 0;
	while(1) {
		off = buf_offset.get_offset(current);
		current = (current + 1) % num_entries;;
		if (!marks[off / entry_size])
			break;
		num++;
	}
	used_entries++;
	marks[off / entry_size] = 1;
	return &buf[off];
}

void rand_buf::free_entry(char *buf) {
	int off = (buf - this->buf) / entry_size;
	if (marks[off] == 0)
		printf("free %p error\n", buf);
	assert(marks[off]);
	marks[off] = 0;
	used_entries--;
}
