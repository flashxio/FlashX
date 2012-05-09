#include "rand_buf.h"

rand_buf::rand_buf(int buf_size, int entry_size): free_refs(buf_size / entry_size) {
	num_entries = buf_size / entry_size;
	rand_permute buf_offset(num_entries, entry_size);
	for (int i = 0; i < num_entries; i++)
		free_refs.push_back(buf_offset.get_offset(i));

	this->entry_size = entry_size;
	printf("there are %d entries in the rand buffer\n", num_entries);
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
	assert(!free_refs.is_empty());
	int off = free_refs.front();
	free_refs.pop_front();
	assert(marks[off / entry_size] == 0);
	marks[off / entry_size] = 1;
	return &buf[off];
}

void rand_buf::free_entry(char *buf) {
	int buf_size = num_entries * entry_size;
	if (!((long) buf >= (long) this->buf
			&& (long) buf < (long) this->buf + buf_size))
		printf("fail to free %p\n", buf);
	assert((long) buf >= (long) this->buf
			&& (long) buf < (long) this->buf + buf_size);
	int off = buf - this->buf;
	free_refs.push_back(off);
	off /= entry_size;
	if (marks[off] == 0)
		printf("free %p error\n", buf);
	assert(marks[off]);
	marks[off] = 0;
}
