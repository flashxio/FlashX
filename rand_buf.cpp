#include "rand_buf.h"
#include "container.cpp"

rand_buf::rand_buf(int buf_size, int entry_size_): entry_size(
		entry_size_), num_entries(buf_size / entry_size_), free_refs(
			"free_refs", buf_size / entry_size_ + 1) {
	rand_permute buf_offset(num_entries, entry_size_);
	for (int i = 0; i < num_entries; i++) {
		int offset = (int) buf_offset.get_offset(i);
		free_refs.push_back(offset);
	}

	printf("there are %d entries in the rand buffer\n", num_entries);
	buf = (char *) numa_alloc_local(buf_size);

	if (buf == NULL){
		fprintf(stderr, "can't allocate buffer\n");
		exit(1);
	}
	/* trigger page faults and bring pages to memory. */
	for (int i = 0; i < buf_size / PAGE_SIZE; i++)
		buf[i * PAGE_SIZE] = 0;
}

char *rand_buf::next_entry() {
	int off = free_refs.pop_front();
	char *ret = &buf[off];
	return ret;
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
}
