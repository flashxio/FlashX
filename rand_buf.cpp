#include "rand_buf.h"

rand_buf::rand_buf(int buf_size, int entry_size,
		int nodeid): free_refs(buf_size / entry_size)
#ifdef MEMCHECK
	  , allocator(entry_size)
#endif
{
	num_entries = buf_size / entry_size;
	rand_permute buf_offset(num_entries, entry_size, 0);
	for (int i = 0; i < num_entries; i++) {
		off_t off = buf_offset.get_offset(i);
		free_refs.push_back(off);
	}

	this->entry_size = entry_size;
	printf("there are %d entries in the rand buffer\n", num_entries);
	if (nodeid >= 0)
		buf = (char *) numa_alloc_onnode(buf_size, nodeid);
	else
		buf = (char *) numa_alloc_local(buf_size);
	marks = (char *) numa_alloc_local(num_entries);
	memset(marks, 0, num_entries);
	printf("%ld: rand_buf start %p, end %p\n", pthread_self(), buf, buf + buf_size);

	if (buf == NULL){
		fprintf(stderr, "can't allocate buffer\n");
		exit(1);
	}
	/* trigger page faults and bring pages to memory. */
	for (int i = 0; i < buf_size / PAGE_SIZE; i++)
		buf[i * PAGE_SIZE] = 0;

	current = 0;
	// TODO I shouldn't use spin lock here.
	pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
}

char *rand_buf::next_entry(int size) {
#ifdef MEMCHECK
	return (char *) allocator.alloc(size);
#else
	if (size > entry_size)
		return (char *) valloc(size);
	pthread_spin_lock(&lock);
	assert(!free_refs.is_empty());
	int off = free_refs.pop_front();
	assert(marks[off / entry_size] == 0);
	marks[off / entry_size] = 1;
	char *ret = &buf[off];
	pthread_spin_unlock(&lock);
	return ret;
#endif
}

void rand_buf::free_entry(char *buf) {
#ifdef MEMCHECK
	allocator.dealloc(buf);
#else
	pthread_spin_lock(&lock);
	int buf_size = num_entries * entry_size;
	if (!((long) buf >= (long) this->buf
			&& (long) buf < (long) this->buf + buf_size)) {
		pthread_spin_unlock(&lock);
		free(buf);
		return;
	}
	off_t off = buf - this->buf;
	free_refs.push_back(off);
	off /= entry_size;
	if (marks[off] == 0)
		printf("%ld: free %p error\n", pthread_self(), buf);
	assert(marks[off]);
	marks[off] = 0;
	pthread_spin_unlock(&lock);
#endif
}
