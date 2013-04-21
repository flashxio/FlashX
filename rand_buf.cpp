#include "rand_buf.h"
#include "container.cpp"

rand_buf::rand_buf(int buf_size, int entry_size,
		int nodeid): free_refs(buf_size / entry_size)
#ifdef MEMCHECK
	  , allocator(entry_size)
#endif
{
	num_entries = buf_size / entry_size;
	for (int i = 0; i < num_entries; i++) {
		off_t off = i * entry_size;
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
	pthread_key_create(&local_buf_key, NULL);
	pthread_key_create(&local_free_key, NULL);
}

char *rand_buf::next_entry(int size) {
#ifdef MEMCHECK
	return (char *) allocator.alloc(size);
#else
	if (size > entry_size)
		return (char *) valloc(size);
	fifo_queue<off_t> *local_buf_refs
		= (fifo_queue<off_t> *) pthread_getspecific(local_buf_key);
	if (local_buf_refs == NULL) {
		local_buf_refs = new fifo_queue<off_t>(LOCAL_BUF_SIZE);
		pthread_setspecific(local_buf_key, local_buf_refs);
	}
	if (local_buf_refs->is_empty()) {
		off_t offs[LOCAL_BUF_SIZE];
		int num = free_refs.fetch(offs, LOCAL_BUF_SIZE);
		assert(num > 0);
		int num_added = local_buf_refs->add(offs, num);
		assert(num_added == num);
	}
	off_t off = local_buf_refs->pop_front();
	return &buf[off];
#if 0
	pthread_spin_lock(&lock);
	assert(!free_refs.is_empty());
	int off = free_refs.pop_front();
	assert(marks[off / entry_size] == 0);
	marks[off / entry_size] = 1;
	char *ret = &buf[off];
	pthread_spin_unlock(&lock);
	return ret;
#endif
#endif
}

void rand_buf::free_entry(char *buf) {
#ifdef MEMCHECK
	allocator.dealloc(buf);
#else
	const int buf_size = num_entries * entry_size;
	if (!((long) buf >= (long) this->buf
			&& (long) buf < (long) this->buf + buf_size)) {
		free(buf);
		return;
	}
	fifo_queue<off_t> *local_free_refs
		= (fifo_queue<off_t> *) pthread_getspecific(local_free_key);
	if (local_free_refs == NULL) {
		local_free_refs = new fifo_queue<off_t>(LOCAL_BUF_SIZE);
		pthread_setspecific(local_free_key, local_free_refs);
	}
	if (local_free_refs->is_full()) {
		off_t offs[LOCAL_BUF_SIZE];
		int num = local_free_refs->fetch(offs, LOCAL_BUF_SIZE);
		int num_added = free_refs.add(offs, num);
		assert(num == num_added);
	}
	off_t off = buf - this->buf;
	local_free_refs->push_back(off);
#if 0
	pthread_spin_lock(&lock);
	off_t off = buf - this->buf;
	free_refs.push_back(off);
	off /= entry_size;
	if (marks[off] == 0)
		printf("%ld: free %p error\n", pthread_self(), buf);
	assert(marks[off]);
	marks[off] = 0;
	pthread_spin_unlock(&lock);
#endif
#endif
}
