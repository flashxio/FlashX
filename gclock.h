#ifndef __GCLOCK_H__
#define __GCLOCK_H__

#include "concurrency.h"
#include "cache.h"

class frame: public thread_safe_page
{
	/* equivalent to the number of hits */
	atomic_integer wcount;
	/* equivalent to the number of references */
	atomic_integer pinning;
	frame *prev, *next;

public:
	frame(): thread_safe_page() {
		prev = next = this;
	}

	frame(off_t offset, char *data): thread_safe_page(offset, data) {
		prev = next = this;
	}

	void *volatileGetValue() {
		/*
		 * maybe a memory fence here,
		 * but it's not needed in the x86 machine.
		 */
		return get_data();
	}

	bool CASValue(void *expect,void *update) {
		return __sync_bool_compare_and_swap(&data, expect, update);
	}

	int getWC() {
		return wcount.get();
	}

	void incrWC(int num = 1) {
		wcount.inc(num);
	}

	int decrWC(int num = 1) {
		return wcount.dec(num);
	}

	bool tryEvict() {
		return pinning.CAS(0, -1);
	}

	void evictUnshared() {
		pinning.CAS(1, -1);
	}

	int pinCount() {
		return pinning.get();
	}

	bool pin() {
		int x;
		do {
			x = pinning.get();
			if (x <= -1)
				return false;
		} while (!pinning.CAS(x, x + 1));
		return true;
	}

	void dec_ref() {
		unpin();
	}

	void unpin() {
		pinning.dec(1);
	}

	/* for linked pages */
	void add_front(frame *pg) {
		frame *next = this->next;
		pg->next = next;
		pg->prev = this;
		this->next = pg;
		next->prev = pg;
	}

	void add_back(frame *pg) {
		frame *prev = this->prev;
		pg->next = this;
		pg->prev = prev;
		this->prev = pg;
		prev->next = pg;
	}

	void remove_from_list() {
		frame *prev = this->prev;
		frame *next = this->next;
		prev->next = next;
		next->prev = prev;
		this->next = this;
		this->prev = this;
	}

	bool is_empty() {
		return this->next == this;
	}

	frame *front() {
		return next;
	}

	frame *back() {
		return prev;
	}

};

class gclock_buffer
{
public:
	/**
	 * Add a frame to the buffer, and evict one if the buffer is full.
	 */
	virtual frame *add(frame *entry) = 0;
};

/**
 * This is a lock-free GClock implementation.
 */
class LF_gclock_buffer: public gclock_buffer
{
	atomic_array<frame *> pool;
	atomic_unsigned_integer free;
	atomic_unsigned_integer clock_hand;
	unsigned int size;

	frame *swap(frame *entry);
public:
	LF_gclock_buffer(unsigned int size): pool(size), free(size) {
		this->size = size;
	}

	virtual frame *add(frame *entry);

	void moveClockHand(int curr, int start) {
		int delta;
		if (curr < start)
			delta = curr + size - start + 1;
		else
			delta = curr - start + 1;
		clock_hand.inc(delta);
	}

	void print() {
		printf("**************************\n");
		printf("there are %d frames in the lf_buffer\n", size - free.get());
		for (unsigned i = 0; i < size - free.get(); i++) {
			printf("\toffset: %ld, hits: %d\n",
					pool.get(i)->get_offset(), pool.get(i)->getWC());
		}
	}
};

/**
 * This is a GClock implementation. It can be used in the following way:
 * each core has a buffer attached to it and each thread accesses its local
 * buffer (the one belongs to the core where the thread is running).
 * Each core has only one thread attached to it, so only one thread can
 * access the gclock buffer.
 *
 * The GClock implementation is to solve two problems:
 * 1. GClock tends to invalidate cache lines rather frequently because
 * it needs to decrease the cache hit counter when it scans the buffer.
 *
 * Solution: use a scan counter. Instead of decreasing the hit counter,
 * increase the scan counter whenever the entire buffer is increased. We need to
 * decrease the hit counter eventually, but we can decrease the number of 
 * cache invalidations significantly.
 *
 * 2. GClock may need to scan many pages in order to evict a page. We need to 
 * roughly sort the pages, so we can avoid scan unnecessary pages.
 *
 * Solution: TODO We can decrease the number of cache misses.
 */
class enhanced_gclock_buffer: public gclock_buffer
{
	int scan_nrounds;	// how many times has the buffer been scanned
	bool start_dec;		// start to decrease the hit count

	frame **pool;
	unsigned free;
	unsigned clock_hand;
	unsigned size;

	frame *swap(frame *entry);
public: 
	enhanced_gclock_buffer(int size): gclock_buffer() {
		pool = (frame **) numa_alloc_local(sizeof(frame *) * size);
		memset(pool, 0, sizeof(frame *) * size);
		free = size;
		clock_hand = 0;
		this->size = size;

		start_dec = false;
		scan_nrounds = 0;
	}

	~enhanced_gclock_buffer() {
		numa_free(pool, sizeof(frame *) * size);
	}

	frame *add(frame *entry);

	void print() {
		printf("**************************\n");
		printf("there are %d frames in the buffer\n", size - free);
		for (unsigned i = 0; i < size - free; i++) {
			printf("\toffset: %ld, hits: %d\n",
					pool[i]->get_offset(), pool[i]->getWC());
		}
	}
};

#endif
