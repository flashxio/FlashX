#ifndef __GCLOCK_H__
#define __GCLOCK_H__

#include "container.h"
#include "concurrency.h"
#include "cache.h"

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

	long tot_nwrites;
	long tot_naccesses;

	frame *swap(frame *entry);
public:
	LF_gclock_buffer(unsigned int size): pool(size), free(size) {
		printf("lock-free glock buffer\n");
		this->size = size;
		tot_nwrites = 0;
		tot_naccesses = 0;
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

	void print(bool stat_only = false) {
		printf("**************************\n");
		printf("there are %d frames in the lf_buffer, clock_hand: %d\n",
				size - free.get(), clock_hand.get());
		printf("there were %ld writes and %ld page accesses\n", tot_nwrites, tot_naccesses);
		if (!stat_only)
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
 */
class enhanced_gclock_buffer: public gclock_buffer
{
	int scan_nrounds;	// how many times has the buffer been scanned
	bool start_dec;		// start to decrease the hit count

	long tot_nwrites;
	long tot_naccesses;
	long tot_nrounds;
	long tot_write_nrounds;

	frame **pool;
	unsigned free;
	unsigned clock_hand;
	unsigned size;

	frame *swap(frame *entry);
public: 
	enhanced_gclock_buffer(int size): gclock_buffer() {
		printf("enhanced glock buffer\n");
		pool = (frame **) numa_alloc_local(sizeof(frame *) * size);
		memset(pool, 0, sizeof(frame *) * size);
		free = size;
		clock_hand = 0;
		this->size = size;

		start_dec = false;
		scan_nrounds = 0;

		tot_nwrites = 0;
		tot_naccesses = 0;
		tot_nrounds = 0;
		tot_write_nrounds = 0;
	}

	~enhanced_gclock_buffer() {
		numa_free(pool, sizeof(frame *) * size);
	}

	frame *add(frame *entry);

	void print(bool stat_only = false) {
		printf("**************************\n");
		printf("there are %d frames in the buffer\n", size - free);
		printf("there were %ld writes and %ld page accesses, %ld rounds, %ld write rounds\n",
				tot_nwrites, tot_naccesses, tot_nrounds, tot_write_nrounds);
		if (!stat_only)
			for (unsigned i = 0; i < size - free; i++) {
				printf("\toffset: %ld, hits: %d\n",
						pool[i]->get_offset(), pool[i]->getWC());
			}
	}
};

class LP_gclock_buffer: public enhanced_gclock_buffer
{
	pthread_spinlock_t lock;
public:
	LP_gclock_buffer(int size): enhanced_gclock_buffer(size) {
		printf("lock-protected glock buffer\n");
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	frame *add(frame *entry) {
		pthread_spin_lock(&lock);
		frame *ret = enhanced_gclock_buffer::add(entry);
		pthread_spin_unlock(&lock);
		return ret;
	}
};

/*
 * 2. GClock may need to scan many pages in order to evict a page. We need to 
 * roughly sort the pages, so we can avoid scan unnecessary pages.
 *
 * Solution: we split the page set according to the number of cache hits.
 */
class enhanced_gclock_buffer1: public gclock_buffer
{
	int scan_nrounds;	// how many times has the buffer been scanned
	bool start_dec;		// start to decrease the hit count

	long tot_nwrites;
	long tot_naccesses;
	long tot_nrounds;
	long tot_write_nrounds;

	class range_queue: public linked_page_queue {
		unsigned int max_hits;
		unsigned int min_hits;
	public:
		range_queue() {
			max_hits = 0x7fffffff;
			min_hits = 0;
		}

		range_queue(unsigned int min_hits,
				unsigned int max_hits = 0x7fffffff) {
			this->min_hits = min_hits;
			this->max_hits = max_hits;
		}

		void set_min_hits(unsigned int hits) {
			this->min_hits = hits;
		}

		void set_max_hits(unsigned int hits) {
			this->max_hits = hits;
		}

		unsigned get_min_hits() {
			return min_hits;
		}

		unsigned get_max_hits() {
			return max_hits;
		}
	};

	/* 
	 * All pages in a queue have WC values higher than or equal to
	 * the lower bound of the queue. But it's possible that the WC
	 * value of a page is higher than the upper bound of the queue.
	 */
	range_queue *queues;
	int num_queues;

	unsigned free;
	linked_page_queue::iterator clock_hand;
	unsigned size;

	frame *swap(frame *entry);
	void merge_all_queues();
	void merge_queues(int nrounds);
	void add2range(frame *e, unsigned int hits);
public: 
	enhanced_gclock_buffer1(int size,
		const int *ranges, int num_ranges);

	~enhanced_gclock_buffer1() {
		delete [] queues;
	}

	frame *add(frame *entry);

	void print(bool stat_only = false) {
		printf("**************************\n");
		long offset = -1;
		if (clock_hand.curr())
			offset = clock_hand.curr()->get_offset();
		printf("there are %d frames in the buffer1, scan_nrounds: %d, clock_hand: %ld\n",
				size - free, scan_nrounds, offset);
		printf("there were %ld writes and %ld page accesses and %ld rounds, %ld write rounds\n",
				tot_nwrites, tot_naccesses, tot_nrounds, tot_write_nrounds);
		printf("there are %d range queues: \n", num_queues);
		if (!stat_only)
			for (int i = 0; i < num_queues; i++) {
				printf("[%d, %d), size: %d\n", queues[i].get_min_hits(),
						queues[i].get_max_hits(), queues[i].size());
				for (linked_page_queue::iterator it = queues[i].begin(); it.has_next(); ) {
					frame *pg = it.next();
					printf("\toffset: %ld, hits: %d\n", pg->get_offset(), pg->getWC());
				}
			}
	}

	void sanity_check();
};

#endif
