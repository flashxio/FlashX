#include "gclock.h"

const int MAX_SCAN_NROUNDS = 1024;

frame *enhanced_gclock_buffer::add(frame *entry)
{
	if (free == 0)
		return swap(entry);

	/*
	 * The real hits is the frame's wcount - scan_nrounds,
	 * so we need to change the wcount accordingly when we insert the frame.
	 */
	entry->incrWC(scan_nrounds);
	pool[size - free] = entry;
	free--;
	return NULL;
}

/**
 * This evicts a frame from the buffer and adds the new frame.
 */
frame *enhanced_gclock_buffer::swap(frame *entry)
{
	int num_writes = 0;
	int num_accesses = 0;
	unsigned int num_pinning = 0;
	for (; ;) {
		if (clock_hand == size) {	// if we have scanned all pages in the buffer
			if (start_dec) {// if we have decreased the hit count of all pages.
				start_dec = false;
				scan_nrounds = 0;
				tot_write_nrounds++;
			}
			scan_nrounds++;
			if (scan_nrounds >= MAX_SCAN_NROUNDS)
				start_dec = true;
			clock_hand = 0;
			tot_nrounds++;
		}

		frame *e = pool[clock_hand++];
		if (e == NULL)
			continue;

		num_accesses++;
		if (start_dec) {
			e->decrWC(scan_nrounds);
			num_writes++;
		}

		int pin_count = e->pinCount();
		if (pin_count == -1) {	// evicted?
			pool[clock_hand - 1] = entry;
			/*
			 * If it's in the decreasing mode, scan_nrounds is
			 * virtually 0.
			 */
			if (!start_dec)
				entry->incrWC(scan_nrounds);
			tot_nwrites += num_writes;
			tot_naccesses += num_accesses;
			return e;
		}

		if (pin_count > 0) {	// pinned?
			if (++num_pinning >= size)
				/* 
				 * If all pages are pinned, we have to wait until
				 * some pages can be unpinned.
				 */
				// TODO wait
				;
			continue;
		}

//		printf("hits: %d, nrounds: %d\n", e->getWC(), scan_nrounds);
		if ((e->getWC() - scan_nrounds <= 0 && !start_dec)
				|| (e->getWC() <= 0 && start_dec)) {
			if (e->tryEvict()) {
				tot_nwrites += num_writes;
				tot_naccesses += num_accesses;
				pool[clock_hand - 1] = entry;
				if (!start_dec)
					entry->incrWC(scan_nrounds);
				return e;
			}
		}
	}	// end for
}

enhanced_gclock_buffer1::enhanced_gclock_buffer1(int size,
		const int *ranges, int num_ranges)
{
	/* initialize the ranges. */
	queues = new range_queue[num_ranges];
	for (int i = 0; i < num_ranges - 1; i++) {
		queues[i].set_min_hits(ranges[i]);
		queues[i].set_max_hits(ranges[i + 1]);
	}
	queues[num_ranges - 1].set_min_hits(ranges[num_ranges - 1]);
	num_queues = num_ranges;

	free = size;
	clock_hand = queues[0].begin();
	this->size = size;

	start_dec = false;
	scan_nrounds = 0;

	tot_naccesses = 0;
	tot_nwrites = 0;
	tot_nrounds = 0;
	tot_write_nrounds = 0;
}

frame *enhanced_gclock_buffer1::add(frame *entry)
{
	if (free == 0)
		return swap(entry);

	/*
	 * The real hits is the frame's wcount - scan_nrounds,
	 * so we need to change the wcount accordingly when we insert the frame.
	 */
	entry->incrWC(scan_nrounds);
	queues[0].push_back(entry);
	free--;
	return NULL;
}

/**
 * Merge all queues to the first queue.
 */
void enhanced_gclock_buffer1::merge_all_queues()
{
	for (int i = 1; i < num_queues; i++) {
		if (queues[i].empty())
			continue;

		queues[0].merge(queues + i);
	}
}

/**
 * Merge the queues whose lower bound is `nrounds' * n to the first queue.
 */
void enhanced_gclock_buffer1::merge_queues(int nrounds)
{
	for (int i = 1; i < num_queues; i++) {
		if (queues[i].empty())
			continue;
		if (queues[i].get_min_hits() % nrounds == 0)
			queues[0].merge(queues + i);
	}
}

/**
 * Add the frame to the queue corresponding to the hits.
 */
void enhanced_gclock_buffer1::add2range(frame *e, unsigned int hits)
{
	for (int i = 1; i < num_queues; i++) {
		if (queues[i].get_min_hits() <= hits
				&& queues[i].get_max_hits() > hits) {
			queues[i].push_back(e);
			break;
		}
	}
}

/**
 * sanity check of the state of the buffer.
 */
void enhanced_gclock_buffer1::sanity_check()
{
	unsigned real_size = 0;
	for (int i = 0; i < num_queues; i++)
		real_size += queues[i].size();
	assert(real_size == size);
	assert(clock_hand.owner() == &queues[0]);
}

/**
 * This evicts a frame from the buffer and adds the new frame.
 * There are multiple page lists to hold pages with different
 * number of cache hits. We scan the pages in the first list,
 * which potentially contains pages with a small number of cache hits.
 * The pages in the lists are loosely organized: pages are reorganized
 * only when lists are scanned.
 *
 * When `scan_nrounds' reaches the upper bound of a range, all pages
 * in the next range needs to be added to the first range. For example,
 * assume we have range [0, 2), [2, 4), [4, 8), [8, ]. When `scan_nranges'
 * reaches 2 * n, we need to add all pages in [2, 4) to the first range. 
 * When `scan_nranges' reaches 4 * n, we need to add all pages in [4, 8)
 * to the first range.
 *
 * When `scan_nranges' reaches the maximal value, we should place all pages 
 * in the first range because we need to change the hit count of all pages
 * anyway.
 *
 * It means we only need to scan the pages in the first range. In the process
 * of scanning the first range, we need to place the pages to the right range.
 */
frame *enhanced_gclock_buffer1::swap(frame *entry)
{
	unsigned int num_pinning = 0;

	/* the number of pages to be accessed before a page can be evicted. */
	int num_accesses = 0;
	/* The number of writes before a page can be evicted. */
	int num_writes = 0;
	for (; ;) {
		/* We have scanned all pages in the first range */
		while (!clock_hand.has_next()) {
			if (start_dec) {// if we have decreased the hit count of all pages.
				start_dec = false;
				scan_nrounds = 0;
				tot_write_nrounds++;
			}
			scan_nrounds++;
			if (scan_nrounds >= MAX_SCAN_NROUNDS) {
				merge_all_queues();
				start_dec = true;
				assert(queues[0].size() == size);
			}
			else
				merge_queues(scan_nrounds);
			clock_hand = queues[0].begin();
			tot_nrounds++;
		}
		sanity_check();

		/* Get the next frame in the list. */
		frame *e = clock_hand.next();
		assert(e);
		num_accesses++;

		// TODO I need to look back later.
		if (start_dec) {
			e->decrWC(scan_nrounds);
			num_writes++;
		}

		int pin_count = e->pinCount();
		if (pin_count == -1) {	// evicted?
			/* Use `entry' to replace `e' in the list. */
			clock_hand.set(entry);
			/*
			 * If it's in the decreasing mode, scan_nrounds is
			 * virtually 0.
			 */
			if (!start_dec)
				entry->incrWC(scan_nrounds);
			tot_naccesses += num_accesses;
			tot_nwrites += num_writes;
			return e;
		}

		if (pin_count > 0) {	// pinned?
			if (++num_pinning >= size)
				/* 
				 * If all pages are pinned, we have to wait until
				 * some pages can be unpinned.
				 */
				// TODO wait
				;
			continue;
		}

		int real_hits = start_dec ? e->getWC() : (e->getWC() - scan_nrounds);
		if (real_hits <= 0) {
			if (e->tryEvict()) {
				clock_hand.set(entry);

				if (!start_dec)
					entry->incrWC(scan_nrounds);
				tot_naccesses += num_accesses;
				tot_nwrites += num_writes;
				return e;
			}
		}
		/*
		 * pages are split into multiple queues according to their WC values
		 * instead of real hit counts.
		 */
		else if (e->getWC() >= queues[0].get_max_hits()) {
			clock_hand.remove();
			add2range(e, e->getWC());
		}
	}	// end for
}

frame *LF_gclock_buffer::add(frame *entry) {
	do {
		int v = free.get();
		if (v == 0)
			return swap(entry);
		if (free.CAS(v, v - 1))
			break;
	} while (true);
	int idx = clock_hand.get();
	while (!pool.CAS(idx % size, NULL, entry))
		idx++;
	clock_hand.inc(1);
	return NULL;
}

frame *LF_gclock_buffer::swap(frame *entry) {
	unsigned int num_pinning = 0;
	unsigned int start = clock_hand.get();
	int num_writes = 0;
	int num_accesses = 0;
	for (unsigned int i = start % size; ; i = (i + 1) % size) {
		frame *e = pool.get(i);
		if (e == NULL)
			continue;

		num_accesses++;
		int pin_count = e->pinCount();
		if (pin_count == -1) {	// evicted?
			if (pool.CAS(i, e, entry)) {
				moveClockHand(i, start);
				return e;
			}
			continue;
		}

		if (pin_count > 0) {	// pinned?
			if (++num_pinning >= size)
				/* 
				 * If all pages are pinned, we have to wait until
				 * some pages can be unpinned.
				 */
				// TODO wait
				;
			continue;
		}

		num_writes++;
		if (e->decrWC() <= 0) {
			if (e->tryEvict() && pool.CAS(i, e, entry)) {
				moveClockHand(i, start);
				tot_nwrites += num_writes;
				tot_naccesses += num_accesses;
				return e;
			}
		}
	}	// end for
}
