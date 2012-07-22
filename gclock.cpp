#include "gclock.h"

const int MAX_SCAN_NROUNDS = 32;

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
	unsigned int num_pinning = 0;
	for (; ;) {
		if (clock_hand == size) {	// if we have scanned all pages in the buffer
			if (start_dec) {// if we have decreased the hit count of all pages.
				start_dec = false;
				scan_nrounds = 0;
			}
			scan_nrounds++;
			if (scan_nrounds >= MAX_SCAN_NROUNDS)
				start_dec = true;
			clock_hand = 0;
		}

		frame *e = pool[clock_hand++];
		if (e == NULL)
			continue;

		if (start_dec)
			e->decrWC(scan_nrounds);

		int pin_count = e->pinCount();
		if (pin_count == -1) {	// evicted?
			pool[clock_hand - 1] = entry;
			/*
			 * If it's in the decreasing mode, scan_nrounds is
			 * virtually 0.
			 */
			if (!start_dec)
				entry->incrWC(scan_nrounds);
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
				pool[clock_hand - 1] = entry;
				if (!start_dec)
					entry->incrWC(scan_nrounds);
				return e;
			}
		}
	}	// end for
}

frame *LF_gclock_buffer::swap(frame *entry) {
	unsigned int num_pinning = 0;
	unsigned int start = clock_hand.get();
	for (unsigned int i = start % size; ; i = (i + 1) % size) {
		frame *e = pool.get(i);
		if (e == NULL)
			continue;

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

		if (e->decrWC() <= 0) {
			if (e->tryEvict() && pool.CAS(i, e, entry)) {
				moveClockHand(i, start);
				return e;
			}
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
