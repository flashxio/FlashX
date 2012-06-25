#ifndef __GCLOCK_CACHE_H__
#define __GCLOCK_CACHE_H__

#include <deque>

#include "cache.h"

/**
 * this is a test file to measure the cache hit rate of GCLOCK
 * on the random zipfian workload.
 */

class gclock_cache: public page_cache {
	thread_safe_page *pages;		// page arrays to contain all physical pages
	int npages;			// the number of all physical pages
	/* 
	 * the clock header of the CLOCK algorithm,
	 * it is an index of the page array.
	 */
	int clock_hdr;
	std::map<off_t, thread_safe_page *> page_map;
	std::deque<thread_safe_page *> free_pages;

	void evict_pages() {
		while(true) {
			int idx = clock_hdr;
			clock_hdr = (clock_hdr + 1) % npages;

			/* if the page is used by others, skip it. */
			assert(pages[idx].get_ref() >= 0);
			if (pages[idx].get_ref())
				continue;

			/*
			 * decrease the hits of the page.
			 * if the hit is 1, return the page.
			 */
			int hits = pages[idx].get_hits();
			if (hits == 1) {
				pages[idx].reset_hits();
				free_pages.push_back(&pages[idx]);
				return;
			}
			pages[idx].set_hits(hits - 1);
		}
//		/*
//		 * we should be able to evict at least one page.
//		 * if we are here, that means all pages are in use.
//		 * we can't do anything.
//		 */
//		assert(false);
	}

public:
	gclock_cache(long cache_size) {
		printf("gclock cache is used\n");
		npages = cache_size / PAGE_SIZE;
		pages = new thread_safe_page[npages];
		for (int i = 0; i < npages; i++) {
			pages[i] = thread_safe_page(-1, i * PAGE_SIZE);
			free_pages.push_back(&pages[i]);
		}
		printf("there are %ld free pages\n", free_pages.size());
		this->npages = npages;
		this->clock_hdr = 0;
	}

	page *search(off_t offset, off_t &old_off) {
		std::map<off_t, thread_safe_page *>::iterator it = page_map.find(offset);
		if (!(it == page_map.end())) {
			thread_safe_page *pg = it->second;
			pg->inc_ref();
			pg->hit();
			return pg;
		}

		/*
		 * we need to get a free page.
		 * evicted pages will be in free_pages.
		 */
		if (free_pages.empty())
			evict_pages();

		thread_safe_page *pg = free_pages.front();
		assert(pg->get_ref() == 0);
		free_pages.pop_front();
		page_map.insert(std::pair<off_t, thread_safe_page *>(offset, pg));

		pg->reset_hits();
		pg->set_data_ready(false);
		old_off = pg->get_offset();
		pg->set_offset(offset);
		pg->inc_ref();
		pg->hit();
		return pg;
	}
};

#endif
