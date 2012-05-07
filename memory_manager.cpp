#include "memory_manager.h"

const long SHRINK_NPAGES = 1024;
const long INCREASE_SIZE = 1024 * 1024 * 128;

/**
 * get `npages' pages for `request_cache'.
 * In the case of shrinking caches, it makes no sense
 * to shrink the cache that is requesting free pages.
 */
bool memory_manager::get_free_pages(int npages,
		char **pages, page_cache *request_cache) {
	pthread_spin_lock(&lock);
	if (num_free_pages < npages) {
		pthread_spin_unlock(&lock);
		if (size < max_size) {
			char *pages = (char *) numa_alloc_local(INCREASE_SIZE);
			pthread_spin_lock(&lock);
			for (int i = 0; i < INCREASE_SIZE / PAGE_SIZE; i++) {
				linked_page *header = (linked_page *) (pages
						+ PAGE_SIZE * i);
				*header = linked_page();
				list.add_back(header);
			}
			size += INCREASE_SIZE;
			num_free_pages += INCREASE_SIZE / PAGE_SIZE;
		}
		else {
			long size = 0;
			page_cache *cache = NULL;
			/* shrink the cache of the largest size */
			for (unsigned int i = 0; i < caches.size(); i++) {
				if (size < caches[i]->size()) {
					size = caches[i]->size();
					cache = caches[i];
				}
			}
			/* 
			 * if we are going to shrink the cache that requests
			 * free pages, it just fails.
			 */
			if (request_cache == cache) {
				return false;
			}
			int num_shrink = SHRINK_NPAGES;
			if (num_shrink < npages)
				num_shrink = npages;
			char *buf[num_shrink];
			if (!cache->shrink(num_shrink, buf)) {
				return false;
			}
			pthread_spin_lock(&lock);
			for (int i = 0; i < num_shrink; i++) {
				linked_page *header = (linked_page *) buf[i];
				*header = linked_page();
				list.add_back(header);
			}
			num_free_pages += num_shrink;
		}
	}
	/* now it's guaranteed that we have enough free pages. */
	for (int i = 0; i < npages; i++) {
		linked_page *p = list.front();
		p->remove_from_list();
		pages[i] = (char *) p;
	}
	num_free_pages -= npages;
	pthread_spin_unlock(&lock);
	return true;
}

void memory_manager::free_pages(int npages, char **pages) {
	pthread_spin_lock(&lock);
	for (int i = 0; i < npages; i++) {
		linked_page *lp = (linked_page *) pages[i];
		*lp = linked_page();
		list.add_front(lp);
	}
	num_free_pages += npages;
	pthread_spin_unlock(&lock);
}
