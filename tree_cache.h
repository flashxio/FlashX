#ifndef __TREE_CACHE_H__
#define __TREE_CACHE_H__

#include "cache.h"

/**
 * The tree page cache
 */
class tree_cache: public page_cache
{
	page_buffer<page> *buf;
	std::map<off_t, struct page *> map;

public:
	tree_cache(long size, long page_buf) {
		printf("tree cache is used\n");
		buf = new page_buffer<page>(size / PAGE_SIZE, page_buf);
	}

	~tree_cache() {
		delete buf;
	}
	
	page *search(const off_t offset) {
		std::map<off_t, struct page *>::iterator it = map.find(offset);
		if (it == map.end()) {
			page *ret = get_empty_page(offset);
			return ret;
		}
		return (*it).second;
	}

	/**
	 * get an empty page from the cahce for `offset'.
	 * we add the page back to the map to tell other threads
	 * that the request for this page has been issued.
	 * If the cache is full, evict a page and return the evicted page
	 */
	page *get_empty_page(off_t offset) {
		struct page *page = buf->get_empty_page();
		if (page->initialized()) {
			map.erase(page->get_offset());
		}
		page->set_offset(offset);
		page->set_data_ready(false);
		map[offset] = page;
		return page;
	}

#if 0
	int add_page(void *data, const off_t offset) {
		class page page(offset, data);
		class page removed;
		struct page *ret = buf.push_back(page, removed);
		if (removed.is_valid()) {
			/* if the page is valid,
			 * we should remove the page in the map. */
			map.erase(removed.get_offset());
		}
		map[ret->get_offset()] = ret;
		return 0;
	}
#endif
};

#endif
