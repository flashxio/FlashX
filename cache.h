#ifndef __CACHE_H__
#define __CACHE_H__

#include <unistd.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>

#include <map>

#define PAGE_SIZE 4096
#define ROUND_PAGE(off) (((long) off) & (~(PAGE_SIZE - 1)))

class page
{
	off_t offset;
	void *data;
	char flags;
	short refcnt;
public:
	page():offset(-1), data(NULL), flags(0) { }
	page(off_t off, void *d): offset(off), data(d), flags(0) { }
//	page(const page &p) { offset = p.offset; data = p.data; }
	/* when we set the offset,
	 * it means the data in the page becomes invalid. */
	void set_offset(off_t off) {
		offset = off;
		__sync_fetch_and_and(&flags, ~0x1);
	}

	off_t get_offset() const { return offset; }
	void *get_data() const { return data; }
	bool is_valid() const { return offset != -1; }
	bool data_ready() const { return flags & 0x1; }
	bool wait_ready() const {
		while (!data_ready()) {}
	}
	void set_data_ready() {
		__sync_fetch_and_or(&flags, 0x1);
	}
	void inc_ref() {
		__sync_fetch_and_add(&refcnt, 1);
	}
	void dec_ref() {
		__sync_fetch_and_sub(&refcnt, 1);
	}
	short get_ref() {
		return refcnt;
	}
};

class page_cache
{
	pthread_spinlock_t _lock;
public:
	page_cache() {
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

	int lock() {
		return pthread_spin_lock(&_lock);
	}

	int unlock() {
		return pthread_spin_unlock(&_lock);
	}

	~page_cache() {
		pthread_spin_destroy(&_lock);
	}
	virtual page *get_page(off_t offset) = 0;
	virtual page *get_empty_page(off_t off) = 0;
};

/**
 * This data structure is to implement LRU.
 */
class page_buffer
{
	long size;			// the number of pages that can be buffered
	struct page *buf;	// a circular buffer to keep pages.
	char *pages;		// the pointer to the space for all pages.
	int beg_idx;		// the index of the beginning of the buffer
	int end_idx;		// the index of the end of the buffer

public:
	page_buffer(long size) {
		this->size = size;
		buf = new struct page[size];
		pages = (char *) valloc(size * PAGE_SIZE);
		for (int i = 0; i < size; i++) {
			buf[i] = page(-1, &pages[i * PAGE_SIZE]);
		}
		beg_idx = 0;
		end_idx = 0;
	}

	~page_buffer() {
		delete [] buf;
		free(pages);
	}

	bool is_full() const {
		return end_idx - beg_idx == -1
			|| end_idx - beg_idx == size - 1;
	}

	/**
	 * return an empty page.
	 * I expected the page will be filled with data,
	 * so I change the begin and end index of the circular buffer.
	 */
	struct page *get_empty_page() {
		if (is_full()) {
			beg_idx = (beg_idx + 1) % size;
		}
		struct page *ret = &buf[end_idx];
		end_idx = (end_idx + 1) % size;
		return ret;
	}

#if 0
	/* push a page to the buffer.
	 * if the buffer is full, the first page is removed,
	 * and return the user */
	struct page *push_back(const struct page &page, struct page &beg) {
		struct page *ret;
		if (is_full()) {
			beg = buf[beg_idx];
			beg_idx = (beg_idx + 1) % size;
		}
		buf[end_idx] = page;
		ret = &buf[end_idx];
		end_idx = (end_idx + 1) % size;
		return ret;
	}
#endif
};

/**
 * The tree page cache
 */
class tree_cache: public page_cache
{
	page_buffer buf;
	std::map<off_t, struct page *> map;

public:
	tree_cache(long size): map(), buf(size / PAGE_SIZE) { }
	
	page *get_page(const off_t offset) {
		std::map<off_t, struct page *>::iterator it = map.find(offset);
		if (it == map.end())
			return NULL;
		return (*it).second;
	}

	/**
	 * get an empty page from the cahce for `offset'.
	 * we add the page back to the map to tell other threads
	 * that the request for this page has been issued.
	 * If the cache is full, evict a page and return the evicted page
	 */
	page *get_empty_page(off_t offset) {
		struct page *page = buf.get_empty_page();
		if (page->is_valid()) {
			map.erase(page->get_offset());
		}
		page->set_offset(offset);
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
