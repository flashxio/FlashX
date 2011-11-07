#ifndef __CACHE_H__
#define __CACHE_H__

#include <unistd.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>

#include <map>

#define PAGE_SIZE 4096
#define ROUND_PAGE(off) ((off) & (~(PAGE_SIZE - 1)))

class page_cache
{
public:
	virtual ssize_t get_from_cache(char *buf, off_t offset, ssize_t size) = 0;
	virtual struct page *get_empty_page(off_t off) = 0;
};

class page
{
	off_t offset;
	void *data;
	char flags;
public:
	page():offset(-1), data(NULL) { }
	page(off_t off, void *d): offset(off), data(d) { }
//	page(const page &p) { offset = p.offset; data = p.data; }
	void set_offset(off_t off) { offset = off; }
	off_t get_offset() const { return offset; }
	void *get_data() const { return data; }
	bool is_valid() const { return offset != -1; }
	bool data_ready() const { return flags & 0x1; }
	void set_data_ready() { flags = flags | 0x1; }
};

/**
 * This data structure is to implement LRU.
 */
class page_buffer
{
	int size;			// the number of pages that can be buffered
	struct page *buf;	// a circular buffer to keep pages.
	char *pages;		// the pointer to the space for all pages.
	int beg_idx;		// the index of the beginning of the buffer
	int end_idx;		// the index of the end of the buffer

public:
	page_buffer(int size) {
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
	tree_cache(int size): map(), buf(size / PAGE_SIZE) { }
	
	ssize_t get_from_cache(char *buf, const off_t offset,
			const ssize_t size) {
		off_t page_off = ROUND_PAGE(offset);
		std::map<off_t, struct page *>::iterator it = map.find(page_off);
		if (it == map.end())
			return -1;
		struct page *page = (*it).second;
		memcpy(buf, (char *) page->get_data() + (offset - page_off), size);
		return size;
	}

	/**
	 * get an empty page from the cahce for `offset'.
	 * we add the page back to the map to tell other threads
	 * that the request for this page has been issued.
	 * If the cache is full, evict a page and return the evicted page
	 */
	struct page *get_empty_page(off_t offset) {
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
