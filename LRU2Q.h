#ifndef __2QLRU__
#define __2QLRU__

#include <deque>
#include "cache.h"

#define RECLAIM_NPAGES 32

class linked_page: public page {
	linked_page *prev, *next;
public:
	linked_page() {
		prev = this;
		next = this;
	}

	linked_page(off_t off, long data): page(off, data) {
		prev = this;
		next = this;
	}

	void add_front(linked_page *pg) {
		linked_page *next = this->next;
		pg->next = next;
		pg->prev = this;
		this->next = pg;
		next->prev = pg;
	}

	void add_back(linked_page *pg) {
		linked_page *prev = this->prev;
		pg->next = this;
		pg->prev = prev;
		this->prev = pg;
		prev->next = pg;
	}

	void remove_from_list() {
		linked_page *prev = this->prev;
		linked_page *next = this->next;
		prev->next = next;
		next->prev = prev;
		this->next = this;
		this->prev = this;
	}

	bool is_empty() {
		return this->next == this;
	}

	linked_page *front() {
		return next;
	}

	linked_page *back() {
		return prev;
	}
};

class linked_page_queue {
	linked_page head;
	int _size;
public:
	linked_page_queue() { }

	void push_back(linked_page *pg) {
		head.add_back(pg);
		_size++;
	}

	void pop_front() {
		linked_page *pg = front();
		remove(pg);
	}

	void remove(linked_page *pg) {
		if (pg->is_empty())
			return;
		// TODO do I need to check whether the page is in the queue.
		pg->remove_from_list();
		_size--;
	}

	bool empty() {
		return head.is_empty();
	}

	linked_page *front() {
		return head.front();
	}

	linked_page *back() {
		return head.back();
	}

	int size() {
		return _size;
	}
};

class LRU2Q_cache: public page_cache {
	int npages;
	linked_page *pages;
	std::map<off_t, linked_page *> page_map;
	linked_page_queue free_pages;
	linked_page_queue active_queue;
	linked_page_queue inactive_queue;

	void evict_pages() {
		/*
		 * if the inactive list is larger than active list,
		 * we need to start to move pages from the active
		 * list to the inactive list.
		 */
		if (inactive_queue.size() < active_queue.size()) {
			for (int i = 0; i < RECLAIM_NPAGES;) {
				linked_page *pg = active_queue.front();
				active_queue.pop_front();
				if (pg->referenced()) {
					pg->set_referenced(false);
					active_queue.push_back(pg);
				}
				else {
					inactive_queue.push_back(pg);
					pg->set_active(false);
					i++;
				}
			}
		}

		for (int i = 0; i < RECLAIM_NPAGES;) {
			/* reclaim the oldest pages in the inactive queue. */
			linked_page *pg = inactive_queue.front();
			inactive_queue.pop_front();
			if (pg->referenced()) {
				inactive_queue.push_back(pg);
				pg->set_referenced(false);
			}
			else {
				pg->set_referenced(false);
				free_pages.push_back(pg);
				i++;
				page_map.erase(pg->get_offset());
			}
		}
	}

public:
	LRU2Q_cache(long cache_size) {
		printf("LRU2Q cache is used\n");
		npages = cache_size / PAGE_SIZE;
		pages = new linked_page[npages];
		for (int i = 0; i < npages; i++) {
			pages[i] = linked_page(-1, i * PAGE_SIZE);
			free_pages.push_back(&pages[i]);
		}
	}

	page *search(off_t offset, off_t &old_off) {
		std::map<off_t, linked_page *>::iterator it = page_map.find(offset);
		if (!(it == page_map.end())) {
			linked_page *pg = it->second;
			/* if the page is in the inactive list */
			if (!pg->active()) {
				inactive_queue.remove(pg);
				active_queue.push_back(pg);
				pg->set_active(true);
				pg->set_referenced(false);
			}
			else
				pg->set_referenced(true);
			return pg;
		}

		/* we need to get a free page. */
		if (free_pages.empty())
			evict_pages();
		linked_page *pg = free_pages.front();
		free_pages.pop_front();
		page_map.insert(std::pair<off_t, linked_page *>(offset, pg));

		/* we need to add the page to the inactive queue */
		inactive_queue.push_back(pg);
		pg->set_referenced(true);
		pg->set_active(false);
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
