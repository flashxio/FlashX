#ifndef __2QLRU__
#define __2QLRU__

#include <deque>
#include "cache.h"

#define RECLAIM_NPAGES 32

class LRU2Q_cache: public page_cache {
	class linked_page: public frame {
	public:
		linked_page() {
		}

		linked_page(off_t offset, char *data): frame(offset, data) {
		}
	};

	/**
	 * This class is to support removing any page in the queue.
	 * It's just for showing the idea of LRU2Q, so performance isn't
	 * very important here.
	 */
	class indexed_page_queue: public linked_page_queue {
		/* 
		 * I need to hide this method so no one will be able to 
		 * get the iterator.
		 */
		iterator begin() {
			return linked_page_queue::begin();
		}

		/* This is to acclerate the deletion of any page. */
		std::map<frame *, linked_obj *const> page_map;
	public:
		indexed_page_queue() {
		}

		linked_obj *const push_back(frame *pg) {
			linked_obj *const obj = linked_page_queue::push_back(pg);
			page_map.insert(std::pair<frame *, linked_obj *const>(pg, obj));
			return obj;
		}

		void pop_front() {
			if (size() <= 0)
				return;
			frame *pg = front();
			remove(pg);
		}

		void remove(frame *pg) {
			std::map<frame *, linked_obj *const>::iterator it
				= page_map.find(pg);
			if (it == page_map.end()) {
				fprintf(stderr, "page %p (offset %ld) doesn't exist in the queue\n",
						pg, pg->get_offset());
				return;
			}
			page_map.erase(pg);
			linked_page_queue::remove(it->second);
		}
	};

	int npages;
	linked_page *pages;
	std::map<off_t, linked_page *> page_map;
	indexed_page_queue free_pages;
	indexed_page_queue active_queue;
	indexed_page_queue inactive_queue;

	void evict_pages() {
		/*
		 * if the inactive list is larger than active list,
		 * we need to start to move pages from the active
		 * list to the inactive list.
		 */
		if (inactive_queue.size() < active_queue.size()) {
			for (int i = 0; i < RECLAIM_NPAGES;) {
				linked_page *pg = (linked_page *) active_queue.front();
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
			linked_page *pg = (linked_page *) inactive_queue.front();
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

	memory_manager *manager;
public:
	LRU2Q_cache(long cache_size) {
		printf("LRU2Q cache is used\n");
		npages = cache_size / PAGE_SIZE;
		pages = new linked_page[npages];
		manager = memory_manager::create(cache_size, -1);
		manager->register_cache(this);
		for (int i = 0; i < npages; i++) {
			char *page = NULL;
			manager->get_free_pages(1, &page, this);
			pages[i] = linked_page(-1, page);
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
			pg->inc_ref();
			return pg;
		}

		/* we need to get a free page. */
		if (free_pages.empty())
			evict_pages();
		linked_page *pg = (linked_page *) free_pages.front();
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

	long size() {
		return ((long) npages) * PAGE_SIZE;
	}
};

#endif
