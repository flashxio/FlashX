#ifndef __MEMORY_MANAGER__
#define __MEMORY_MANAGER__

#include <vector>

#include "cache.h"

/**
 * manage free pages in the cache.
 * It also allocates pages from the operating system.
 */
class memory_manager
{
	class linked_page {
		linked_page *prev, *next;
		public:
		linked_page() {
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

	long max_size;
	/* the size we have currently allocated. */
	long size;

	linked_page list;
	long num_free_pages;
	pthread_spinlock_t lock;

	std::vector<page_cache *> caches;
public:
	memory_manager(long max_size) {
		this->max_size = max_size;
		size = 0;
		num_free_pages = 0;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	~memory_manager() {
		pthread_spin_destroy(&lock);
	}

	void register_cache(page_cache *cache) {
		caches.push_back(cache);
	}

	void unregister_cache(page_cache *cache) {
	}

	bool get_free_pages(int npages, char **pages, page_cache *cache);

	void free_pages(int npages, char **pages) {
	}

	long average_cache_size() {
		return max_size / caches.size();
	}

	long get_max_size() {
		return max_size;
	}
};

#endif
