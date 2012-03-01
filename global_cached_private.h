#ifndef __GLOBAL_CACHED_PRIVATE_H__
#define __GLOBAL_CACHED_PRIVATE_H__

#include "direct_private.h"
#include "cache.h"
#include "tree_cache.h"
#include "associative_cache.h"
#include "cuckoo_cache.h"
#include "LRU2Q.h"

enum {
	TREE_CACHE,
	ASSOCIATIVE_CACHE,
	CUCKOO_CACHE,
	LRU2Q_CACHE,
};

class global_cached_private: public direct_private
{
	int num_waits;
	long cache_size;
	static page_cache *global_cache;
public:
	global_cached_private(const char *names[], int num, long size, int idx,
			int entry_size): direct_private(names, num, size, idx, entry_size) {
		num_waits = 0;
		cache_size = 0;
	}

	static page_cache *create_cache(int cache_type, long cache_size) {
		page_cache *global_cache;
		switch (cache_type) {
			case TREE_CACHE:
				global_cache = new tree_cache(cache_size, 0);
				break;
			case ASSOCIATIVE_CACHE:
				global_cache = new associative_cache(cache_size);
				break;
			case CUCKOO_CACHE:
				global_cache = new cuckoo_cache(cache_size);
				break;
			case LRU2Q_CACHE:
				global_cache = new LRU2Q_cache(cache_size);
				break;
			default:
				fprintf(stderr, "wrong cache type\n");
				exit(1);
		}
		return global_cache;
	}

	global_cached_private(const char *names[], int num, long size, int idx, long cache_size,
			int entry_size, int cache_type): direct_private(names, num, size, idx, entry_size) {
		num_waits = 0;
		this->cache_size = cache_size;
		if (global_cache == NULL) {
			page::allocate_cache(cache_size);
			global_cache = create_cache(cache_type, cache_size);
		}
	}

	virtual page_cache *get_global_cache() {
		return global_cache;
	}

	int preload(off_t start, long size) {
		if (size > cache_size) {
			fprintf(stderr, "we can't preload data larger than the cache size\n");
			exit(1);
		}

		/* open the file. It's a hack, but it works for now. */
		thread_init();

		assert(ROUND_PAGE(start) == start);
		for (long offset = start; offset < start + size; offset += PAGE_SIZE) {
			off_t old_off = -1;
			thread_safe_page *p = (thread_safe_page *) (get_global_cache()->search(ROUND_PAGE(offset), old_off));
			if (!p->data_ready()) {
				ssize_t ret = read_private::access((char *) p->get_data(),
						ROUND_PAGE(offset), PAGE_SIZE, READ);
				if (ret < 0) {
					perror("read");
					return ret;
				}
				p->set_io_pending(false);
				p->set_data_ready(true);
			}
		}
		/* close the file as it will be opened again in the real workload. */
		thread_end();
		return 0;
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method) {
		ssize_t ret;
		off_t old_off = -1;
		thread_safe_page *p = (thread_safe_page *) (get_global_cache()->search(ROUND_PAGE(offset), old_off));

		/*
		 * the page isn't in the cache,
		 * so the cache evict a page and return it to us.
		 */
		if (old_off != ROUND_PAGE(offset) && old_off != -1) {
			/* 
			 * if the new page we get is dirty,
			 * we need to write its data back to the file first.
			 */
			if (p->is_dirty()) {
				unsigned long *l = (unsigned long *) p->get_data();
				unsigned long start = old_off / sizeof(long);
				if (*l != start)
					printf("start: %ld, l: %ld\n", start, *l);
				p->lock();
				read_private::access((char *) p->get_data(),
						old_off, PAGE_SIZE, WRITE);
				p->set_dirty(false);
				p->unlock();
			}
		}

		if (!p->data_ready()) {
			/* if the page isn't io pending,
			 * set it io pending, and return
			 * original result. otherwise,
			 * just return the original value.
			 *
			 * This is an ugly hack, but with this
			 * atomic operation, I can avoid using
			 * locks. */
			if(!p->test_and_set_io_pending()) {
				p->lock();
				/* 
				 * we have to make sure the page is clean
				 * before reading data to the page.
				 */
				while(p->is_dirty()) {}
				ret = read_private::access((char *) p->get_data(),
						ROUND_PAGE(offset), PAGE_SIZE, READ);
				p->unlock();
				if (ret < 0) {
					perror("read");
					exit(1);
				}
				p->set_io_pending(false);
				p->set_data_ready(true);
			}
			else {
				num_waits++;
				p->wait_ready();
			}
		}
		else {
#ifdef STATISTICS
			cache_hits++;
#endif
		}
		int page_off = offset - ROUND_PAGE(offset);
		p->lock();
		if (access_method == WRITE) {
			unsigned long *l = (unsigned long *) ((char *) p->get_data() + page_off);
			unsigned long start = offset / sizeof(long);
			if (*l != start)
				printf("write: start: %ld, l: %ld, offset: %ld\n",
						start, *l, p->get_offset() / sizeof(long));
			memcpy((char *) p->get_data() + page_off, buf, size);
			p->set_dirty(true);
		}
		else 
			/* I assume the data I read never crosses the page boundary */
			memcpy(buf, (char *) p->get_data() + page_off, size);
		p->unlock();
		p->dec_ref();
		ret = size;
		return ret;
	}

#ifdef STATISTICS
	void print_stat() {
		direct_private::print_stat();
		printf("there are %d waits in thread %d\n", num_waits, idx);
	}
#endif
};
page_cache *global_cached_private::global_cache;

#endif
