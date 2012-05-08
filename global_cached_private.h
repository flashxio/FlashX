#ifndef __GLOBAL_CACHED_PRIVATE_H__
#define __GLOBAL_CACHED_PRIVATE_H__

#include "read_private.h"
#include "cache.h"
#include "tree_cache.h"
#include "associative_cache.h"
#include "cuckoo_cache.h"
#include "hash_index_cache.h"
#include "LRU2Q.h"

enum {
	TREE_CACHE,
	ASSOCIATIVE_CACHE,
	HASH_INDEX_CACHE,
	CUCKOO_CACHE,
	LRU2Q_CACHE,
};

class global_cached_private: public thread_private
{
	int num_waits;
	long cache_size;
	static page_cache *global_cache;
	/* the underlying IO. */
	read_private *underlying;
public:
	inline global_cached_private(read_private *underlying,
			int idx, int entry_size): thread_private(idx, entry_size) {
		this->underlying = underlying;
		num_waits = 0;
		cache_size = 0;
	}

	static page_cache *create_cache(int cache_type,
			long cache_size, memory_manager *manager) {
		page_cache *global_cache;
		switch (cache_type) {
			case TREE_CACHE:
				global_cache = new tree_cache(cache_size, 0);
				break;
			case ASSOCIATIVE_CACHE:
				global_cache = new associative_cache(manager);
				break;
			case HASH_INDEX_CACHE:
				global_cache = new hash_index_cache(manager);
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

	global_cached_private(read_private *, int, long, int, int,
			memory_manager *);

	virtual page_cache *get_global_cache() {
		return global_cache;
	}

	int preload(off_t start, long size);
	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method);
	ssize_t access(io_request *requests, int num, int access_method);

	ssize_t get_size() {
		return underlying->get_size();
	}

	virtual int thread_init() {
		thread_private::thread_init();
		return underlying->thread_init();
	}

	virtual bool support_bulk() {
		return underlying->support_bulk();
	}

	virtual void cleanup() {
		underlying->cleanup();
	}

#ifdef STATISTICS
	void print_stat() {
		thread_private::print_stat();
		printf("there are %d waits in thread %d\n", num_waits, idx);
	}
#endif
};

#endif
