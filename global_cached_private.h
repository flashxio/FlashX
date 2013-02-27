#ifndef __GLOBAL_CACHED_PRIVATE_H__
#define __GLOBAL_CACHED_PRIVATE_H__

#include "read_private.h"
#include "cache.h"
#include "associative_cache.h"
#include "hash_index_cache.h"
#include "LRU2Q.h"

enum {
	TREE_CACHE,
	ASSOCIATIVE_CACHE,
	HASH_INDEX_CACHE,
	CUCKOO_CACHE,
	LRU2Q_CACHE,
	GCLOCK_CACHE,
};

class global_cached_io: public io_interface
{
	int num_waits;
	long cache_size;
	static page_cache *global_cache;
	/* the underlying IO. */
	io_interface *underlying;
	callback *cb;

	/**
	 * If a thread wants to issue a request but only allows non-blocking
	 * operations, the request should be added to the queue. All requests
	 * will be issued to the underlying IO in the user's thread when 
	 * the next user IO comes.
	 */
	thread_safe_FIFO_queue<io_request *> pending_requests;

	int cache_hits;

	/**
	 * It's another version of read() and write(), but it's responsible
	 * for deleting `req'.
	 */
	ssize_t __read(io_request *req, thread_safe_page *p);
	ssize_t __write(io_request *req, thread_safe_page *p,
		std::vector<thread_safe_page *> &dirty_pages);
public:
	global_cached_io(io_interface *underlying);
	global_cached_io(io_interface *, long, int, int node_id);

	static page_cache *create_cache(int cache_type, long cache_size,
			int node_id) {
		page_cache *global_cache;
		switch (cache_type) {
#if 0
			// These are just for testing in the single thread.
			case TREE_CACHE:
				global_cache = new tree_cache(cache_size, 0);
				break;
			case CUCKOO_CACHE:
				global_cache = new cuckoo_cache(cache_size);
				break;
			case GCLOCK_CACHE:
				global_cache = new gclock_cache(cache_size);
				break;
#endif
			case LRU2Q_CACHE:
				global_cache = new LRU2Q_cache(cache_size);
				break;
			case ASSOCIATIVE_CACHE:
				global_cache = new associative_cache(cache_size, MAX_CACHE_SIZE, node_id);
				break;
			case HASH_INDEX_CACHE:
				global_cache = new hash_index_cache(cache_size, node_id);
				break;
			default:
				fprintf(stderr, "wrong cache type\n");
				exit(1);
		}
		return global_cache;
	}

	virtual page_cache *get_global_cache() {
		return global_cache;
	}

	int preload(off_t start, long size);
	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method);
	/**
	 * A request can access data of arbitrary size and from arbitrary offset.
	 */
	ssize_t access(io_request *requests, int num);

	/**
	 * One read can access multiple pages while one write can only write
	 * data in a page because if there is a large write (into multiple pages),
	 * the only IO it can cause is to read the first and the last pages if
	 * the offset of the write isn't aligned to a page size. Otherwise,
	 * we can just simply write data to a page without issuing any IO requests.
	 * TODO the only case that we can optimize writes is that a write needs
	 * to touch two pages and the beginning and the end of the write aren't
	 * aligned with a page size.
	 */
	ssize_t read(io_request &req, thread_safe_page *pages[], int npages);
	ssize_t write(io_request &req, thread_safe_page *p,
		std::vector<thread_safe_page *> &dirty_pages);

	void queue_request(io_request *req) {
		pending_requests.addByForce(&req, 1);
	}

	int handle_pending_requests();

	ssize_t get_size() {
		return underlying->get_size();
	}

	/* When a thread begins, this method will be called. */
	virtual int init() {
		int ret = underlying->init();
		global_cache->init(underlying);
		return ret;
	}

	virtual bool set_callback(callback *cb) {
		if (underlying->support_aio())
			this->cb = cb;
		return underlying->support_aio();
	}
	
	virtual callback *get_callback() {
		return cb;
	}

	virtual bool support_aio() {
		return underlying->support_aio();
	}

	virtual void cleanup() {
		underlying->cleanup();
	}

#ifdef STATISTICS
	void print_stat() {
		static int tot_hits = 0;
		static int seen_threads = 0;
		seen_threads++;
		tot_hits += cache_hits;
		if (seen_threads == nthreads)
			printf("there are %d cache hits\n", tot_hits);
		printf("there are %d waits\n", num_waits);
	}
#endif
};

#endif
