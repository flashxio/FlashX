#ifndef __GLOBAL_CACHED_PRIVATE_H__
#define __GLOBAL_CACHED_PRIVATE_H__

#include "read_private.h"
#include "cache.h"
#include "NUMA_cache.h"

/**
 * This slab allocator allocates IO requests, and all of them are
 * extended requests.
 */
class request_allocator: public obj_allocator<io_request>
{
	class req_initiator: public obj_initiator<io_request>
	{
	public:
		void init(io_request *req) {
			req->init();
		}
	};
public:
	request_allocator(long increase_size,
			long max_size = MAX_SIZE): obj_allocator<io_request>(
				increase_size, max_size, new req_initiator()) {
	}

	virtual int alloc_objs(io_request **reqs, int num) {
		int ret = obj_allocator<io_request>::alloc_objs(reqs, num);
		// Make sure all requests are extended requests.
		for (int i = 0; i < ret; i++) {
			if (!reqs[i]->is_extended_req()) {
				io_request tmp(true);
				*reqs[i] = tmp;
			}
		}
		return ret;
	}

	virtual io_request *alloc_obj() {
		io_request *req = obj_allocator<io_request>::alloc_obj();
		if (!req->is_extended_req()) {
			io_request tmp(true);
			*req = tmp;
		}
		return req;
	}
};

class global_cached_io: public io_interface
{
	int num_waits;
	long cache_size;
	page_cache *global_cache;
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

	request_allocator req_allocator;

	// These are used for implementing sync IO.
	// A thread can only be blocked on one sync IO, and this class can only
	// be used in a single thread.
	pthread_mutex_t sync_mutex;
	pthread_cond_t sync_cond;
	io_request *wait_req;
	int status;

	// This only counts the requests that use the slow path.
	long curr_req_id;

	long num_accesses;

	int cache_hits;
	int num_fast_process;

	/**
	 * It's another version of read() and write(), but it's responsible
	 * for deleting `req'.
	 */
	ssize_t __read(io_request *req, thread_safe_page *p);
	ssize_t __write(io_request *req, thread_safe_page *p,
		std::vector<thread_safe_page *> &dirty_pages);
public:
	global_cached_io(io_interface *, page_cache *cache);

	page_cache *get_global_cache() {
		return global_cache;
	}

	obj_allocator<io_request> *get_req_allocator() {
		return &req_allocator;
	}

	int preload(off_t start, long size);
	io_status access(char *buf, off_t offset, ssize_t size, int access_method);
	/**
	 * A request can access data of arbitrary size and from arbitrary offset.
	 */
	void access(io_request *requests, int num, io_status *status = NULL);
	virtual void flush_requests() {
		underlying->flush_requests();
	}

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
	ssize_t read(io_request &req, thread_safe_page *pages[], int npages, io_request *orig);

	void process_cached_reqs(io_request *cached_reqs[],
			thread_safe_page *cached_pages[], int num_cached_reqs);

	void queue_requests(io_request *reqs[], int num) {
		pending_requests.addByForce(reqs, num);
		// the pending requests are processed in the app thread.
		// If the app thread is waiting for a sync thread to complete,
		// we need to wake it up to process pending requests.
		pthread_mutex_lock(&sync_mutex);
		if (wait_req)
			pthread_cond_signal(&sync_cond);
		pthread_mutex_unlock(&sync_mutex);
	}

	int handle_pending_requests();

	/* When a thread begins, this method will be called. */
	virtual int init() {
		int ret = underlying->init();
		get_global_cache()->init(underlying);
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

	int get_cache_hits() const {
		return cache_hits;
	}

	int get_num_fast_process() const {
		return num_fast_process;
	}

	// These two methods notify of application threads the completion of requests.
	// For global cache, they call invoke() callback directly.
	// For parted global cache, they should send reply messages to remote threads.
	virtual void notify_completion(io_request *req);
	virtual void notify_completion(io_request *reqs[], int num);

	void finalize_partial_request(io_request &partial, io_request *orig);
	void finalize_request(io_request &req);

private:
	// This method can only be called in a single thread.
	void wait4req(io_request *req) {
		pthread_mutex_lock(&sync_mutex);
		wait_req = req;
		while (wait_req) {
			if (!pending_requests.is_empty()) {
				pthread_mutex_unlock(&sync_mutex);
				handle_pending_requests();
				pthread_mutex_lock(&sync_mutex);
			}
			else
				pthread_cond_wait(&sync_cond, &sync_mutex);
		}
		pthread_mutex_unlock(&sync_mutex);
	}

public:
	void wakeup_on_req(io_request *req, int status) {
		pthread_mutex_lock(&sync_mutex);
		assert(req);
		if (wait_req == req) {
			wait_req = NULL;
			this->status = status;
			pthread_cond_signal(&sync_cond);
		}
		else if (wait_req && !pending_requests.is_empty()) {
			// The thread that owns the global cached io is blocked.
			// And there are pending requests, let's wake up the thread.
			pthread_cond_signal(&sync_cond);
		}
		pthread_mutex_unlock(&sync_mutex);
	}

#ifdef STATISTICS
	void print_stat(int nthreads) {
		underlying->print_stat(nthreads);
		static int tot_hits = 0;
		static int seen_threads = 0;
		static int tot_fast_process = 0;
		seen_threads++;
		tot_hits += cache_hits;
		tot_fast_process += num_fast_process;
		if (seen_threads == nthreads) {
			printf("there are %d cache hits\n", tot_hits);
			global_cache->print_stat();
		}
		printf("there are %d waits\n", num_waits);
		printf("There are %d requests processed in the fast path\n", tot_fast_process);
	}
#endif
};

#endif
