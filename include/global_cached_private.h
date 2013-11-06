#ifndef __GLOBAL_CACHED_PRIVATE_H__
#define __GLOBAL_CACHED_PRIVATE_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "io_interface.h"
#include "cache.h"
#include "container.h"

class request_allocator;
class req_ext_allocator;

class global_cached_io: public io_interface
{
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

	request_allocator *req_allocator;
	req_ext_allocator *ext_allocator;

	thread_safe_FIFO_queue<io_request> complete_queue;

	// This only counts the requests that use the slow path.
	long curr_req_id;

	long num_accesses;
	size_t num_bytes;		// The number of accessed bytes
	int cache_hits;
	int num_fast_process;
	int num_evicted_dirty_pages;

	// Count the number of async requests.
	atomic_integer num_completed_areqs;
	atomic_integer num_issued_areqs;
#ifdef STATISTICS
	atomic_integer num_from_underlying;
#endif

	/**
	 * It's another version of read() and write(), but it's responsible
	 * for deleting `req'.
	 */
	ssize_t __read(io_request *req, thread_safe_page *p);
	ssize_t __write(io_request *req, thread_safe_page *p,
		std::vector<thread_safe_page *> &dirty_pages);
	int multibuf_completion(io_request *request,
			std::vector<thread_safe_page *> &dirty_pages);
	void process_completed_requests(io_request requests[], int num);
	int process_completed_requests(int num);
	void process_all_completed_requests();

	void wait4req(io_request *req);
public:
	global_cached_io(thread *t, io_interface *, page_cache *cache);

	~global_cached_io();

	page_cache *get_global_cache() {
		return global_cache;
	}

	request_allocator *get_req_allocator() {
		return req_allocator;
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
		get_thread()->activate();
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

	virtual int get_file_id() const {
		return underlying->get_file_id();
	}

	virtual void cleanup() {
		wait4complete(num_pending_ios());
		underlying->cleanup();
	}

	int get_cache_hits() const {
		return cache_hits;
	}

	int get_num_fast_process() const {
		return num_fast_process;
	}

	virtual int wait4complete(int num);
	virtual void notify_completion(io_request *reqs[], int num);
	virtual int num_pending_ios() const {
		assert(num_issued_areqs.get() >= num_completed_areqs.get());
		return num_issued_areqs.get() - num_completed_areqs.get();
	}

	void finalize_partial_request(io_request &partial, io_request *orig);
	void finalize_request(io_request &req);

	void write_dirty_page(thread_safe_page *p, off_t off, io_request *orig);

	void wakeup_on_req(io_request *req, int status) {
		assert(req->is_sync());
		assert(req->is_complete());
		get_thread()->activate();
	}

#ifdef STATISTICS
	void print_stat(int nthreads) {
		underlying->print_stat(nthreads);
		static size_t tot_bytes = 0;
		static size_t tot_accesses = 0;
		static size_t tot_hits = 0;
		static size_t tot_fast_process = 0;
		static int seen_threads = 0;
		seen_threads++;
		tot_bytes += num_bytes;
		tot_accesses += num_accesses;
		tot_hits += cache_hits;
		tot_fast_process += num_fast_process;
		printf("global_cached_io: %d requests are completed from the underlying io\n",
				num_from_underlying.get());
		printf("global_cached_io: There are %d evicted dirty pages\n", num_evicted_dirty_pages);
		if (seen_threads == nthreads) {
			printf("global_cached_io: in total, there are %ld accessed bytes, %ld pages\n",
					tot_bytes, tot_accesses);
			printf("and there are %ld cache hits and %ld processed in the fast path\n",
					tot_hits, tot_fast_process);
			global_cache->print_stat();
		}
	}
#endif

	virtual void print_state() {
		printf("global cached io %d has %d pending reqs and %d reqs from underlying\n",
				get_io_idx(), num_pending_ios(), num_from_underlying.get());
		printf("%d completed pending reqs, %d queued completed reqs from underlying\n",
				pending_requests.get_num_entries(), complete_queue.get_num_entries());
		underlying->print_state();
	}
};

#endif
