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

typedef std::pair<thread_safe_page *, original_io_request *> page_req_pair;

class original_io_request: public io_request
{
	struct page_status
	{
		thread_safe_page *pg;
		// Point to the next request that queues to the same page.
		original_io_request *next;
		bool completed;

		page_status() {
			pg = NULL;
			next = NULL;
			completed = false;
		}
	};

	/* 
	 * This is to protect the object from being removed
	 * while others are still using it.
	 */
	atomic_number<short> refcnt;
	atomic_number<ssize_t> completed_size;

	embedded_array<page_status> status_arr;

	io_interface *orig_io;

	off_t get_first_page_offset() const {
		off_t mask = PAGE_SIZE - 1;
		mask = ~mask;
		return get_offset() & mask;
	}

	page_status &get_page_status(thread_safe_page *pg) {
		off_t first_pg_off = get_first_page_offset();
		int idx = (pg->get_offset() - first_pg_off) / PAGE_SIZE;
		return status_arr[idx];
	}

	page_status &get_page_status(off_t off) {
		off_t first_pg_off = get_first_page_offset();
		int idx = (off - first_pg_off) / PAGE_SIZE;
		return status_arr[idx];
	}
public:
	original_io_request() {
		orig_io = NULL;
	}

	bool is_initialized() const {
		return status_arr.get_capacity() > 0;
	}

	void init() {
		io_request::init();
		refcnt = atomic_number<short>(0);
		completed_size = atomic_number<ssize_t>(0);
		orig_io = NULL;
	}

	void init(const io_request &req) {
		// Once an IO request is created, I can't change its type. I have to
		// use this ugly way to change it.
		data_loc_t loc(req.get_file_id(), req.get_offset());
		if (req.get_req_type() == io_request::BASIC_REQ) {
			io_request::init(req);
		}
		else if (req.get_req_type() == io_request::USER_COMPUTE) {
			// TODO
			assert(0);
		}
		else
			assert(0);

		refcnt = atomic_number<short>(0);
		completed_size = atomic_number<ssize_t>(0);
		orig_io = NULL;
		status_arr.resize(get_num_covered_pages());
	}

	int inc_ref() {
		return refcnt.inc(1);
	}

	int dec_ref() {
		return refcnt.dec(1);
	}

	void wait4unref() {
		while (refcnt.get() > 0) {}
	}

	thread_safe_page *complete_req(thread_safe_page *p, bool lock);

	bool complete_page(thread_safe_page *pg) {
		get_page_status(pg).completed = true;
		int size = get_overlap_size(pg);
		ssize_t ret = completed_size.inc(size);
		return ret == get_size();
	}

	bool complete_range(off_t off, size_t size) {
		assert(ROUND_PAGE(off) == off);
		assert(size % PAGE_SIZE == 0);
		for (unsigned i = 0; i < size / PAGE_SIZE; i++) {
			get_page_status(off + i * PAGE_SIZE).completed = true;
		}
		ssize_t ret = completed_size.inc(size);
		return ret == get_size();
	}

	bool is_complete() const {
		return completed_size.get() == get_size();
	}

	original_io_request *get_next_req_on_page(thread_safe_page *pg) {
		return get_page_status(pg).next;
	}

	void set_next_req_on_page(thread_safe_page *pg, original_io_request *req) {
		get_page_status(pg).next = req;
	}

	io_interface *get_orig_io() const {
		return orig_io;
	}

	void set_orig_io(io_interface *io) {
		orig_io = io;
	}
};

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
	thread_safe_FIFO_queue<page_req_pair> pending_requests;

	request_allocator *req_allocator;
	req_ext_allocator *ext_allocator;

	// It contains the completed asynchronous user requests.
	thread_safe_FIFO_queue<original_io_request *> complete_queue;

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
	ssize_t __read(original_io_request *req, thread_safe_page *p);
	ssize_t __write(original_io_request *orig, thread_safe_page *p,
		std::vector<thread_safe_page *> &dirty_pages);
	int multibuf_completion(io_request *request);

	void wait4req(original_io_request *req);
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
	ssize_t read(io_request &req, thread_safe_page *pages[], int npages,
			original_io_request *orig);

	void process_cached_reqs(io_request *cached_reqs[],
			thread_safe_page *cached_pages[], int num_cached_reqs);

	void queue_requests(page_req_pair reqs[], int num) {
		pending_requests.addByForce(reqs, num);
		get_thread()->activate();
	}

	int handle_pending_requests();

	/* When a thread begins, this method will be called. */
	virtual int init() {
		int ret = underlying->init();
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

	/**
	 * Process the completed requests issued to the disks.
	 * These requests may be part of the users' requests.
	 */
	void process_disk_completed_requests(io_request requests[], int num);
	/**
	 * Process all completed users' requests.
	 */
	int process_completed_requests();
	/**
	 * Process all completed users' requests as well as pending requests.
	 */
	void process_all_completed_requests();

	bool has_pending_requests() {
		return !pending_requests.is_empty();
	}

	virtual int num_pending_ios() const {
		assert(num_issued_areqs.get() >= num_completed_areqs.get());
		return num_issued_areqs.get() - num_completed_areqs.get();
	}

	void finalize_partial_request(io_request &partial, original_io_request *orig);
	void finalize_partial_request(thread_safe_page *p, original_io_request *orig);

	void write_dirty_page(thread_safe_page *p, const page_id_t &pg_id,
			original_io_request *orig);

	void wakeup_on_req(original_io_request *req, int status) {
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
#ifdef STATISTICS
		printf("global cached io %d has %d pending reqs and %d reqs from underlying\n",
				get_io_idx(), num_pending_ios(), num_from_underlying.get());
#endif
		printf("%d completed pending reqs, %d queued completed reqs from underlying\n",
				pending_requests.get_num_entries(), complete_queue.get_num_entries());
		underlying->print_state();
	}
};

#endif
