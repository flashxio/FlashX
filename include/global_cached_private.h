#ifndef __GLOBAL_CACHED_PRIVATE_H__
#define __GLOBAL_CACHED_PRIVATE_H__

/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <atomic>

#include "io_interface.h"
#include "cache.h"
#include "container.h"
#include "comp_io_scheduler.h"

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
		off_t idx = (pg->get_offset() - first_pg_off) / PAGE_SIZE;
		return status_arr[idx];
	}

	page_status &get_page_status(off_t off) {
		off_t first_pg_off = get_first_page_offset();
		off_t idx = (off - first_pg_off) / PAGE_SIZE;
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
		completed_size = atomic_number<ssize_t>(0);
		orig_io = NULL;
	}

	void init(const io_request &req) {
		// Once an IO request is created, I can't change its type. I have to
		// use this ugly way to change it.
		data_loc_t loc(req.get_file_id(), req.get_offset());
		if (req.get_req_type() == io_request::BASIC_REQ
				|| req.get_req_type() == io_request::USER_COMPUTE) {
			*(io_request *) this = req;
		}
		else
			assert(0);

		completed_size = atomic_number<ssize_t>(0);
		orig_io = NULL;
		status_arr.resize(get_num_covered_pages());
		memset(status_arr.data(), 0,
				sizeof(page_status) * get_num_covered_pages());
	}

	thread_safe_page *complete_req(thread_safe_page *p, bool lock);

	bool complete_page(thread_safe_page *pg) {
		get_page_status(pg).completed = true;
		int size = get_overlap_size(pg);
		ssize_t ret = completed_size.inc(size);
		return ret == get_size();
	}

	bool complete_range(off_t off, size_t size) {
		ssize_t ret = completed_size.inc(size);
		off_t pg_begin = ROUND_PAGE(off);
		off_t pg_end = ROUNDUP_PAGE(off + size);
		while (pg_begin < pg_end) {
			assert(!get_page_status(pg_begin).completed);
			get_page_status(pg_begin).completed = true;
			pg_begin += PAGE_SIZE;
		}
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

	void compute();

	friend class original_req_byte_array;
};

/**
 * This is a page byte array based on the original I/O request.
 */
class original_req_byte_array: public page_byte_array
{
	off_t off;
	size_t valid: 1;
	size_t size: 63;
	embedded_array<thread_safe_page *> pages;

	int get_num_covered_pages() const {
		off_t begin_pg = ROUND_PAGE(off);
		off_t end_pg = ROUNDUP_PAGE(off + size);
		return (end_pg - begin_pg) / PAGE_SIZE;
	}

	original_req_byte_array(original_req_byte_array &arr) {
		this->off = arr.off;
		this->size = arr.size;
		this->pages = arr.pages;
		arr.valid = 0;
		this->valid = 1;
	}
public:
	original_req_byte_array(original_io_request &req) {
		off = req.get_offset();
		size = req.get_size();
		valid = 1;
		int num_pages = req.get_num_covered_pages();
		pages.resize(num_pages);
		for (int i = 0; i < num_pages; i++) {
			pages[i] = req.status_arr[i].pg;
			req.status_arr[i].pg = NULL;
		}
	}

	~original_req_byte_array() {
		if (valid) {
			int num_pages = get_num_covered_pages();
			for (int i = 0; i < num_pages; i++) {
				if (pages[i])
					pages[i]->dec_ref();
			}
		}
	}

	virtual off_t get_offset_in_first_page() const {
		return off % PAGE_SIZE;
	}

	virtual thread_safe_page *get_page(int pg_idx) const {
		return pages[pg_idx];
	}

	virtual size_t get_size() const {
		return size;
	}

	void lock() {
		// TODO
		assert(0);
	}

	void unlock() {
		// TODO
		assert(0);
	}

	page_byte_array *clone() {
		return new original_req_byte_array(*this);
	}
};

class global_cached_io: public io_interface
{
	/**
	 * This class represents a request that is being processed by global_cached_io.
	 */
	class partial_request
	{
		io_request req;
		off_t begin_pg_offset;
		off_t end_pg_offset;
		off_t curr_pg_offset;
		original_io_request *orig;
	public:
		partial_request() {
			begin_pg_offset = 0;
			end_pg_offset = 0;
			curr_pg_offset = 0;
			orig = NULL;
		}

		void init(const io_request &req) {
			this->req = req;
			begin_pg_offset = ROUND_PAGE(req.get_offset());
			end_pg_offset = ROUNDUP_PAGE(req.get_offset() + req.get_size());
			curr_pg_offset = begin_pg_offset;
			orig = NULL;
		}

		bool is_empty() const {
			return curr_pg_offset == end_pg_offset;
		}

		page_id_t get_curr_page_id() const {
			page_id_t pg_id(req.get_file_id(), curr_pg_offset);
			return pg_id;
		}

		void move_next() {
			curr_pg_offset += PAGE_SIZE;
		}

		const io_request &get_request() const {
			return req;
		}

		void init_orig(original_io_request *orig, io_interface *io) {
			this->orig = orig;
			orig->init(req);
			io_interface *orig_io = orig->get_io();
			orig->set_io(io);
			orig->set_orig_io(orig_io);
		}

		original_io_request *get_orig() const {
			return orig;
		}

		size_t get_remaining_size() const {
			off_t curr = curr_pg_offset;
			if (curr < req.get_offset())
				curr = req.get_offset();
			return req.get_offset() + req.get_size() - curr;
		}
	};

	long cache_size;
	page_cache *global_cache;
	/* the underlying IO. */
	io_interface *underlying;
	callback *cb;

	request_allocator *req_allocator;
	req_ext_allocator *ext_allocator;

	// This contains the original requests issued by the application.
	// An original request is placed in this queue when the I/O on a page
	// covered by the request is complete.
	// We need this queue because there is no guarantee that the requests
	// queued on the page are from the same global_cached_io. Once the I/O
	// on the page completes, all requests on the page will be sent to the
	// global_cached_io where the requests were generated.
	thread_safe_FIFO_queue<page_req_pair> pending_requests;

	// It contains the completed asynchronous user requests.
	thread_safe_FIFO_queue<original_io_request *> complete_queue;
	// It contains the completed requests issued to the underlying IO by
	// global_cached_io.
	thread_safe_FIFO_queue<io_request> completed_disk_queue;

	// This contains the requests issued to the underlying IO.
	// There is only one thread that can access the request buffer,
	// it doesn't need to be thread-safe.
	std::vector<io_request> underlying_requests;

	// This contains the requests in the fast process path.
	fifo_queue<std::pair<io_request, thread_safe_page *> > cached_requests;
	// This contains the requests from the application.
	fifo_queue<io_request> user_requests;
	// This is a buffer of requests generated by user tasks.
	fifo_queue<io_request> user_comp_requests;
	// This contains a request from the application. It contains a request
	// in progress.
	partial_request processing_req;
	comp_io_scheduler *comp_io_sched;

	long num_pg_accesses;
	size_t num_bytes;		// The number of accessed bytes
	size_t cache_hits;
	size_t num_fast_process;
	size_t num_evicted_dirty_pages;

	// Count the number of async requests.
	// The number of async requests that have been completed.
	atomic_number<size_t> num_completed_areqs;
	// The number of async requests issued by the application.
	atomic_number<size_t> num_issued_areqs;
	// The number of async requests that have been processed.
	// It is useful when we want to count the number of requests being
	// processed.
	atomic_number<size_t> num_processed_areqs;

	atomic_number<size_t> num_to_underlying;
	atomic_number<size_t> num_from_underlying;
	atomic_number<size_t> num_underlying_pages;

	/**
	 * It's another version of read() and write(), but it's responsible
	 * for deleting `req'.
	 */
	ssize_t __read(original_io_request *req, thread_safe_page *p);
	ssize_t __write(original_io_request *orig, thread_safe_page *p,
		std::vector<thread_safe_page *> &dirty_pages);
	int multibuf_completion(io_request *request);

	void wait4req(original_io_request *req);

	int get_num_underlying_reqs() const {
		return num_to_underlying.get() - num_from_underlying.get();
	}

	void send2underlying(io_request &req) {
		if (params.is_merge_reqs()) {
			underlying_requests.push_back(req);
		}
		else {
			io_status status;
			num_to_underlying.inc(1);
			num_underlying_pages.inc(req.get_num_bufs());
			underlying->access(&req, 1, &status);
			if (status == IO_FAIL) {
				abort();
			}
		}
	}
public:
	global_cached_io(thread *t, io_interface *, page_cache *cache,
			comp_io_scheduler *sched = NULL);

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
	virtual void flush_requests();

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

	// Finish processing cached I/O requests.
	void process_cached_reqs();
	// Process a request from the application.
	void process_user_req(std::vector<thread_safe_page *> &dirty_pages,
			io_status *status);
	// Process the remaining requests issued by the application.
	void process_user_reqs(queue_interface<io_request> &queue);

	// This function performs post-computation steps, after we perform the user
	// computation.
	void complete_user_compute(user_compute *compute);

	void queue_requests(page_req_pair reqs[], int num) {
		int num_added = pending_requests.add(reqs, num);
		assert(num_added == num);
		get_thread()->activate();
	}

	int handle_pending_requests();

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
		// wait4complete may generate more requests because of user compute
		// tasks. We have to make sure all requests are completed.
		while (num_pending_ios() > 0 || !comp_io_sched->is_empty())
			wait4complete(num_pending_ios());
		underlying->cleanup();
		assert(num_processed_areqs.get() == num_completed_areqs.get());
		assert(num_processed_areqs.get() == num_issued_areqs.get());
		assert(get_num_underlying_reqs() == 0);
		assert(num_underlying_pages.get() == 0);
		assert(pending_requests.is_empty());
		assert(complete_queue.is_empty());
		assert(completed_disk_queue.is_empty());
		assert(underlying_requests.empty());
		assert(cached_requests.is_empty());
		assert(user_requests.is_empty());
		assert(processing_req.is_empty());
	}

	size_t get_cache_hits() const {
		return cache_hits;
	}

	size_t get_num_fast_process() const {
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
	 * Process all queued requests.
	 */
	void process_all_requests();

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
		static std::atomic_ulong tot_bytes;
		static std::atomic_ulong tot_accesses;
		static std::atomic_ulong tot_hits;
		static std::atomic_ulong tot_fast_process;
		static std::atomic_int seen_threads;
		tot_bytes += num_bytes;
		tot_accesses += num_pg_accesses;
		tot_hits += cache_hits;
		tot_fast_process += num_fast_process;
		printf("global_cached_io: %ld reqs, %ld bytes, %ld reqs to the underlying io\n",
				num_processed_areqs.get(), num_bytes, num_from_underlying.get());
		printf("global_cached_io: There are %ld evicted dirty pages\n", num_evicted_dirty_pages);
		if (seen_threads.fetch_add(1) == nthreads - 1) {
			printf("global_cached_io: in total, there are %ld accessed bytes, %ld pages\n",
					tot_bytes.load(), tot_accesses.load());
			printf("and there are %ld cache hits and %ld processed in the fast path\n",
					tot_hits.load(), tot_fast_process.load());
			global_cache->print_stat();
		}
	}
#endif

	virtual void print_state() {
#ifdef STATISTICS
		printf("global cached io %d has %d pending reqs and %ld reqs from underlying\n",
				get_io_id(), num_pending_ios(), num_from_underlying.get());
#endif
		printf("%d completed pending reqs, %d queued completed reqs from underlying\n",
				pending_requests.get_num_entries(), complete_queue.get_num_entries());
		underlying->print_state();
	}
};

#endif
