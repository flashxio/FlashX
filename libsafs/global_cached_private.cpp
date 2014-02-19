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

/**
 * This file implements global cache.
 * There are three types of requests:
 * original request: it is a copy of the request passed from access() and 
 *		is allocated in the heap.
 * partial request: it represents part of a request, and it should have 
 *		a point to the original request. This exists in two cases: when
 *		writing an old dirty page; when reading a page with a pending IO.
 * underlying request: it is a request sent to the underlying IO and is
 *		allocated on the stack. It has a point to the original request for
 *		a read request and the partial request for a write request. Only
 *		this type of requests can be multi-buf requests, and it should have
 *		a point to a page if it's a single-buf request.
 *
 * When a request is issued to the cache with access(), we first make a copy
 * of the request as there is no guarantee that the original request from
 * the invoker will be available after access() returns.
 *
 * However, it's often that we need to break a request into smaller ones for
 * different reasons when issuing them to the underlying IO.
 *
 * In some case, a request tries to access a page that another request has
 * been issued to the underlying IO for the page, the request will be added
 * to the page.
 */

#include <vector>
#include <algorithm>

#include "global_cached_private.h"
#include "slab_allocator.h"

const int COMPLETE_QUEUE_SIZE = 10240;
const int REQ_BUF_SIZE = 64;
const int OBJ_ALLOC_INC_SIZE = 1024 * 1024;
const int MAX_COMP_REQS = 1024 * 1024;

/**
 * The initial size of the queue for pending IO requests
 */
const int INIT_GCACHE_PENDING_SIZE = 1000;

void thread_safe_page::add_req(original_io_request *req)
{
	req->set_next_req_on_page(this, reqs);
	reqs = req;
}

original_io_request *thread_safe_page::reset_reqs()
{
	original_io_request *ret = reqs;
	reqs = NULL;
	return ret;
}

original_io_request *thread_safe_page::get_io_req() const
{
	return reqs;
}

class req_ext_allocator: public obj_allocator<io_req_extension>
{
	class ext_initiator: public obj_initiator<io_req_extension>
	{
	public:
		void init(io_req_extension *ext) {
			if (ext->is_valid()) {
				ext->init();
			}
			else {
				new (ext) io_req_extension();
			}
		}
	} initiator;
public:
	req_ext_allocator(int node_id, long max_size = params.get_max_obj_alloc_size(
				)): obj_allocator<io_req_extension>(
				std::string("req_ext_allocator-") + itoa(node_id), node_id,
				OBJ_ALLOC_INC_SIZE, max_size, &initiator) {
	}
	
	virtual io_req_extension *alloc_obj() {
		io_req_extension *ext = obj_allocator<io_req_extension>::alloc_obj();
		assert(ext);
		return ext;
	}
};

/**
 * This slab allocator allocates IO requests, and all of them are
 * extended requests.
 */
class request_allocator: public obj_allocator<original_io_request>
{
	class req_initiator: public obj_initiator<original_io_request>
	{
	public:
		void init(original_io_request *req) {
			if (req->is_initialized())
				req->init();
			else
				new (req) original_io_request();
		}
	} initiator;
public:
	request_allocator(int node_id, long max_size = params.get_max_obj_alloc_size(
				)): obj_allocator<original_io_request>(
				std::string("gcached_req_allocator-") + itoa(node_id), node_id,
				OBJ_ALLOC_INC_SIZE, max_size, &initiator) {
	}

	virtual original_io_request *alloc_obj() {
		original_io_request *req = obj_allocator<original_io_request>::alloc_obj();
		assert(req);
		return req;
	}
};

thread_safe_page *original_io_request::complete_req(thread_safe_page *p,
		bool lock)
{
	if (get_req_type() == io_request::BASIC_REQ) {
		int page_off;
		thread_safe_page *ret = NULL;
		char *req_buf;
		int req_size;

		if (within_1page()) {
			page_off = get_offset() - ROUND_PAGE(get_offset());
			req_buf = get_buf();
			req_size = get_size();
		}
		else {
			io_request extracted;
			extract(p->get_offset(), PAGE_SIZE, extracted);
			page_off = extracted.get_offset() - ROUND_PAGE(extracted.get_offset());
			req_buf = extracted.get_buf();
			req_size = extracted.get_size();
		}

		if (lock)
			p->lock();
		if (get_access_method() == WRITE) {
			memcpy((char *) p->get_data() + page_off, req_buf, req_size);
			if (!p->set_dirty(true))
				ret = p;
		}
		else 
			/* I assume the data I read never crosses the page boundary */
			memcpy(req_buf, (char *) p->get_data() + page_off, req_size);
		if (lock)
			p->unlock();
		return ret;
	}
	else {
		p->inc_ref();
		get_page_status(p).pg = p;
		return NULL;
	}
}

void original_io_request::compute()
{
	assert(this->get_req_type() == io_request::USER_COMPUTE);
	original_req_byte_array byte_arr(this);
	get_compute()->complete_request(byte_arr);
	int num_pages = get_num_covered_pages();
	for (int i = 0; i < num_pages; i++) {
		assert(status_arr[i].pg);
		status_arr[i].pg->dec_ref();
	}
}

/**
 * It returns the page that is dirtied by the function for the first time.
 */
static inline thread_safe_page *__complete_req(original_io_request *orig,
		thread_safe_page *p)
{
	return orig->complete_req(p, true);
}

static inline thread_safe_page *__complete_req_unlocked(
		original_io_request *orig, thread_safe_page *p)
{
	return orig->complete_req(p, false);
}

/**
 * To test whether the request is issued by this_io.
 */
bool is_local_req(io_request *req, io_interface *this_io)
{
	return req->get_io() == this_io;
}

void notify_completion(io_interface *this_io, io_request *req)
{
	if (is_local_req(req, this_io)) {
		if (this_io->get_callback()
				&& req->get_req_type() == io_request::BASIC_REQ)
			this_io->get_callback()->invoke(&req, 1);
	}
	else {
		io_interface *io = req->get_io();
		// We have to set the io instance to its original one.
		req->set_io(io);
		io->notify_completion(&req, 1);
	}
}

static void notify_completion_func(io_interface *io, io_request *reqs[], int num)
{
	io->notify_completion(reqs, num);
}

void notify_completion(io_interface *this_io, io_request *reqs[], int num)
{
	io_request *local_reqs[num];
	int num_local = 0;
	io_request *remote_reqs[num];
	int num_remote = 0;
	for (int i = 0; i < num; i++) {
		if (is_local_req(reqs[i], this_io)) {
			// We only need to invoke the local basic requests.
			// ignore all local user-compute requests.
			if (reqs[i]->get_req_type() == io_request::BASIC_REQ)
				local_reqs[num_local++] = reqs[i];
		}
		else {
			remote_reqs[num_remote++] = reqs[i];
		}
	}

	if (this_io->get_callback() && num_local > 0)
		this_io->get_callback()->invoke(local_reqs, num_local);

	if (num_remote == 0)
		return;

	process_reqs_on_io(remote_reqs, num_remote, notify_completion_func);
}

void global_cached_io::finalize_partial_request(io_request &partial,
		original_io_request *orig)
{
	if (orig->complete_range(partial.get_offset(), partial.get_size())) {
		orig->set_io(orig->get_orig_io());
		// It's important to notify the IO interface that issues the request.
		// In the case of parted global cache, the IO interface that processes
		// the reqeust isn't the IO interface that issue the request.
		// The request may be handled differently.
		if (orig->is_sync()) {
			// The I/O request with user compute can't be a synchronous request.
			assert(orig->get_req_type() == io_request::BASIC_REQ);
			assert(orig->get_io() == this);
			((global_cached_io *) orig->get_io())->wakeup_on_req(orig, IO_OK);
			// Now we can delete it.
			req_allocator->free(orig);
		}
		else {
			complete_queue.push_back(orig);
		}
	}
}

void global_cached_io::finalize_partial_request(thread_safe_page *p,
		original_io_request *orig)
{
	if (orig->within_1page())
		finalize_partial_request(*orig, orig);
	else {
		io_request partial;
		orig->extract(p->get_offset(), PAGE_SIZE, partial);
		finalize_partial_request(partial, orig);
	}
}

void process_page_reqs_on_io(page_req_pair reqs[], int num)
{
	class get_io_func {
	public:
		io_interface *operator()(page_req_pair req) {
			return req.second->get_io();
		}
	} io_func;

	class process_req_func {
	public:
		void operator()(io_interface *io, page_req_pair reqs[], int num) {
			static_cast<global_cached_io *>(io)->queue_requests(reqs, num);
		}
	} proc_func;

	process_reqs_on_io(reqs, num, io_func, proc_func);
}

void queue_requests(std::vector<page_req_pair> &pending_reqs)
{
	if (pending_reqs.empty())
		return;

	int num_orig = pending_reqs.size();
	for (int i = 0; i < num_orig; i++) {
		thread_safe_page *p = pending_reqs[i].first;
		original_io_request *req = pending_reqs[i].second->get_next_req_on_page(p);
		while (req) {
			original_io_request *next = req->get_next_req_on_page(p);
			req->set_next_req_on_page(p, NULL);
			pending_reqs.push_back(page_req_pair(p, req));
			req = next;
		}
		pending_reqs[i].second->set_next_req_on_page(p, NULL);
	}

	process_page_reqs_on_io(pending_reqs.data(), pending_reqs.size());
}

int global_cached_io::multibuf_completion(io_request *request)
{
	/*
	 * Right now the global cache only support normal access().
	 */
	std::vector<page_req_pair> pending_reqs;
	// The pages that are set dirty for the first time.
	off_t off = request->get_offset();
	for (int i = 0; i < request->get_num_bufs(); i++) {
		thread_safe_page *p = request->get_page(i);
		/*
		 * The pages in the buffer of the request are sorted according
		 * to their offsets.
		 */
		assert(p);
		p->lock();
		assert(p->is_io_pending());
		if (request->get_access_method() == READ)
			p->set_data_ready(true);
		else {
			p->set_dirty(false);
			p->set_old_dirty(false);
		}
		p->set_io_pending(false);
		original_io_request *pending_req = p->reset_reqs();
		p->unlock();
		if (pending_req)
			pending_reqs.push_back(page_req_pair(p, pending_req));
		if (request->get_access_method() == WRITE) {
			// The reference count of a dirty page is always 1 + # original
			// requests, so we can decrease the extra reference here.
			p->dec_ref();
			assert(p->get_ref() >= 0);
		}
		off += PAGE_SIZE;
	}
	::queue_requests(pending_reqs);
	ext_allocator->free(request->get_extension());

	return -1;
}

void global_cached_io::process_disk_completed_requests(io_request requests[],
		int num)
{
	num_from_underlying.inc(num);
	std::vector<page_req_pair> pending_reqs;
	for (int i = 0; i < num; i++) {
		io_request *request = &requests[i];
		num_underlying_pages.dec(request->get_num_bufs());

		if (request->get_num_bufs() > 1) {
			multibuf_completion(request);
			continue;
		}

		thread_safe_page *p = (thread_safe_page *) request->get_priv();
		assert(request->get_size() <= PAGE_SIZE);

		p->lock();
		// If we write data to part of a page, we need to first read
		// the entire page to memory first.
		if (request->get_access_method() == READ) {
			p->set_data_ready(true);
		}
		// We just evict a page with dirty data and write the original
		// dirty data in the page to a file.
		else {
			assert(p->is_old_dirty());
			assert(!p->data_ready());
			p->set_old_dirty(false);
		}
		p->set_io_pending(false);
		original_io_request *old = p->reset_reqs();
		p->unlock();

		if (request->get_access_method() == WRITE) {
			// The reference count of a dirty page is always 1 + # original
			// requests, so we can decrease the extra reference here.
			p->dec_ref();
			assert(p->get_ref() >= 0);
		}
		// TODO I can process read requests.

		if (old)
			pending_reqs.push_back(page_req_pair(p, old));
		// These requests can't be deleted yet.
		// They will be deleted when these write requests are finally served.
		ext_allocator->free(request->get_extension());
	}
	if (!pending_reqs.empty()) {
		::queue_requests(pending_reqs);
	}
}

void global_cached_io::process_incomplete_computes(user_comp_req_queue &requests)
{
	// If the request queue is already full, don't do anything.
	if (requests.get_num_entries() >= MAX_COMP_REQS)
		return;

	// We go through all user computations in the queue one by one at most once.
	// For each user computation, we fetch as many requests as possible.
	// We exit the loop if we run out of space to allocate more requests
	// or if we have run through all user computations.
	int size = incomplete_computes.get_num_entries();
	for (int i = 0; i < size; i++) {
		user_compute *compute = incomplete_computes.front();
		int ret = compute->fetch_requests(this, requests);
		// We should limit the number of pending requests generated by
		// user compute.
		if (requests.get_num_entries() >= MAX_COMP_REQS) {
			assert(ret > 0);
			break;
		}

		// We have fetched all requests generated by the computation,
		// we can remove it from the queue now.
		incomplete_computes.pop_front();
		if (compute->has_completed()) {
			compute->set_flag(user_compute::IN_QUEUE, false);
			compute->dec_ref();
			assert(compute->get_ref() == 0);
			compute_allocator *alloc = compute->get_allocator();
			alloc->free(compute);
		}
		else {
			// The computation hasn't been completed, we need to add it
			// back to the queue for further processing.
			incomplete_computes.push_back(compute);
		}
	}
}

/**
 * This function performs post-computation steps, after we perform the user
 * computation. The post-computation steps include:
 *	fetch requests generated by the computation;
 *	deallocate the user computation if it's no longer being used.
 * It returns true if all requests generated by the user computation has been
 * fetched.
 */
bool global_cached_io::complete_user_compute(user_compute *compute,
		user_comp_req_queue &requests)
{
	compute->fetch_requests(this, requests);

	// If we can't fetch all requests generated by the computation,
	// the computation has completed, we need to keep it and it will run
	// when the required data is available.
	if (compute->has_requests() || !compute->has_completed()) {
		// We have to make sure the computation has requested new data
		// successfully, otherwise, it may not be executed again.
		if (!compute->test_flag(user_compute::IN_QUEUE)) {
			compute->inc_ref();
			incomplete_computes.push_back(compute);
			compute->set_flag(user_compute::IN_QUEUE, true);
		}
		return false;
	}

	// If no one else is referencing the user compute, it means it's not
	// in the queue. We can delete it now.
	if (compute->get_ref() == 0) {
		assert(!compute->test_flag(user_compute::IN_QUEUE));
		compute_allocator *alloc = compute->get_allocator();
		alloc->free(compute);
	}
	return true;
}

int global_cached_io::process_completed_requests(user_comp_req_queue &requests)
{
	if (complete_queue.is_empty())
		return 0;

	io_request *reqp_buf[REQ_BUF_SIZE];
	int num_reqs = 0;
	int num_completed = 0;
	while (!complete_queue.is_empty()) {
		original_io_request *reqp = complete_queue.pop_front();
		reqp_buf[num_reqs++] = reqp;
		num_completed++;
		if (reqp->get_req_type() == io_request::USER_COMPUTE) {
			// This is a user-compute request.
			assert(reqp->get_req_type() == io_request::USER_COMPUTE);
			reqp->compute();
			user_compute *compute = reqp->get_compute();
			complete_user_compute(compute, requests);
		}

		if (num_reqs == REQ_BUF_SIZE) {
			::notify_completion(this, reqp_buf, num_reqs);
			for (int i = 0; i < num_reqs; i++) {
				req_allocator->free((original_io_request *) reqp_buf[i]);
			}
			num_reqs = 0;
		}
	}
	num_completed_areqs.inc(num_completed);
	if (num_reqs > 0) {
		::notify_completion(this, reqp_buf, num_reqs);
		for (int i = 0; i < num_reqs; i++) {
			req_allocator->free((original_io_request *) reqp_buf[i]);
		}
	}
	return num_completed;
}

global_cached_io::global_cached_io(thread *t, io_interface *underlying,
		page_cache *cache): io_interface(t), pending_requests(
			std::string("pending_req_queue-") + itoa(underlying->get_node_id()),
			underlying->get_node_id(), INIT_GCACHE_PENDING_SIZE, INT_MAX),
	complete_queue(std::string("gcached_complete_queue-") + itoa(
				underlying->get_node_id()), underlying->get_node_id(),
			COMPLETE_QUEUE_SIZE, INT_MAX),
	completed_disk_queue(std::string("gcached_complete_disk_queue-") + itoa(
				underlying->get_node_id()), underlying->get_node_id(),
			COMPLETE_QUEUE_SIZE, INT_MAX),
	cached_requests(underlying->get_node_id(), COMPLETE_QUEUE_SIZE, true),
	user_requests(underlying->get_node_id(), COMPLETE_QUEUE_SIZE, true),
	incomplete_computes(underlying->get_node_id(), MAX_COMP_REQS, true)
{
	assert(t == underlying->get_thread());
	ext_allocator = new req_ext_allocator(underlying->get_node_id());
	req_allocator = new request_allocator(underlying->get_node_id());
	cb = NULL;

	// Initialize the stat values.
	cache_hits = 0;
	num_pg_accesses = 0;
	num_bytes = 0;
	num_fast_process = 0;
	num_evicted_dirty_pages = 0;

	this->underlying = underlying;
	this->cache_size = cache->size();
	global_cache = cache;
	assert(processing_req.is_empty());
}

global_cached_io::~global_cached_io()
{
	delete underlying;
	delete req_allocator;
}

/**
 * Write the data in the request to the page.
 * @orig: the very original request issued by the user. It may span
 * multiple pages.
 */
ssize_t global_cached_io::__write(original_io_request *orig, thread_safe_page *p,
		std::vector<thread_safe_page *> &dirty_pages)
{
	ssize_t ret = 0;
	p->lock();
	assert(!p->is_old_dirty());
	if (!p->data_ready()) {
		if(!p->is_io_pending()) {
			assert(!p->is_dirty());
			assert(orig->has_overlap(p->get_offset(), PAGE_SIZE));

			// We are going to write to part of a page, therefore,
			// we need to first read the page.
			if (orig->get_offset() > p->get_offset()
					|| orig->get_offset() + orig->get_size()
					< p->get_offset() + PAGE_SIZE) {
				off_t off = orig->get_offset();
				data_loc_t pg_loc(orig->get_file_id(), ROUND_PAGE(off));
				io_req_extension *ext = ext_allocator->alloc_obj();

				io_request read_req(ext, pg_loc, READ, this, p->get_node_id());
				read_req.add_page(p);
				read_req.set_priv(p);
				assert(p->get_io_req() == NULL);
				p->set_io_pending(true);
				p->add_req(orig);
				p->unlock();
				send2underlying(read_req);
			}
			else {
				// This is an optimization. If we can overwrite the entire page,
				// we don't need to read the page first. However, we have to
				// make sure data is written to a page without anyone else
				// having IO operations on it.
				p->set_data_ready(true);
				thread_safe_page *dirty = __complete_req_unlocked(orig, p);
				if (dirty)
					dirty_pages.push_back(dirty);
				p->unlock();
				ret = PAGE_SIZE;
				finalize_partial_request(p, orig);
				// TODO I may need to move page dereference further down.
				// dirty page now doesn't have a reference.
				p->dec_ref();
			}
		}
		else {
			// If there is an IO pending, it means a read request
			// has been issuded. It can't be a write request, otherwise,
			// the data in the page will be ready.
			assert(orig->get_access_method() == WRITE);
			p->add_req(orig);
			p->unlock();
		}
	}
	else {
		// The data in the page is ready. We can write data to the page directly.
		//
		// If data is ready, there shouldn't be an IO pending.
		// In other words, if the thread for writing dirty pages is writing
		// a page, the page will be referenced and therefore, can't be returned
		// from the cache.
		// TODO we should delay the write if the page is being written back.
//		assert(!p->is_io_pending());
		p->unlock();

		thread_safe_page *dirty = __complete_req(orig, p);
		if (dirty)
			dirty_pages.push_back(dirty);
		ret = orig->get_size();
		finalize_partial_request(p, orig);
		// TODO I may need to move page dereference further down.
		// dirty page now doesn't have a reference.
		p->dec_ref();
	}
	return ret;
}

ssize_t global_cached_io::__read(original_io_request *orig, thread_safe_page *p)
{
	ssize_t ret = 0;
	p->lock();
	if (!p->data_ready()) {
		if(!p->is_io_pending()) {
			assert(p->get_io_req() == NULL);
			p->set_io_pending(true);
			assert(!p->is_dirty());

			io_req_extension *ext = ext_allocator->alloc_obj();

			data_loc_t pg_loc(p->get_file_id(), p->get_offset());
			io_request req(ext, pg_loc, READ, this, get_node_id());
			req.set_priv(p);
			req.add_page(p);
			p->add_req(orig);
			p->unlock();
			send2underlying(req);
		}
		else {
			assert(orig->get_access_method() == READ);
			p->add_req(orig);
			p->unlock();
		}
	}
	else {
		// If the data in the page is ready, we don't need to change any state
		// of the page and just read data.
		p->unlock();
		ret = orig->get_size();
		__complete_req(orig, p);
		finalize_partial_request(p, orig);
		p->dec_ref();
	}
	return ret;
}

/**
 * In this method, we are going to issue multi-page read requests.
 * However, we may still break the input request if the data in a page
 * is ready or the page is in the state of IO pending.
 * @req: potentially part of a request.
 */
ssize_t global_cached_io::read(io_request &req, thread_safe_page *pages[],
		int npages, original_io_request *orig)
{
	ssize_t ret = 0;

	assert(npages <= params.get_RAID_block_size());

	io_req_extension *ext = ext_allocator->alloc_obj();
	io_request multibuf_req(ext, INVALID_DATA_LOC, req.get_access_method(), this,
			get_node_id());

	assert(npages > 0);
	int file_id = pages[0]->get_file_id();
	/*
	 * The pages in `pages' should be sorted with their offsets.
	 * We are going to grab multiple locks below. As long as we always
	 * lock pages in the order of page offset, there won't be deadlock.
	 */
	for (int i = 0; i < npages; i++) {
		thread_safe_page *p = pages[i];
		assert(file_id == p->get_file_id());
		p->lock();
		if (!p->data_ready() && !p->is_io_pending()) {
			assert(p->get_io_req() == NULL);
			p->add_req(orig);
			p->set_io_pending(true);
			assert(!p->is_dirty());
			if (multibuf_req.is_empty()) {
				data_loc_t loc(p->get_file_id(), p->get_offset());
				multibuf_req.set_data_loc(loc);
			}
			/* We don't need to worry buffer overflow here. */
			multibuf_req.add_page(p);
			multibuf_req.set_priv(p);
			p->unlock();
		}
		else if (!p->data_ready() && p->is_io_pending()) {
			p->add_req(orig);
			p->unlock();

			/*
			 * If we have got some partial of the request, we need to submit
			 * the partial request.
			 */
			if (!multibuf_req.is_empty()) {
				send2underlying(multibuf_req);

				io_req_extension *ext = ext_allocator->alloc_obj();
				io_request tmp(ext, INVALID_DATA_LOC, req.get_access_method(),
						this, get_node_id());
				multibuf_req = tmp;
			}
		}
		/* 
		 * We have data ready in the page, we are still going to break
		 * the request.
		 */
		else {
			p->unlock();
			/*
			 * If we have collected part of the request, issue the partial
			 * request.
			 */
			if (!multibuf_req.is_empty()) {
				send2underlying(multibuf_req);

				io_req_extension *ext = ext_allocator->alloc_obj();
				io_request tmp(ext, INVALID_DATA_LOC, req.get_access_method(),
						this, get_node_id());
				multibuf_req = tmp;
			}
			io_request complete_partial;
			orig->extract(p->get_offset(), PAGE_SIZE, complete_partial);
			ret += complete_partial.get_size();
			__complete_req(orig, p);
			finalize_partial_request(complete_partial, orig);
			p->dec_ref();
		}
	}
	if (!multibuf_req.is_empty()) {
		send2underlying(multibuf_req);
	}
	else {
		ext_allocator->free(multibuf_req.get_extension());
	}
	return ret;
}

int global_cached_io::handle_pending_requests()
{
	int tot = 0;
	std::vector<thread_safe_page *> dirty_pages;
	while (!pending_requests.is_empty()) {
		page_req_pair reqs[MAX_FETCH_REQS];
		int num = pending_requests.fetch(reqs, MAX_FETCH_REQS);
		for (int i = 0; i < num; i++) {
			// It may be the head of a request list. All requests
			// in the list should point to the same page.
			original_io_request *req = reqs[i].second;
			thread_safe_page *p = reqs[i].first;
			assert(!p->is_old_dirty());
			/**
			 * Right now all pending requests are writes.
			 * All writes are single-buf requests.
			 */
			assert(req->get_next_req_on_page(p) == NULL);
			assert(req->get_io() == this);
			if (req->get_access_method() == WRITE)
				__write(req, p, dirty_pages);
			else
				__read(req, p);
		}
		tot += num;
	}
	// It's not very likely we can get dirty pages here because this is
	// the place where we just finish writing old dirty pages to the disk.
	// The only possible reason is that we happen to overwrite the entire
	// page.
	get_global_cache()->mark_dirty_pages(dirty_pages.data(),
			dirty_pages.size(), underlying);
	return tot;
}

void merge_pages2req(io_request &req, page_cache *cache)
{
	if (!params.is_cache_large_write())
		return;

	thread_safe_page *p;
	off_t off = req.get_offset();
	off_t forward_off = off + PAGE_SIZE;
	off_t block_off = ROUND(off, params.get_RAID_block_size() * PAGE_SIZE);
	off_t block_end_off = block_off + params.get_RAID_block_size() * PAGE_SIZE;
	page_id_t pg_id(req.get_file_id(), forward_off);
	while (forward_off < block_end_off
			&& (p = (thread_safe_page *) cache->search(pg_id))) {
		p->lock();
		if (!p->is_dirty()) {
			p->unlock();
			p->dec_ref();
			break;
		}
		if (!p->is_io_pending()) {
			p->set_io_pending(true);
			req.add_page(p);
		}
		else {
			p->unlock();
			p->dec_ref();
			break;
		}
		p->unlock();
		forward_off += PAGE_SIZE;
		pg_id = page_id_t(req.get_file_id(), forward_off);
	}
	if (off >= PAGE_SIZE) {
		off_t backward_off = off - PAGE_SIZE;
		pg_id = page_id_t(req.get_file_id(), backward_off);
		while (backward_off >= block_off
				&& (p = (thread_safe_page *) cache->search(pg_id))) {
			p->lock();
			if (!p->is_dirty()) {
				p->unlock();
				p->dec_ref();
				break;
			}
			if (!p->is_io_pending()) {
				p->set_io_pending(true);
				req.add_page_front(p);
				req.set_data_loc(pg_id);
			}
			else {
				p->unlock();
				p->dec_ref();
				break;
			}
			p->unlock();
			if (backward_off >= PAGE_SIZE) {
				backward_off -= PAGE_SIZE;
				pg_id = page_id_t(req.get_file_id(), backward_off);
			}
			else
				break;
		}
	}
	assert(req.inside_RAID_block());
}

/**
 * Write the dirty page. If possible, we merge it with pages adjacent to
 * it and write a larger request.
 */
void global_cached_io::write_dirty_page(thread_safe_page *p,
		const page_id_t &pg_id, original_io_request *orig)
{
	p->lock();
	assert(!p->is_io_pending());
	p->set_io_pending(true);

	io_req_extension *ext = ext_allocator->alloc_obj();
	ext->set_priv(p);
	io_request req(ext, pg_id, WRITE, this, p->get_node_id());
	assert(p->get_ref() > 0);
	req.add_page(p);
	p->add_req(orig);
	/*
	 * I need to add another reference.
	 * Normally, the reference count of a page should be the same as the number
	 * of original I/O requests pending on the page. In the case of writing
	 * dirty pages, more dirty pages are merged and write together in the same
	 * request. The merged dirty pages don't have original I/O requests, so
	 * their reference counts are 1 + # original I/O requests.
	 * To add another reference to the page that triggers the write, the
	 * reference count of all dirty pages is 1 + # original I/O requests.
	 * It just simplifies the code of handling the completion of the write
	 * request.
	 */
	p->inc_ref();
	p->unlock();

	merge_pages2req(req, get_global_cache());
	// The writeback data should have no overlap with the original request
	// that triggered this writeback.
	assert(!req.has_overlap(orig->get_offset(), orig->get_size()));

	if (orig->is_sync())
		req.set_low_latency(true);

	/*
	 * We have tried to merge the write request to make it as large as
	 * possible. There is no reason to merge it with other reqeusts again.
	 * Besides, we are flushing the old dirty data, it's unlikely to merge
	 * it with other flushes.
	 */
	io_status status;
	num_to_underlying.inc(1);
	num_underlying_pages.inc(req.get_num_bufs());
	underlying->access(&req, 1, &status);
	if (status == IO_FAIL) {
		abort();
	}
}

class simple_page_byte_array: public page_byte_array
{
	io_request *req;
	thread_safe_page *p;
public:
	simple_page_byte_array(io_request *req, thread_safe_page *p) {
		this->req = req;
		this->p = p;
	}

	virtual void lock() {
		// TODO
		assert(0);
	}

	virtual void unlock() {
		// TODO
		assert(0);
	}

	virtual off_t get_offset_in_first_page() const {
		return req->get_offset() % PAGE_SIZE;
	}

	virtual thread_safe_page *get_page(int idx) const {
		assert(idx == 0);
		return p;
	}

	virtual size_t get_size() const {
		return req->get_size();
	}
};

thread_safe_page *complete_cached_req(io_request *req, thread_safe_page *p)
{
	if (req->get_req_type() == io_request::BASIC_REQ) {
		int page_off;
		thread_safe_page *ret = NULL;
		char *req_buf;
		int req_size;

		page_off = req->get_offset() - ROUND_PAGE(req->get_offset());
		req_buf = req->get_buf();
		req_size = req->get_size();

		p->lock();
		if (req->get_access_method() == WRITE) {
			memcpy((char *) p->get_data() + page_off, req_buf, req_size);
			if (!p->set_dirty(true))
				ret = p;
		}
		else 
			/* I assume the data I read never crosses the page boundary */
			memcpy(req_buf, (char *) p->get_data() + page_off, req_size);
		p->unlock();
		return ret;
	}
	else {
		simple_page_byte_array arr(req, p);
		user_compute *compute = req->get_compute();
		compute->complete_request(arr);
		return NULL;
	}
}

void global_cached_io::process_cached_reqs(user_comp_req_queue &requests)
{
	if (cached_requests.is_empty())
		return;

	int num_async_reqs = 0;
	int num_reqs_in_buf = 0;
	io_request req_buf[REQ_BUF_SIZE];
	io_request *reqp_buf[REQ_BUF_SIZE];
	while (!cached_requests.is_empty()) {
		num_fast_process++;
		std::pair<io_request, thread_safe_page *> pair
			= cached_requests.pop_front();
		thread_safe_page *dirty = complete_cached_req(&pair.first,
				pair.second);
		page_cache *cache = get_global_cache();
		if (dirty)
			cache->mark_dirty_pages(&dirty, 1, underlying);
		pair.second->dec_ref();
		if (!pair.first.is_sync()) {
			req_buf[num_reqs_in_buf] = pair.first;
			reqp_buf[num_reqs_in_buf] = &req_buf[num_reqs_in_buf];
			num_reqs_in_buf++;
			num_async_reqs++;
		}

		if (pair.first.get_req_type() == io_request::USER_COMPUTE) {
			user_compute *compute = pair.first.get_compute();
			complete_user_compute(compute, requests);
		}

		if (num_reqs_in_buf == REQ_BUF_SIZE) {
			::notify_completion(this, reqp_buf, num_reqs_in_buf);
			num_reqs_in_buf = 0;
		}
	}
	// We don't need to notify completion for sync requests.
	// Actually, we don't even need to do anything for sync requests.
	num_completed_areqs.inc(num_async_reqs);
	if (num_reqs_in_buf > 0)
		::notify_completion(this, reqp_buf, num_reqs_in_buf);
}

void global_cached_io::process_user_req(
		std::vector<thread_safe_page *> &dirty_pages, io_status *status)
{
	int pg_idx = 0;
	// In max, we use the number of pages in the RAID block.
	thread_safe_page *pages[params.get_RAID_block_size()];
	int num_pages_hit = 0;
	int num_bytes_completed = 0;
	while (!processing_req.is_empty()) {
		thread_safe_page *p;
		page_id_t pg_id = processing_req.get_curr_page_id();
		page_id_t old_id;
		do {
			p = (thread_safe_page *) (get_global_cache()
					->search(pg_id, old_id));
			// If the cache can't evict a page, it's probably because
			// all pages have been referenced. It's likely that we issued
			// too many requests. Let's stop issuing more requests for now.
			if (p == NULL)
				goto end;
		} while (p == NULL);
		processing_req.move_next();

		num_pg_accesses++;
		if (num_pg_accesses % 100 < params.get_test_hit_rate()) {
			if (!p->data_ready()) {
				p->set_io_pending(false);
				p->set_data_ready(true);
				old_id = page_id_t();
				if (p->is_old_dirty()) {
					p->set_dirty(false);
					p->set_old_dirty(false);
					p->set_io_pending(false);
				}
			}
		}
		/* 
		 * If old_off is -1, it means search() didn't evict a page, i.e.,
		 * it's a cache hit.
		 */
		if (old_id.get_offset() == -1) {
#ifdef STATISTICS
			cache_hits++;
#endif
			num_pages_hit++;
			// Let's optimize for cached single-page requests by stealing
			// them from normal code path of processing them.
			if (processing_req.get_request().within_1page() && p->data_ready()) {
				std::pair<io_request, thread_safe_page *> cached(
						processing_req.get_request(), p);
				cached_requests.push_back(cached);
				break;
			}
		}
		// We delay copying the IO request until here, so we don't
		// need to do it for cached single-page requests..
		if (processing_req.get_orig() == NULL) {
			processing_req.init_orig(req_allocator->alloc_obj(), this);
		}
		/*
		 * Cache may evict a dirty page and return the dirty page
		 * to the user before it is written back to a file.
		 *
		 * We encounter a situation that two threads get the old dirty
		 * evicted page, one thread gets its old offset thanks to
		 * old_off, the other can't, so the other thread has to wait
		 * until the dirty page is written to the file, and we need to
		 * give the page another status to indicate it's an old dirty
		 * page.
		 */

		/* This page has been evicted. */
		if (p->is_old_dirty()) {
			num_evicted_dirty_pages++;
			/*
			 * We got a few contiguous pages for read, so we should split
			 * the request and issue reads for the contiguous pages first.
			 * We always break write requests into pages, so it has to be
			 * read requests.
			 */
			if (pg_idx) {
				io_request req;
				processing_req.get_orig()->extract(pages[0]->get_offset(),
						pg_idx * PAGE_SIZE, req);
				read(req, pages, pg_idx, processing_req.get_orig());
				pg_idx = 0;
			}

			/* The page is evicted in this thread */
			assert(old_id.get_offset() != pg_id.get_offset());
			if (old_id.get_offset() != -1) {
				/*
				 * Only one thread can come here because only one thread
				 * can evict the dirty page and the thread gets its old
				 * offset, and only this thread can write back the old
				 * dirty page.
				 */
				write_dirty_page(p, old_id, processing_req.get_orig());
				continue;
			}
			else {
				// At this moment, the page is being written back to the file
				// by another thread. We should queue the request to tht page,
				// so when the dirty page completes writing back, we can proceed
				// writing.
				p->lock();
				if (p->is_old_dirty()) {
					p->add_req(processing_req.get_orig());
					p->unlock();
					// the request has been added to the page, when the old dirty
					// data is written back to the file, the write request will be
					// reissued to the file.
					continue;
				}
				else
					p->unlock();
			}
		}

		/*
		 * Large access only makes sense for reading. As large writes
		 * essentially overwrite entire pages in the memory, so we may
		 * only need to read the first and the last pages.
		 */
		if (processing_req.get_orig()->get_access_method() == WRITE) {
			num_bytes_completed += __write(processing_req.get_orig(), p,
					dirty_pages);
		}
		else {
			// Right now, we don't care in which nodes the pages are.
			pages[pg_idx++] = p;
			assert(pg_idx <= params.get_RAID_block_size());
			if ((pages[0]->get_offset() + PAGE_SIZE * pg_idx)
					% (params.get_RAID_block_size() * PAGE_SIZE) == 0) {
				io_request req;
				processing_req.get_orig()->extract(pages[0]->get_offset(),
						pg_idx * PAGE_SIZE, req);
				num_bytes_completed += read(req, pages, pg_idx,
						processing_req.get_orig());
				pg_idx = 0;
			}
		}
	}

end:
	/*
	 * The only reason that pg_idx > 0 is that there is a large read request.
	 */
	if (pg_idx) {
		io_request req;
		processing_req.get_orig()->extract(pages[0]->get_offset(),
				pg_idx * PAGE_SIZE, req);
		read(req, pages, pg_idx, processing_req.get_orig());
	}

	// If all pages accessed by the request are in the cache, the request
	// can be completed by the time when the functions returns.
	if (status) {
		if (num_pages_hit == processing_req.get_request().get_num_covered_pages()
				// It's possible that a request is completed in the slow path.
				// The requested pages may become ready in the slow path;
				// or we write the entire page.
				|| num_bytes_completed == processing_req.get_request().get_size())
			*status = IO_OK;
		else {
			assert(processing_req.get_orig());
			*status = IO_PENDING;
			status->set_priv_data((long) processing_req.get_orig());
		}
	}
}

void global_cached_io::process_user_reqs(queue_interface<io_request> &queue)
{
	std::vector<thread_safe_page *> dirty_pages;

	// If we haven't finished processing a request, we need to continue
	// processing it.
	if (!processing_req.is_empty())
		process_user_req(dirty_pages, NULL);

	// Now we can process the remaining requests.
	// If processing_req isn't empty, it's likely because there are too many
	// referenced pages in the cache and we can't evict a page from a page set.
	// So we can stop processing the remaining requests for now.
	while (processing_req.is_empty() && !queue.is_empty()
			// Limit the number of async requests being processed.
			&& num_processed_areqs.get() - num_completed_areqs.get(
				) < get_max_num_pending_ios()
			// TODO the maximal number should be configurable.
			&& num_underlying_pages.get() < 100) {
		io_request req = queue.pop_front();
		num_processed_areqs.inc(1);
		// We don't allow the user's requests to be extended requests.
		assert(!req.is_extended_req());
		assert(!req.is_sync());
		// TODO right now it only supports single-buf requests.
		assert(req.get_num_bufs() == 1);

		if (req.is_flush()) {
			num_completed_areqs.inc(1);
			continue;
		}
		processing_req.init(req);
		num_bytes += req.get_size();
		process_user_req(dirty_pages, NULL);
	}

	get_global_cache()->mark_dirty_pages(dirty_pages.data(),
				dirty_pages.size(), underlying);

	flush_requests();
}

void global_cached_io::access(io_request *requests, int num, io_status *status)
{
	if (num == 0)
		return;

	ASSERT_EQ(get_thread(), thread::get_curr_thread());
	num_issued_areqs.inc(num);

	bool syncd = false;
	std::vector<thread_safe_page *> dirty_pages;

	// The user compute will be referenced by IO requests. I need to
	// increase their references now.
	for (int i = 0; i < num; i++) {
		if (requests[i].get_req_type() == io_request::USER_COMPUTE) {
			user_compute *compute = requests[i].get_compute();
			compute->inc_ref();
		}
	}

	if (!processing_req.is_empty()) {
		process_user_req(dirty_pages, NULL);
		if (!processing_req.is_empty()) {
			user_requests.add(requests, num);
			goto end;
		}
	}

	for (int i = 0; i < num; i++) {
		// We don't allow the user's requests to be extended requests.
		assert(!requests[i].is_extended_req());
		// TODO right now it only supports single-buf requests.
		assert(requests[i].get_num_bufs() == 1);

		if (requests[i].is_flush()) {
			assert(!requests[i].is_sync());
			syncd = true;
			num_completed_areqs.inc(1);
			num_processed_areqs.inc(1);
			continue;
		}
		else if (requests[i].is_sync()) {
			syncd = true;
		}
		assert(processing_req.is_empty());
		num_processed_areqs.inc(1);
		processing_req.init(requests[i]);
		num_bytes += requests[i].get_size();
		io_status *stat_p = NULL;
		if (status)
			stat_p = &status[i];
		process_user_req(dirty_pages, stat_p);
		// We can't process all requests. Let's queue the remaining requests.
		if (!processing_req.is_empty() && i < num - 1) {
			user_requests.add(&requests[i + 1], num - i - 1);
			break;
		}
	}

end:
	get_global_cache()->mark_dirty_pages(dirty_pages.data(),
				dirty_pages.size(), underlying);

	if (syncd)
		flush_requests();
}

io_status global_cached_io::access(char *buf, off_t offset,
		ssize_t size, int access_method)
{
	data_loc_t loc(this->get_file_id(), offset);
	io_request req(buf, loc, size, access_method, this, this->get_node_id(), true);
	io_status status;
	access(&req, 1, &status);
	if (status == IO_PENDING) {
		original_io_request *orig = (original_io_request *) status.get_priv_data();
		assert(orig);
		wait4req(orig);
	}
	// TODO IO may fail, I need to return an error in case it fails.
	status = IO_OK;
	status.set_priv_data(size);
	return status;
}

int global_cached_io::preload(off_t start, long size) {
	if (size > cache_size) {
		fprintf(stderr, "we can't preload data larger than the cache size\n");
		exit(1);
	}

	assert(ROUND_PAGE(start) == start);
	for (long offset = start; offset < start + size; offset += PAGE_SIZE) {
		page_id_t pg_id(get_file_id(), ROUND_PAGE(offset));
		page_id_t old_id;
		thread_safe_page *p = (thread_safe_page *) (get_global_cache()->search(
					pg_id, old_id));
		// This is mainly for testing. I don't need to really read data from disks.
		if (!p->data_ready()) {
			p->set_io_pending(false);
			p->set_data_ready(true);
		}
		p->dec_ref();
	}
	return 0;
}

void global_cached_io::process_all_requests()
{
	// We first process the completed requests from the disk.
	// It will add completed user requests and pending requests to queues
	// for further processing.
	while (!completed_disk_queue.is_empty()) {
		int num = completed_disk_queue.get_num_entries();
		stack_array<io_request> reqs(num);
		int ret = completed_disk_queue.fetch(reqs.data(), num);
		process_disk_completed_requests(reqs.data(), ret);
	}

	// Process the requests that are pending on the pages.
	// It may add completed user requests to queues for further processing. 
	if (!pending_requests.is_empty())
		handle_pending_requests();

	// The number of requests generated by user compute.
	int num_new_reqs;
	do {
		// Let's first try to process incomplete computation.
		// It may generate more user requests to the cache and the underlying IO.
		int orig_num = user_comp_requests.get_num_entries();
		process_incomplete_computes(user_comp_requests);
		num_new_reqs = user_comp_requests.get_num_entries() - orig_num;
		num_issued_areqs.inc(num_new_reqs);

		// Process buffered the requests generated by user compute.
		// We should give the requests generated by user compute higher priority.
		// User compute may keep some extra resources, so we should complete it
		// as quickly as possible.
		// It may add completed user requests to queues for further processing. 
		process_user_reqs(user_comp_requests);
		// If user compute has generated some requests and we have consumed
		// all requests, we can try to generate more requests.
	} while (user_comp_requests.is_empty() && num_new_reqs > 0);

	if (user_comp_requests.is_empty()) {
		// Process buffered user requests.
		// It may add completed user requests to queues for further processing. 
		process_user_reqs(user_requests);
	}

	int orig_num = user_comp_requests.get_num_entries();

	// Process the completed requests served in the cache directly.
	// The user compute may generate more user requests.
	process_cached_reqs(user_comp_requests);

	// Process completed user requests.
	// The user compute may generate more user requests.
	process_completed_requests(user_comp_requests);
	num_new_reqs = user_comp_requests.get_num_entries() - orig_num;
	num_issued_areqs.inc(num_new_reqs);

	// Process buffered the requests generated by user compute.
	// We got more requests generated by user compute. We should try to
	// issue them to the underlying IO if possible.
	// It may add completed user requests to queues for further processing. 
	process_user_reqs(user_comp_requests);

	// Processing the pending requests on the pages might issue
	// more I/O requests.
	flush_requests();
}

void global_cached_io::wait4req(original_io_request *req)
{
	while (!req->is_complete()) {
		process_all_requests();
		if (req->is_complete())
			break;
		get_thread()->wait();
	}
}

/**
 * We wait for at least the specified number of requests to complete.
 */
int global_cached_io::wait4complete(int num_to_complete)
{
	flush_requests();
	size_t prev_completed_areqs = num_completed_areqs.get();
	int pending = num_pending_ios();
	num_to_complete = min(pending, num_to_complete);

	process_all_requests();
	while (num_completed_areqs.get()
			- prev_completed_areqs < (size_t) num_to_complete) {
		// We only wait when there are pending requests in the underlying IO.
		if (num_to_underlying.get() - num_from_underlying.get() > 0) {
			get_thread()->wait();
		}
		process_all_requests();
	}
	return num_completed_areqs.get() - prev_completed_areqs;
}

void global_cached_io::notify_completion(io_request *reqs[], int num)
{
	stack_array<io_request> req_copies(num);
	for (int i = 0; i < num; i++) {
		req_copies[i] = *reqs[i];
		assert(req_copies[i].get_io());
	}
	// By default, global_cached_io processes completed disk requests in
	// the application thread.
	completed_disk_queue.add(req_copies.data(), num);

	get_thread()->activate();
}

static bool cross_RAID_block_bound(off_t off, size_t size)
{
	off_t end = off + size;
	off_t block_begin = ROUND(off, params.get_RAID_block_size() * PAGE_SIZE);
	return end - block_begin > params.get_RAID_block_size() * PAGE_SIZE;
}

static bool merge_req(io_request &merged, const io_request &req)
{
	assert(merged.get_offset() <= req.get_offset());
	assert(merged.get_file_id() == req.get_file_id());

	if (merged.get_offset() + merged.get_size() != req.get_offset()
			|| merged.get_access_method() != req.get_access_method()
			|| merged.is_sync() != req.is_sync()
			|| merged.is_high_prio() != req.is_high_prio()
			|| merged.is_low_latency() != req.is_low_latency())
		return false;

	for (int i = 0; i < req.get_num_bufs(); i++) {
		const io_buf &buf = req.get_io_buf(i);
		thread_safe_page *p = buf.get_page();
		// Assume two requests are overlapped, the overlapped part
		// must point to the same pages.
		// Actually, the real case is that two requests should never
		// overlap. The check in flush_requests() has confirmed it.
		if (merged.contain_offset(p->get_offset())) {
			int idx = (p->get_offset() - merged.get_offset()) / PAGE_SIZE;
			assert(p == merged.get_page(idx));
		}
		else
			merged.add_io_buf(buf);
	}
	// We don't need to reference the page for a multi-buf request.
	merged.set_priv(NULL);
	return true;
}

static inline void access(io_interface *underlying, io_request &req)
{
	io_status status;
	underlying->access(&req, 1, &status);
	if (status == IO_FAIL) {
		abort();
	}
}

void global_cached_io::flush_requests()
{
	if (!params.is_merge_reqs())
		assert(underlying_requests.empty());
	else if (underlying_requests.size() > 0) {
		io_request req = underlying_requests[0];
		int num_sent = 0;
		int num_pages = 0;
		for (unsigned i = 1; i < underlying_requests.size(); i++) {
			io_request *under_req = &underlying_requests[i];

			// If the two requests are not connected
			if (under_req->get_offset()
					!= req.get_offset() + req.get_size()
					// The merged request will cross the boundary of
					// a RAID block.
					|| cross_RAID_block_bound(req.get_offset(),
						req.get_size() + under_req->get_size())
					// We can't merge the two requests.
					|| !merge_req(req, *under_req)) {
				num_sent++;
				num_pages += req.get_num_bufs();
				::access(underlying, req);

				req = *under_req;
			}
			else {
				// If we merge the requests, we need to free the I/O
				// extension of the second request.
				ext_allocator->free(under_req->get_extension());
			}
		}
		underlying_requests.clear();
		num_sent++;
		num_pages += req.get_num_bufs();
		num_to_underlying.inc(num_sent);
		num_underlying_pages.inc(num_pages);
		::access(underlying, req);
	}
	underlying->flush_requests();
}
