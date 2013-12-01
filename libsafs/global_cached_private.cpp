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

#include <vector>
#include <algorithm>

#include "global_cached_private.h"
#include "slab_allocator.h"

// TODO I assume the block size of the RAID array is 16 pages.
const int RAID_BLOCK_SIZE = 16 * PAGE_SIZE;
const int COMPLETE_QUEUE_SIZE = 10240;

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
	req_ext_allocator(int node_id, long increase_size,
			long max_size = INT_MAX): obj_allocator<io_req_extension>(
				std::string("req_ext_allocator-") + itoa(node_id), node_id,
				increase_size, max_size, &initiator) {
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
			new (req) original_io_request();
		}
	} initiator;
public:
	request_allocator(int node_id, long increase_size,
			long max_size = INT_MAX): obj_allocator<original_io_request>(
				std::string("gcached_req_allocator-") + itoa(node_id), node_id,
				increase_size, max_size, &initiator) {
	}
};

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

static thread_safe_page *generic_complete_req(io_request *req,
		thread_safe_page *p, bool lock)
{
	int page_off;
	thread_safe_page *ret = NULL;
	char *req_buf;
	int req_size;

	if (req->within_1page()) {
		page_off = req->get_offset() - ROUND_PAGE(req->get_offset());
		req_buf = req->get_buf();
		req_size = req->get_size();
	}
	else {
		io_request extracted;
		req->extract(p->get_offset(), PAGE_SIZE, extracted);
		page_off = extracted.get_offset() - ROUND_PAGE(extracted.get_offset());
		req_buf = extracted.get_buf();
		req_size = extracted.get_size();
	}

	if (lock)
		p->lock();
	if (req->get_access_method() == WRITE) {
		memcpy((char *) p->get_data() + page_off, req_buf, req_size);
		if (!p->set_dirty(true))
			ret = p;
	}
	else 
		/* I assume the data I read never crosses the page boundary */
		memcpy(req_buf, (char *) p->get_data() + page_off, req_size);
	if (lock)
		p->unlock();
	// TODO this is a bug. If the page is returned, we shouldn't
	// dereference it here.
	p->dec_ref();
	return ret;
}

/**
 * It returns the page that is dirtied by the function for the first time.
 */
static inline thread_safe_page *__complete_req(io_request *orig,
		thread_safe_page *p)
{
	return generic_complete_req(orig, p, true);
}

static inline thread_safe_page *__complete_req_unlocked(io_request *orig,
		thread_safe_page *p)
{
	return generic_complete_req(orig, p, false);
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
		if (this_io->get_callback())
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
		if (is_local_req(reqs[i], this_io))
			local_reqs[num_local++] = reqs[i];
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
	orig->inc_ref();
	if (orig->complete_range(partial.get_offset(), partial.get_size())) {
		orig->set_io(orig->get_orig_io());
		// It's important to notify the IO interface that issues the request.
		// In the case of parted global cache, the IO interface that processes
		// the reqeust isn't the IO interface that issue the request.
		// The request may be handled differently.
		if (orig->is_sync()) {
			assert(orig->get_io() == this);
			((global_cached_io *) orig->get_io())->wakeup_on_req(orig, IO_OK);
			orig->dec_ref();
			orig->wait4unref();
			// Now we can delete it.
			req_allocator->free(orig);
		}
		else {
			num_completed_areqs.inc(1);
			orig->dec_ref();
			complete_queue.push_back(orig);
		}
	}
	else
		orig->dec_ref();
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
	original_io_request *orig = (original_io_request *) request->get_orig();
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
		if (pending_req)
			pending_reqs.push_back(page_req_pair(p, pending_req));
		if (request->get_access_method() == READ) {
			thread_safe_page *dirty = __complete_req_unlocked(orig, p);
			assert(dirty == NULL);
			p->unlock();
		}
		else {
			p->unlock();
			// The page isn't flushed by the page eviction policy.
			// It's flushed because we want to flush data with a large request.
			// The flushed data is contiguous with the old dirty data in
			// the page evicted by the eviction policy.
			if (!orig->contain_offset(p->get_offset()))
				p->dec_ref();
			else
				// The private data contains the page that triggered the writeback.
				assert(request->get_priv() == p);
			assert(p->get_ref() >= 0);
		}
		off += PAGE_SIZE;
	}

	if (request->get_access_method() == READ) {
		/*
		 * For a multi-buf request, the private data actually points to
		 * the very original request.
		 */
		io_request partial;
		orig->extract(request->get_offset(), request->get_num_bufs() * PAGE_SIZE,
				partial);
		this->finalize_partial_request(partial, orig);
	}
	else {
		thread_safe_page *p = (thread_safe_page *) request->get_priv();
		assert(orig->get_next_req_on_page(p) == NULL);
		pending_reqs.push_back(page_req_pair(p, orig));
		// These requests can't be deleted yet.
		// They will be deleted when these write requests are finally served.
	}
	::queue_requests(pending_reqs);
	ext_allocator->free(request->get_extension());

	return -1;
}

void global_cached_io::process_disk_completed_requests(io_request requests[],
		int num)
{
#ifdef STATISTICS
	num_from_underlying.inc(num);
#endif
	std::vector<page_req_pair> pending_reqs;
	for (int i = 0; i < num; i++) {
		io_request *request = &requests[i];

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
		bool data_ready = p->data_ready();
		p->unlock();

		// If the data on the page is ready, it won't become unready.
		// The only place where data is set unready is where the page is evicted.
		// Since we have a reference of the page, it won't be evicted.
		// When data is ready, we can execuate any operations on the page.
		if (data_ready) {
			original_io_request *orig = (original_io_request *) request->get_orig();
			thread_safe_page *dirty = __complete_req(orig, p);
			assert(dirty == NULL);
			io_request partial;
			orig->extract(request->get_offset(),
					request->get_num_bufs() * PAGE_SIZE, partial);
			this->finalize_partial_request(partial, orig);
		}
		else {
			original_io_request *orig = (original_io_request *) request->get_orig();
			assert(orig->get_next_req_on_page(p) == NULL);
			pending_reqs.push_back(page_req_pair(p, orig));
		}
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

int global_cached_io::process_completed_requests()
{
	int num = complete_queue.get_num_entries();
	if (num > 0) {
		stack_array<original_io_request *> reqs(num);
		int ret = complete_queue.fetch(reqs.data(), num);
		::notify_completion(this, (io_request **) reqs.data(), ret);
		for (int i = 0; i < ret; i++)
			req_allocator->free(reqs[i]);
		return ret;
	}
	else
		return 0;
}

global_cached_io::global_cached_io(thread *t, io_interface *underlying,
		page_cache *cache): io_interface(t), pending_requests(
			std::string("pending_req_queue-") + itoa(underlying->get_node_id()),
			underlying->get_node_id(), INIT_GCACHE_PENDING_SIZE),
	complete_queue(std::string("gcached_complete_queue-") + itoa(
				underlying->get_node_id()), underlying->get_node_id(),
			COMPLETE_QUEUE_SIZE, INT_MAX)
{
	assert(t == underlying->get_thread());
	ext_allocator = new req_ext_allocator(underlying->get_node_id(),
			sizeof(io_req_extension) * 1024);
	req_allocator = new request_allocator(underlying->get_node_id(),
			sizeof(original_io_request) * 1024);
	cb = NULL;
	cache_hits = 0;
	num_accesses = 0;
	this->underlying = underlying;
	this->cache_size = cache->size();
	global_cache = cache;
	num_evicted_dirty_pages = 0;
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
				ext->set_orig(orig);
				ext->set_priv(p);
				ext->add_buf((char *) p->get_data(), PAGE_SIZE, 0);

				io_request read_req(ext, pg_loc, READ, this, p->get_node_id());
				p->set_io_pending(true);
				p->unlock();
				io_status status;
				underlying->access(&read_req, 1, &status);
				if (status == IO_FAIL) {
					perror("read");
					exit(1);
				}
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
	}
	return ret;
}

ssize_t global_cached_io::__read(original_io_request *orig, thread_safe_page *p)
{
	ssize_t ret = 0;
	p->lock();
	if (!p->data_ready()) {
		if(!p->is_io_pending()) {
			p->set_io_pending(true);
			assert(!p->is_dirty());

			io_req_extension *ext = ext_allocator->alloc_obj();
			ext->set_orig(orig);
			ext->set_priv(p);
			ext->add_buf((char *) p->get_data(), PAGE_SIZE, 0);

			data_loc_t pg_loc(p->get_file_id(), p->get_offset());
			io_request req(ext, pg_loc, READ, this, get_node_id());
			p->unlock();
			io_status status;
			underlying->access(&req, 1, &status);
			if (status == IO_FAIL) {
				perror("read");
				exit(1);
			}
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

	assert(npages <= MAX_NUM_IOVECS);

	io_req_extension *ext = ext_allocator->alloc_obj();
	ext->set_orig(orig);
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
again:
		p->lock();
		if (!p->data_ready() && !p->is_io_pending()) {
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
			/*
			 * If we have got some partial of the request, we need to submit
			 * the partial request.
			 */
			if (!multibuf_req.is_empty()) {
				p->unlock();
				underlying->access(&multibuf_req, 1);

				io_req_extension *ext = ext_allocator->alloc_obj();
				ext->set_orig(orig);
				io_request tmp(ext, INVALID_DATA_LOC, req.get_access_method(),
						this, get_node_id());
				multibuf_req = tmp;
				goto again;
			}
			else {
				p->add_req(orig);
				p->unlock();
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
				underlying->access(&multibuf_req, 1);

				io_req_extension *ext = ext_allocator->alloc_obj();
				ext->set_orig(orig);
				io_request tmp(ext, INVALID_DATA_LOC, req.get_access_method(),
						this, get_node_id());
				multibuf_req = tmp;
			}
			io_request complete_partial;
			orig->extract(p->get_offset(), PAGE_SIZE, complete_partial);
			ret += complete_partial.get_size();
			__complete_req(&complete_partial, p);
			finalize_partial_request(complete_partial, orig);
		}
	}
	if (!multibuf_req.is_empty()) {
		underlying->access(&multibuf_req, 1);
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
	ext->set_orig(orig);
	ext->set_priv(p);
	io_request req(ext, pg_id, WRITE, this, p->get_node_id());
	assert(p->get_ref() > 0);
	req.add_page(p);
	p->unlock();

	merge_pages2req(req, get_global_cache());
	// The writeback data should have no overlap with the original request
	// that triggered this writeback.
	assert(!req.has_overlap(orig->get_offset(), orig->get_size()));

	io_status status;
	if (orig->is_sync())
		req.set_low_latency(true);
	underlying->access(&req, 1, &status);
	if (status == IO_FAIL) {
		perror("write");
		assert(0);
	}
}

void global_cached_io::process_cached_reqs(io_request *cached_reqs[],
		thread_safe_page *cached_pages[], int num_cached_reqs)
{
	io_request *async_reqs[num_cached_reqs];
	int num_async_reqs = 0;
	num_fast_process += num_cached_reqs;
	for (int i = 0; i < num_cached_reqs; i++) {
		io_request *req = cached_reqs[i];
		thread_safe_page *dirty = __complete_req(req, cached_pages[i]);
		page_cache *cache = get_global_cache();
		if (dirty)
			cache->mark_dirty_pages(&dirty, 1, underlying);
		if (!req->is_sync())
			async_reqs[num_async_reqs++] = req;
	}
	// We don't need to notify completion for sync requests.
	// Actually, we don't even need to do anything for sync requests.
	num_completed_areqs.inc(num_async_reqs);
	::notify_completion(this, async_reqs, num_async_reqs);
}

void global_cached_io::access(io_request *requests, int num, io_status *status)
{
	if (num == 0)
		return;

	ASSERT_EQ(get_thread(), thread::get_curr_thread());
	num_issued_areqs.inc(num);

	stack_array<io_request *, 128> cached_reqs(num);
	stack_array<thread_safe_page *, 128> cached_pages(num);
	int num_cached_reqs = 0;

	bool syncd = false;
	std::vector<thread_safe_page *> dirty_pages;
	for (int i = 0; i < num; i++) {
		off_t offset = requests[i].get_offset();
		int size = requests[i].get_size();
		num_bytes += size;
		// We don't allow the user's requests to be extended requests.
		assert(!requests[i].is_extended_req());
		off_t begin_pg_offset = ROUND_PAGE(offset);
		off_t end_pg_offset = ROUNDUP_PAGE(offset + size);
		thread_safe_page *pages[MAX_NUM_IOVECS];
		// TODO right now it only supports single-buf requests.
		assert(requests[i].get_num_bufs() == 1);
		original_io_request *orig = NULL;

		if (requests[i].is_flush()) {
			syncd = true;
			num_completed_areqs.inc(1);
			continue;
		}
		else if (requests[i].is_sync()) {
			syncd = true;
		}

		int pg_idx = 0;
		int num_pages_hit = 0;
		int num_bytes_completed = 0;
		for (off_t tmp_off = begin_pg_offset; tmp_off < end_pg_offset;
				tmp_off += PAGE_SIZE) {
			thread_safe_page *p;
			
			page_id_t pg_id(requests[i].get_file_id(), tmp_off);
			page_id_t old_id;
			do {
				p = (thread_safe_page *) (get_global_cache()
						->search(pg_id, old_id));
				// If the cache can't evict a page, it's probably because
				// all pages have been referenced. Let's flush all requests
				// from the cached IO and process all completed requests,
				// hopefully we can dereference the pages in the cache.
				if (p == NULL) {
					fprintf(stderr, "can't evict a page\n");
					flush_requests();
					process_all_completed_requests();
				}
			} while (p == NULL);

			num_accesses++;
			if (num_accesses % 100 < params.get_test_hit_rate()) {
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
				assert(requests[i].is_valid());
				if (requests[i].within_1page() && p->data_ready()) {
					cached_reqs[num_cached_reqs] = &requests[i];
					cached_pages[num_cached_reqs] = p;
					num_cached_reqs++;
					break;
				}
			}
			// We delay copying the IO request until here, so we don't
			// need to do it for cached single-page requests..
			if (orig == NULL) {
				orig = req_allocator->alloc_obj();
				orig->init(requests[i]);
				io_interface *orig_io = orig->get_io();
				orig->set_io(this);
				orig->set_orig_io(orig_io);
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
					orig->extract(pages[0]->get_offset(), pg_idx * PAGE_SIZE, req);
					read(req, pages, pg_idx, orig);
					pg_idx = 0;
				}

				/* The page is evicted in this thread */
				if (old_id.get_offset() != ROUND_PAGE(offset) && old_id.get_offset() != -1) {
					/*
					 * Only one thread can come here because only one thread
					 * can evict the dirty page and the thread gets its old
					 * offset, and only this thread can write back the old
					 * dirty page.
					 */
					write_dirty_page(p, old_id, orig);
					continue;
				}
				else {
					// At this moment, the page is being written back to the file
					// by another thread. We should queue the request to tht page,
					// so when the dirty page completes writing back, we can proceed
					// writing.
					p->lock();
					if (p->is_old_dirty()) {
						p->add_req(orig);
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
			if (orig->get_access_method() == WRITE) {
				num_bytes_completed += __write(orig, p, dirty_pages);
			}
			else {
				// Right now, we don't care in which nodes the pages are.
				pages[pg_idx++] = p;
				if (pg_idx == MAX_NUM_IOVECS || (pages[0]->get_offset()
							+ PAGE_SIZE * pg_idx) % RAID_BLOCK_SIZE == 0) {
					io_request req;
					orig->extract(pages[0]->get_offset(), pg_idx * PAGE_SIZE, req);
					num_bytes_completed += read(req, pages, pg_idx, orig);
					pg_idx = 0;
				}
			}
		}
		/*
		 * The only reason that pg_idx > 0 is that there is a large read request.
		 */
		if (pg_idx) {
			io_request req;
			orig->extract(pages[0]->get_offset(), pg_idx * PAGE_SIZE, req);
			read(req, pages, pg_idx, orig);
		}

		// If all pages accessed by the request are in the cache, the request
		// can be completed by the time when the functions returns.
		if (status) {
			if (num_pages_hit == (end_pg_offset - begin_pg_offset) / PAGE_SIZE
					// It's possible that a request is completed in the slow path.
					// The requested pages may become ready in the slow path;
					// or we write the entire page.
					|| num_bytes_completed == requests[i].get_size())
				status[i] = IO_OK;
			else {
				assert(orig);
				status[i] = IO_PENDING;
				status[i].set_priv_data((long) orig);
			}
		}
	}
	process_cached_reqs(cached_reqs.data(), cached_pages.data(), num_cached_reqs);
	get_global_cache()->mark_dirty_pages(dirty_pages.data(),
				dirty_pages.size(), underlying);

	if (syncd)
		underlying->flush_requests();
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

void global_cached_io::process_all_completed_requests()
{
	process_completed_requests();
	if (!pending_requests.is_empty()) {
		handle_pending_requests();
		flush_requests();
	}
}

void global_cached_io::wait4req(original_io_request *req)
{
	while (!req->is_complete()) {
		process_all_completed_requests();
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
	int pending = num_pending_ios();
	num_to_complete = min(pending, num_to_complete);

	process_all_completed_requests();
	int iters = 0;
	while (pending - num_pending_ios() < num_to_complete) {
		iters++;
		get_thread()->wait();
		process_all_completed_requests();
	}
	return pending - num_pending_ios();
}

void global_cached_io::notify_completion(io_request *reqs[], int num)
{
	stack_array<io_request> req_copies(num);
	for (int i = 0; i < num; i++) {
		req_copies[i] = *reqs[i];
		assert(req_copies[i].get_io());
	}

	process_disk_completed_requests(req_copies.data(), num);

	get_thread()->activate();
}
