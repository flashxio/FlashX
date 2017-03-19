/*
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

/*
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
#include <system_error>

#include "global_cached_private.h"
#include "slab_allocator.h"

namespace safs
{

static const int COMPLETE_QUEUE_SIZE = 10240;
const int REQ_BUF_SIZE = 64;
const int OBJ_ALLOC_INC_SIZE = 1024 * 1024;

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
			throw std::invalid_argument("wrong request type");

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

	void compute(byte_array_allocator &alloc);

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

	void assign(original_req_byte_array &arr) {
		this->off = arr.off;
		this->size = arr.size;
		this->pages = arr.pages;
		arr.valid = 0;
		this->valid = 1;
	}

	original_req_byte_array(original_req_byte_array &arr) {
		assign(arr);
	}

	original_req_byte_array &operator=(original_req_byte_array &arr) {
		assign(arr);
		return *this;
	}
public:
	original_req_byte_array(byte_array_allocator &alloc): page_byte_array(alloc) {
		off = 0;
		valid = 0;
		size = 0;
	}

	original_req_byte_array(original_io_request &req,
			byte_array_allocator &alloc): page_byte_array(alloc) {
		init(req);
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

	void init(original_io_request &req) {
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

	virtual off_t get_offset() const {
		return off;
	}

	virtual off_t get_offset_in_first_page() const {
		return off % PAGE_SIZE;
	}

	virtual const char *get_page(int pg_idx) const {
		return (const char *) pages[pg_idx]->get_data();
	}

	virtual size_t get_size() const {
		return size;
	}

	void lock() {
		// TODO
		throw unsupported_exception("lock");
	}

	void unlock() {
		// TODO
		throw unsupported_exception("unlock");
	}

	page_byte_array *clone() {
		original_req_byte_array *arr
			= (original_req_byte_array *) get_allocator().alloc();
		*arr = *this;
		return arr;
	}
};

class simple_page_byte_array: public page_byte_array
{
	off_t off;
	size_t size;
	thread_safe_page *p;

	void assign(simple_page_byte_array &arr) {
		this->off = arr.off;
		this->size = arr.size;
		this->p = arr.p;
		arr.p = NULL;
	}

	simple_page_byte_array(simple_page_byte_array &arr) {
		assign(arr);
	}

	simple_page_byte_array &operator=(simple_page_byte_array &arr) {
		assign(arr);
		return *this;
	}
public:
	simple_page_byte_array(byte_array_allocator &alloc): page_byte_array(alloc) {
		off = 0;
		size = 0;
		p = NULL;
	}

	simple_page_byte_array(const io_request &req, thread_safe_page *p,
			byte_array_allocator &alloc): page_byte_array(alloc) {
		init(req, p);
	}

	~simple_page_byte_array() {
		if (p)
			p->dec_ref();
	}

	void init(const io_request &req, thread_safe_page *p) {
		off = req.get_offset();
		size = req.get_size();
		this->p = p;
	}

	virtual void lock() {
		// TODO
		throw unsupported_exception("lock");
	}

	virtual void unlock() {
		// TODO
		throw unsupported_exception("unlock");
	}

	virtual off_t get_offset() const {
		return off;
	}

	virtual off_t get_offset_in_first_page() const {
		return off % PAGE_SIZE;
	}

	virtual const char *get_page(int idx) const {
		assert(idx == 0);
		return (const char *) p->get_data();
	}

	virtual size_t get_size() const {
		return size;
	}

	page_byte_array *clone() {
		simple_page_byte_array *arr
			= (simple_page_byte_array *) get_allocator().alloc();
		*arr = *this;
		return arr;
	}
};

template<class array_type>
class byte_array_allocator_impl: public byte_array_allocator
{
	class array_initiator: public obj_initiator<array_type>
	{
		byte_array_allocator_impl<array_type> *alloc;
	public:
		array_initiator(byte_array_allocator_impl<array_type> *alloc) {
			this->alloc = alloc;
		}

		virtual void init(array_type *obj) {
			new (obj) array_type(*alloc);
		}
	};

	class array_destructor: public obj_destructor<array_type>
	{
	public:
		void destroy(array_type *obj) {
			obj->~array_type();
		}
	};

	obj_allocator<array_type> allocator;
public:
	byte_array_allocator_impl(thread *t): allocator(
			"byte-array-allocator", t->get_node_id(), false, 1024 * 1024,
			params.get_max_obj_alloc_size(),
			typename obj_initiator<array_type>::ptr(new array_initiator(this)),
			typename obj_destructor<array_type>::ptr(new array_destructor())) {
	}

	virtual page_byte_array *alloc() {
		return allocator.alloc_obj();
	}

	virtual void free(page_byte_array *arr) {
		allocator.free((array_type *) arr);
	}
};

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
			new (ext) io_req_extension();
		}
	};

	class ext_destructor: public obj_destructor<io_req_extension>
	{
	public:
		void destroy(io_req_extension *ext) {
			ext->~io_req_extension();
		}
	};
public:
	req_ext_allocator(int node_id, long max_size = params.get_max_obj_alloc_size(
				)): obj_allocator<io_req_extension>(
				std::string("req_ext_allocator-") + itoa(node_id), node_id,
				true, OBJ_ALLOC_INC_SIZE, max_size,
				obj_initiator<io_req_extension>::ptr(new ext_initiator()),
				obj_destructor<io_req_extension>::ptr(new ext_destructor())) {
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
			new (req) original_io_request();
		}
	};

	class req_destructor: public obj_destructor<original_io_request>
	{
	public:
		void destroy(original_io_request *req) {
			req->~original_io_request();
		}
	};
public:
	request_allocator(int node_id, long max_size = params.get_max_obj_alloc_size(
				)): obj_allocator<original_io_request>(
				std::string("gcached_req_allocator-") + itoa(node_id), node_id,
				true, OBJ_ALLOC_INC_SIZE, max_size,
				obj_initiator<original_io_request>::ptr(new req_initiator()),
				obj_destructor<original_io_request>::ptr(new req_destructor())) {
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

void original_io_request::compute(byte_array_allocator &alloc)
{
	assert(this->get_req_type() == io_request::USER_COMPUTE);
	original_req_byte_array byte_arr(*this, alloc);
	get_compute()->run(byte_arr);
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
		if (this_io->have_callback()
				&& req->get_req_type() == io_request::BASIC_REQ)
			this_io->get_callback().invoke(&req, 1);
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

	if (this_io->have_callback() && num_local > 0)
		this_io->get_callback().invoke(local_reqs, num_local);

	if (num_remote == 0)
		return;

	process_reqs_on_io(remote_reqs, num_remote, notify_completion_func);
}

void global_cached_io::partial_request::init_orig(original_io_request *orig,
		io_interface *io)
{
	this->orig = orig;
	orig->init(req);
	io_interface *orig_io = orig->get_io();
	orig->set_io(io);
	orig->set_orig_io(orig_io);
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
			// The sync I/O request should be deleted in wait4req.
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
	safs::queue_requests(pending_reqs);
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
		safs::queue_requests(pending_reqs);
	}
}

int global_cached_io::process_completed_requests()
{
	if (complete_queue.is_empty())
		return 0;

	io_request *reqp_buf[REQ_BUF_SIZE];
	int num_reqs = 0;
	int num_completed = 0;
	while (!complete_queue.is_empty()) {
		original_io_request *reqp = complete_queue.pop_front();
		assert(!reqp->is_sync());
		num_completed++;
		if (reqp->get_req_type() == io_request::USER_COMPUTE) {
			// This is a user-compute request.
			assert(reqp->get_req_type() == io_request::USER_COMPUTE);
			reqp->compute(*orig_array_allocator);
			user_compute *compute = reqp->get_compute();
			comp_io_sched->post_comp_process(compute);
			req_allocator->free(reqp);
		}
		else
			reqp_buf[num_reqs++] = reqp;

		if (num_reqs == REQ_BUF_SIZE) {
			safs::notify_completion(this, reqp_buf, num_reqs);
			for (int i = 0; i < num_reqs; i++) {
				req_allocator->free((original_io_request *) reqp_buf[i]);
			}
			num_reqs = 0;
		}
	}
	num_completed_areqs.inc(num_completed);
	if (num_reqs > 0) {
		safs::notify_completion(this, reqp_buf, num_reqs);
		for (int i = 0; i < num_reqs; i++) {
			req_allocator->free((original_io_request *) reqp_buf[i]);
		}
	}
	return num_completed;
}

global_cached_io::global_cached_io(thread *t, io_interface::ptr underlying,
		page_cache::ptr cache, comp_io_scheduler::ptr sched): io_interface(t,
			underlying->get_header()), pending_requests(
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
	// This is just a buffer to contain requests from user tasks, so it
	// doesn't need to have a large size.
	user_comp_requests(underlying->get_node_id(), 512)
{
	assert(t == underlying->get_thread());
	ext_allocator = std::unique_ptr<req_ext_allocator>(
			new req_ext_allocator(underlying->get_node_id()));
	req_allocator = std::unique_ptr<request_allocator>(
			new request_allocator(underlying->get_node_id()));
	orig_array_allocator
		= std::unique_ptr<byte_array_allocator_impl<original_req_byte_array> >(
				new byte_array_allocator_impl<original_req_byte_array>(t));
	simp_array_allocator
		= std::unique_ptr<byte_array_allocator_impl<simple_page_byte_array> >(
				new byte_array_allocator_impl<simple_page_byte_array>(t));

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

	if (sched == NULL)
		comp_io_sched = comp_io_scheduler::ptr(
				new default_comp_io_scheduler(this->get_node_id()));
	else
		comp_io_sched = sched;
	comp_io_sched->set_io(this);
}

global_cached_io::~global_cached_io()
{
	cleanup();
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

	assert(npages <= get_block_size());

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
		BOOST_VERIFY(file_id == p->get_file_id());
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
	get_global_cache().mark_dirty_pages(dirty_pages.data(),
			dirty_pages.size(), *underlying);
	return tot;
}

void merge_pages2req(io_request &req, page_cache &cache, size_t block_size)
{
	if (!params.is_cache_large_write())
		return;

	thread_safe_page *p;
	off_t off = req.get_offset();
	off_t forward_off = off + PAGE_SIZE;
	off_t block_off = ROUND(off, block_size * PAGE_SIZE);
	off_t block_end_off = block_off + block_size * PAGE_SIZE;
	page_id_t pg_id(req.get_file_id(), forward_off);
	while (forward_off < block_end_off
			&& (p = (thread_safe_page *) cache.search(pg_id))) {
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
				&& (p = (thread_safe_page *) cache.search(pg_id))) {
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
	assert(req.inside_RAID_block(block_size));
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

	merge_pages2req(req, get_global_cache(), get_block_size());
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
	if (status == IO_FAIL)
		throw io_exception("Can't write dirty page");
}

thread_safe_page *complete_cached_req(const io_request &req, thread_safe_page *p,
		byte_array_allocator &alloc)
{
	if (req.get_req_type() == io_request::BASIC_REQ) {
		int page_off;
		thread_safe_page *ret = NULL;
		char *req_buf;
		int req_size;

		page_off = req.get_offset() - ROUND_PAGE(req.get_offset());
		req_buf = req.get_buf();
		req_size = req.get_size();

		p->lock();
		if (req.get_access_method() == WRITE) {
			memcpy((char *) p->get_data() + page_off, req_buf, req_size);
			if (!p->set_dirty(true)) {
				ret = p;
				ret->inc_ref();
			}
		}
		else 
			/* I assume the data I read never crosses the page boundary */
			memcpy(req_buf, (char *) p->get_data() + page_off, req_size);
		p->unlock();
		p->dec_ref();
		return ret;
	}
	else {
		simple_page_byte_array arr(req, p, alloc);
		user_compute *compute = req.get_compute();
		compute->run(arr);
		return NULL;
	}
}

void global_cached_io::process_cached_reqs()
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
		thread_safe_page *dirty = complete_cached_req(pair.first,
				pair.second, *simp_array_allocator);
		page_cache &cache = get_global_cache();
		if (dirty) {
			cache.mark_dirty_pages(&dirty, 1, *underlying);
			dirty->dec_ref();
		}
		if (!pair.first.is_sync()) {
			req_buf[num_reqs_in_buf] = pair.first;
			reqp_buf[num_reqs_in_buf] = &req_buf[num_reqs_in_buf];
			num_reqs_in_buf++;
			num_async_reqs++;
		}

		if (pair.first.get_req_type() == io_request::USER_COMPUTE) {
			user_compute *compute = pair.first.get_compute();
			comp_io_sched->post_comp_process(compute);
		}

		if (num_reqs_in_buf == REQ_BUF_SIZE) {
			safs::notify_completion(this, reqp_buf, num_reqs_in_buf);
			num_reqs_in_buf = 0;
		}
	}
	// We don't need to notify completion for sync requests.
	// Actually, we don't even need to do anything for sync requests.
	num_completed_areqs.inc(num_async_reqs);
	if (num_reqs_in_buf > 0)
		safs::notify_completion(this, reqp_buf, num_reqs_in_buf);
}

void global_cached_io::process_user_req(
		std::vector<thread_safe_page *> &dirty_pages, io_status *status)
{
	int pg_idx = 0;
	// In max, we use the number of pages in the RAID block.
	thread_safe_page *pages[get_block_size()];
	int num_pages_ready = 0;
	int num_bytes_completed = 0;
	while (!processing_req.is_empty()) {
		thread_safe_page *p;
		page_id_t pg_id = processing_req.get_curr_page_id();
		page_id_t old_id;
		do {
			p = (thread_safe_page *) (get_global_cache().search(pg_id, old_id));
			// If the cache can't evict a page, it's probably because
			// all pages have been referenced. It's likely that we issued
			// too many requests. Let's stop issuing more requests for now.
			if (p == NULL)
				goto end;
		} while (p == NULL);
		processing_req.move_next();
		num_pg_accesses++;

		/* 
		 * If old_off is -1, it means search() didn't evict a page, i.e.,
		 * it's a cache hit.
		 */
		if (old_id.get_offset() == -1) {
			cache_hits++;
			if (p->data_ready())
				num_pages_ready++;
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
			assert(pg_idx <= get_block_size());
			if ((pages[0]->get_offset() + PAGE_SIZE * pg_idx)
					% (get_block_size() * PAGE_SIZE) == 0) {
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
		if (num_pages_ready == processing_req.get_request().get_num_covered_pages()
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
				) < (size_t) get_max_num_pending_ios()
			// TODO the maximal number should be configurable.
			&& num_underlying_pages.get() < 1000) {
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

	get_global_cache().mark_dirty_pages(dirty_pages.data(),
				dirty_pages.size(), *underlying);

	flush_requests();
}

void global_cached_io::access(io_request *requests, int num, io_status *status)
{
	if (num == 0)
		return;

	ASSERT_EQ(get_thread(), thread::get_curr_thread());

	bool syncd = false;
	std::vector<thread_safe_page *> dirty_pages;

	int num_async = 0;
	for (int i = 0; i < num; i++) {
		if (requests[i].get_io() == NULL) {
			requests[i].set_io(this);
			requests[i].set_node_id(this->get_node_id());
		}
		if (!requests[i].is_sync())
			num_async++;
		// The user compute will be referenced by IO requests. I need to
		// increase their references now.
		if (requests[i].get_req_type() == io_request::USER_COMPUTE) {
			user_compute *compute = requests[i].get_compute();
			compute->inc_ref();
		}
	}
	num_issued_areqs.inc(num_async);

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
		else
			num_processed_areqs.inc(1);
		assert(processing_req.is_empty());
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
	get_global_cache().mark_dirty_pages(dirty_pages.data(),
				dirty_pages.size(), *underlying);

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
	else {
		// Process the completed requests served in the cache directly.
		// For one-page requests served by the page cache, access() processes them
		// differently. If the sync request is within a single page and is served
		// by the page cache, it is kept in the cached request queue. We need to
		// call this function to complete processing this request.
		process_cached_reqs();
	}
	// TODO IO may fail, I need to return an error in case it fails.
	status = IO_OK;
	status.set_priv_data(size);
	return status;
}

int global_cached_io::preload(off_t start, long size) {
	if (size > cache_size) {
		fprintf(stderr, "we can't preload data larger than the cache size\n");
		return -1;
	}

	assert(ROUND_PAGE(start) == start);
	for (long offset = start; offset < start + size; offset += PAGE_SIZE) {
		page_id_t pg_id(get_file_id(), ROUND_PAGE(offset));
		page_id_t old_id;
		thread_safe_page *p = (thread_safe_page *) (get_global_cache().search(
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

	do {
		// Process the completed requests served in the cache directly.
		process_cached_reqs();
		// Process completed user requests.
		process_completed_requests();

		// The number of requests generated by user compute.
		int num_new_reqs;
		do {
			// Let's first try to process incomplete computation.
			// It may generate more user requests to the cache and the underlying IO.
			int orig_num = user_comp_requests.get_num_entries();
			comp_io_sched->get_requests(user_comp_requests,
					get_remaining_io_slots());
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
		// Processing user tasks and user requests may generate more completed
		// requests. Let's go back and process these completed reqeusts again.
	} while (!cached_requests.is_empty() || !complete_queue.is_empty());
	comp_io_sched->gc_computes();

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
	// Now we can delete it.
	req_allocator->free(req);
}

/**
 * We wait for at least the specified number of requests to complete.
 */
int global_cached_io::wait4complete(int num_to_complete)
{
	flush_requests();
	size_t prev_completed_areqs = num_completed_areqs.get();
	int prev_pending = num_pending_ios();
	num_to_complete = min(prev_pending, num_to_complete);

	process_all_requests();
	// We should use the number of pending I/O requests to measure the number
	// of completed requests because completed requests may still have
	// incomplete user tasks and we need to take into account the number
	// of incomplete tasks.
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

static bool cross_RAID_block_bound(off_t off, size_t size, size_t block_size)
{
	off_t end = off + size;
	off_t block_begin = ROUND(off, block_size * PAGE_SIZE);
	return end - block_begin > (off_t) (block_size * PAGE_SIZE);
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
			BOOST_VERIFY(p == merged.get_page(idx));
		}
		else
			merged.add_io_buf(buf);
	}
	// We don't need to reference the page for a multi-buf request.
	merged.set_priv(NULL);
	return true;
}

static inline void access(io_interface &underlying, io_request &req)
{
	io_status status;
	underlying.access(&req, 1, &status);
	if (status == IO_FAIL)
		throw io_exception("Fail to issue an I/O request");
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
						req.get_size() + under_req->get_size(), get_block_size())
					// We can't merge the two requests.
					|| !merge_req(req, *under_req)) {
				num_sent++;
				num_pages += req.get_num_bufs();
				safs::access(*underlying, req);

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
		safs::access(*underlying, req);
	}
	underlying->flush_requests();
}

void global_cached_io::wakeup_on_req(original_io_request *req, int status)
{
	assert(req->is_sync());
	assert(req->is_complete());
	get_thread()->activate();
}

}
