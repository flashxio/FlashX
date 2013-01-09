#include "global_cached_private.h"
#include "container.cpp"

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

/**
 * Extract a request from the input request.
 * The extract request is within the range [off, off + npages * PAGE_SIZE),
 * where off is aligned with PAGE_SIZE.
 */
static void extract_pages(const io_request &req, off_t off, int npages,
		io_request &extracted)
{
	off_t req_off;
	char *req_buf;
	ssize_t req_size;
	assert((off & (PAGE_SIZE - 1)) == 0);
	// this is the first page in the request.
	if (off == ROUND_PAGE(req.get_offset())) {
		req_off = req.get_offset();
		req_buf = req.get_buf();
		// the remaining size in the page.
		req_size = PAGE_SIZE * npages - (req_off - off);
		if (req_size > req.get_size())
			req_size = req.get_size();
	}
	else {
		req_off = off;
		/* 
		 * We can't be sure if the request buffer is aligned
		 * with the page size.
		 */
		req_buf = req.get_buf() + (off - req.get_offset());
		ssize_t remaining = req.get_size() - (off - req.get_offset());
		req_size = remaining > PAGE_SIZE * npages ? PAGE_SIZE
			* npages : remaining;
	}
	extracted.init(req_buf, req_off, req_size, req.get_access_method(),
			req.get_io());
}

static void __complete_req(io_request *orig, thread_safe_page *p)
{
	io_request extracted;
	extract_pages(*orig, p->get_offset(), 1, extracted);
	int page_off = extracted.get_offset() - ROUND_PAGE(extracted.get_offset());

	p->lock();
	if (orig->get_access_method() == WRITE) {
		memcpy((char *) p->get_data() + page_off, extracted.get_buf(),
				extracted.get_size());
		p->set_dirty(true);
	}
	else 
		/* I assume the data I read never crosses the page boundary */
		memcpy(extracted.get_buf(), (char *) p->get_data() + page_off,
				extracted.get_size());
	p->unlock();
	p->dec_ref();
}

static void __complete_req_unlocked(io_request *orig, thread_safe_page *p)
{
	io_request extracted;
	extract_pages(*orig, p->get_offset(), 1, extracted);
	int page_off = extracted.get_offset() - ROUND_PAGE(extracted.get_offset());

	if (orig->get_access_method() == WRITE) {
		memcpy((char *) p->get_data() + page_off, extracted.get_buf(),
				extracted.get_size());
		p->set_dirty(true);
	}
	else 
		/* I assume the data I read never crosses the page boundary */
		memcpy(extracted.get_buf(), (char *) p->get_data() + page_off,
				extracted.get_size());
	p->dec_ref();
}

class access_page_callback: public callback
{
	page_cache *cache;
public:
	access_page_callback(page_cache *cache) {
		this->cache = cache;
	}
	int invoke(io_request *request);
	int multibuf_invoke(io_request *request);
};

void finalize_partial_request(io_request &partial, io_request *orig)
{
	io_interface *io = partial.get_io();
	orig->inc_complete_count();
	if (orig->complete_size(partial.get_size())) {
		if (io->get_callback())
			io->get_callback()->invoke(orig);
		orig->dec_complete_count();
		orig->wait4unref();
		// Now we can delete it.
		delete orig;
	}
	else
		orig->dec_complete_count();
}

/**
 * This method is to finalize the request. The processing of the request
 * ends here.
 */
void finalize_request(io_request &req)
{
	// It's possible that the request is just a partial request.
	io_interface *io = req.get_io();
	if (req.is_partial()) {
		io_request *original = req.get_orig();
		assert(original);
		assert(original->get_orig() == NULL);
		original->inc_complete_count();
		if (original->complete_size(req.get_size())) {
			if (io->get_callback())
				io->get_callback()->invoke(original);
			original->dec_complete_count();
			original->wait4unref();
			delete original;
		}
		else
			original->dec_complete_count();
	}
	else {
		assert(req.get_orig() == NULL);
		if (io->get_callback())
			io->get_callback()->invoke(&req);
	}
}

int access_page_callback::multibuf_invoke(io_request *request)
{
	assert(request->get_access_method() == READ);
	io_request *orig = request->get_orig();
	assert(orig->get_orig() == NULL);
	assert(orig->get_num_bufs() == 1);
	/*
	 * Right now the global cache only support normal access().
	 */
	io_request *pending_reqs[request->get_num_bufs()];
	thread_safe_page *pages[request->get_num_bufs()];
	off_t off = request->get_offset();
	off_t old_off = -1;
	for (int i = 0; i < request->get_num_bufs(); i++) {
		thread_safe_page *p = (thread_safe_page *) cache->search(off, old_off);
		pages[i] = p;
		// The page must exist in the cache originally.
		assert(old_off == -1);
		p->lock();
		assert(p->get_data() == request->get_buf(i));
		assert(p->is_io_pending());
		assert(!p->data_ready());
		p->set_data_ready(true);
		p->set_io_pending(false);
		pending_reqs[i] = p->reset_reqs();
		__complete_req_unlocked(orig, p);
		p->unlock();
		off += PAGE_SIZE;
	}
	/*
	 * For a multi-buf request, the private data actually points to
	 * the very original request.
	 */
	io_request partial;
	extract_pages(*orig, request->get_offset(), request->get_num_bufs(), partial);
	finalize_partial_request(partial, orig);

	/*
	 * Now we should start to deal with all requests pending to pages
	 * All of these requests should be single buffer requests.
	 */
	for (int i = 0; i < request->get_num_bufs(); i++) {
		io_request *old = pending_reqs[i];
		thread_safe_page *p = pages[i];
		while (old) {
			io_request *next = old->get_next_req();
			__complete_req(old, p);
			finalize_request(*old);
			// Now we can delete it.
			delete old;
			old = next;
		}
		p->dec_ref();
	}

	return -1;
}

int access_page_callback::invoke(io_request *request)
{
	if (request->get_num_bufs() > 1)
		return multibuf_invoke(request);

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
		p->set_old_dirty(false);
	}
	p->set_io_pending(false);
	io_request *old = p->reset_reqs();
	bool data_ready = p->data_ready();
	p->unlock();

	// If the data on the page is ready, it won't become unready.
	// The only place where data is set unready is where the page is evicted.
	// Since we have a reference of the page, it won't be evicted.
	// When data is ready, we can execuate any operations on the page.
	if (data_ready) {
		/* The request should contain the very original request. */
		io_request *orig = request->get_orig();
		assert(orig->get_orig() == NULL);
		__complete_req(orig, p);
		io_request partial;
		extract_pages(*orig, request->get_offset(), request->get_num_bufs(), partial);
		finalize_partial_request(partial, orig);

		int num = 0;
		global_cached_io *cached_io = NULL;
		while (old) {
			/*
			 * It should be guaranteed that there isn't a multi-buf request
			 * in the queue. Because if a page is in IO pending, we won't
			 * issue a multi-buf request for the page.
			 */
			io_request *next = old->get_next_req();
			assert(old->get_num_bufs() == 1);
			__complete_req(old, p);
			io_interface *io = old->get_io();
			cached_io = static_cast<global_cached_io *>(io);

			finalize_request(*old);
			// Now we can delete it.
			delete old;
			old = next;
			num++;
		}
		if (cached_io)
			cached_io->dec_pending(num);
	}
	else {
		io_request *orig = request->get_orig();
		global_cached_io *io = static_cast<global_cached_io *>(
				orig->get_io());
		// We can't invoke write() here because it may block the thread.
		// Instead, we queue the request, so it will be issue to
		// the device by the user thread.
		assert(orig->get_next_req() == NULL);
		orig->set_next_req(old);
		io->queue_request(orig);
		// These requests can't be deleted yet.
		// They will be deleted when these write requests are finally served.
	}
	return 0;
}

global_cached_io::global_cached_io(io_interface *underlying): pending_requests(
		INIT_GCACHE_PENDING_SIZE)
{
	this->underlying = underlying;
	num_waits = 0;
	cache_size = 0;
	cb = NULL;
	cache_hits = 0;
	num_pending_reqs = 0;
}

global_cached_io::global_cached_io(io_interface *underlying, long cache_size,
		int cache_type): pending_requests(INIT_GCACHE_PENDING_SIZE)
{
	cb = NULL;
	cache_hits = 0;
	this->underlying = underlying;
	num_waits = 0;
	num_pending_reqs = 0;
	this->cache_size = cache_size;
	if (global_cache == NULL) {
		global_cache = create_cache(cache_type, cache_size);
	}
	underlying->set_callback(new access_page_callback(global_cache));
}

/**
 * a write request only covers the memory within one page.
 */
ssize_t global_cached_io::__write(io_request *orig, thread_safe_page *p)
{
	ssize_t ret = 0;
//	orig->set_priv((void *) p);
	p->lock();
	assert(!p->is_old_dirty());
	if (!p->data_ready()) {
		if(!p->is_io_pending()) {
			assert(!p->is_dirty());

			// We are going to write to part of a page, therefore,
			// we need to first read the page.
			if (orig->get_size() < PAGE_SIZE) {
				off_t off = orig->get_offset();
				io_request *real_orig = orig->get_orig();
				// If the request doesn't have a private data, it is the real
				// original request.
				if (real_orig == NULL)
					real_orig = orig;
				else
					// `orig' is just part of the original request.
					// we don't need it any more.
					delete orig;
				assert(real_orig->get_orig() == NULL);
				io_request read_req((char *) p->get_data(),
						ROUND_PAGE(off), PAGE_SIZE, READ,
						underlying, real_orig, p);
				p->set_io_pending(true);
				p->unlock();
				inc_pending(1);
				ret = underlying->access(&read_req, 1);
				if (ret < 0) {
					perror("read");
					exit(1);
				}
			}
			else {
				// This is an optimization. If we can overwrite the entire page,
				// we don't need to read the page first. However, we have to
				// make sure data is written to a page without anyone else
				// having IO operations on it.
				__complete_req_unlocked(orig, p);
				p->set_data_ready(true);
				p->unlock();
				ret = PAGE_SIZE;
				finalize_request(*orig);
				// Now we can delete it.
				delete orig;
			}
		}
		else {
			// If there is an IO pending, it means a read request
			// has been issuded. It can't be a write request, otherwise,
			// the data in the page will be ready.
			orig->set_priv(p);
			assert(orig->get_access_method() == WRITE);
			p->add_req(orig);
			p->unlock();
			inc_pending(1);
#ifdef STATISTICS
			cache_hits++;
#endif
		}
	}
	else {
		// The data in the page is ready. We can write data to the page directly.
		//
		// If data is ready, there shouldn't be an IO pending.
		// In other words, if the thread for writing dirty pages is writing
		// a page, the page will be referenced and therefore, can't be returned
		// from the cache.
		assert(!p->is_io_pending());
		p->unlock();

		__complete_req(orig, p);
		ret = orig->get_size();
		finalize_request(*orig);
		// Now we can delete it.
		delete orig;
#ifdef STATISTICS
		cache_hits++;
#endif
	}
	return ret;
}

ssize_t global_cached_io::write(io_request &req, thread_safe_page *p)
{
	io_request *orig = req.get_orig();
	if (orig->get_size() == req.get_size())
		return __write(orig, p);
	else {
		io_request *partial_orig = new io_request(req);
		partial_orig->set_partial(true);
		return __write(partial_orig, p);
	}
}

ssize_t global_cached_io::__read(io_request *orig, thread_safe_page *p)
{
	ssize_t ret = 0;
//	orig->set_priv((void *) p);
	p->lock();
	if (!p->data_ready()) {
		if(!p->is_io_pending()) {
			p->set_io_pending(true);
			assert(!p->is_dirty());

			io_request req((char *) p->get_data(), p->get_offset(),
					/*
					 * it will notify the underlying IO,
					 * which then notifies global_cached_io.
					 */
					PAGE_SIZE, READ, underlying, orig, p);
			p->unlock();
			inc_pending(1);
			assert(orig->get_orig() == NULL);
			ret = underlying->access(&req, 1);
			if (ret < 0) {
				perror("read");
				exit(1);
			}
		}
		else {
			orig->set_priv(p);
			assert(orig->get_access_method() == READ);
			p->add_req(orig);
			p->unlock();
			inc_pending(1);
#ifdef STATISTICS
			cache_hits++;
#endif
		}
	}
	else {
		// If the data in the page is ready, we don't need to change any state
		// of the page and just read data.
		p->unlock();
#ifdef STATISTICS
		cache_hits++;
#endif
		ret = orig->get_size();
		__complete_req(orig, p);
		if (get_callback())
			get_callback()->invoke(orig);
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
		int npages)
{
	ssize_t ret = 0;

	assert(npages <= MAX_NUM_IOVECS);
	io_request *orig = req.get_orig();
	assert(orig->get_orig() == NULL);
	io_request multibuf_req(-1, underlying, req.get_access_method(), orig);

	/*
	 * The pages in `pages' should be sorted with their offsets.
	 * We are going to grab multiple locks below. As long as we always
	 * lock pages in the order of page offset, there won't be deadlock.
	 */
	for (int i = 0; i < npages; i++) {
		thread_safe_page *p = pages[i];
again:
		p->lock();
		if (!p->data_ready() && !p->is_io_pending()) {
			p->set_io_pending(true);
			assert(!p->is_dirty());
			if (multibuf_req.is_empty())
				multibuf_req.set_offset(p->get_offset());
			/* We don't need to worry buffer overflow here. */
			multibuf_req.add_buf((char *) p->get_data(), PAGE_SIZE);
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
				multibuf_req.clear();
				goto again;
			}
			else {
				/*
				 * All pending requests on a page have to be a single-buf request.
				 * Furthermore, the pending requests must only cover one page.
				 */
				// TODO I shouldn't allocate memory within locks.
				io_request *partial_orig = new io_request();
				extract_pages(*orig, p->get_offset(), 1, *partial_orig);
				partial_orig->set_partial(true);
				partial_orig->set_orig(orig);
				partial_orig->set_priv(p);
				p->add_req(partial_orig);
				p->unlock();
				inc_pending(1);
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
				multibuf_req.clear();
			}
#ifdef STATISTICS
			cache_hits++;
#endif
			io_request complete_partial;
			extract_pages(*orig, p->get_offset(), 1, complete_partial);
			ret += complete_partial.get_size();
			__complete_req(&complete_partial, p);
			finalize_partial_request(complete_partial, orig);
		}
	}
	if (!multibuf_req.is_empty()) {
		underlying->access(&multibuf_req, 1);
	}
	return ret;
}

int global_cached_io::handle_pending_requests()
{
	int tot = 0;
	while (!pending_requests.is_empty()) {
		io_request *reqs[MAX_FETCH_REQS];
		int num = pending_requests.fetch(reqs, MAX_FETCH_REQS);
		for (int i = 0; i < num; i++) {
			// It may be the head of a request list. All requests
			// in the list should point to the same page.
			io_request *req = reqs[i];
			thread_safe_page *p = (thread_safe_page *) req->get_priv();
			assert(p);
			assert(!p->is_old_dirty());
			int num_same = 0;
			while (req) {
				/**
				 * Right now all pending requests are writes.
				 * All writes are single-buf requests.
				 */
				assert(req->get_num_bufs() == 1);
				io_request *next = (io_request *) req->get_next_req();
				assert(req->get_priv() == p);
				req->set_next_req(NULL);
				if (req->get_access_method() == WRITE)
					__write(req, p);
				else
					__read(req, p);
				req = next;
				num_same++;
			}
			// These requests have been counted as pending requests.
			// When a request is passed to __write(), if it is sent
			// to the device, it will be double counted as pending.
			// If not, the request has been served. In either case,
			// we need to subtract these requests.
			dec_pending(num_same);
		}
		tot += num;
	}
	return tot;
}

ssize_t global_cached_io::access(io_request *requests, int num)
{
	if (!pending_requests.is_empty()) {
		handle_pending_requests();
	}

	for (int i = 0; i < num; i++) {
		ssize_t ret;
		off_t offset = requests[i].get_offset();
		int size = requests[i].get_size();
		off_t begin_pg_offset = ROUND_PAGE(offset);
		off_t end_pg_offset = ROUNDUP_PAGE(offset + size);
		thread_safe_page *pages[MAX_NUM_IOVECS];
		// TODO right now it only supports single-buf requests.
		assert(requests[i].get_num_bufs() == 1);
		io_request *orig = new io_request(requests[i]);

		int pg_idx = 0;
		for (off_t tmp_off = begin_pg_offset; tmp_off < end_pg_offset;
				tmp_off += PAGE_SIZE) {
			off_t old_off = -1;
			thread_safe_page *p = (thread_safe_page *) (get_global_cache()
					->search(tmp_off, old_off));
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
				/*
				 * We got a few contiguous pages for read, so we should split
				 * the request and issue reads for the contiguous pages first.
				 * We always break write requests into pages, so it has to be
				 * read requests.
				 */
				if (pg_idx) {
					io_request req;
					extract_pages(*orig, pages[0]->get_offset(), pg_idx, req);
					req.set_orig(orig);
					req.set_partial(orig->get_size() > req.get_size());
					read(req, pages, pg_idx);
					pg_idx = 0;
				}

				// Extract the partial access.
				io_request *orig1;
				// If the request accesses more than one page.
				if (end_pg_offset - begin_pg_offset > PAGE_SIZE) {
					orig1 = new io_request();
					extract_pages(*orig, tmp_off, 1, *orig1);
					orig1->set_orig(orig);
					orig1->set_priv(p);
					assert(orig->get_size() > orig1->get_size());
					orig1->set_partial(true);
				}
				else {
					orig1 = orig;
					orig1->set_priv(p);
				}

				assert(!p->data_ready());
				/* The page is evicted in this thread */
				if (old_off != ROUND_PAGE(offset) && old_off != -1) {
					/*
					 * Only one thread can come here because only one thread
					 * can evict the dirty page and the thread gets its old
					 * offset, and only this thread can write back the old
					 * dirty page.
					 */
					p->lock();
					assert(!p->is_io_pending());
					p->set_io_pending(true);
					io_request req((char *) p->get_data(), old_off,
							PAGE_SIZE, WRITE, underlying, orig1, p);
					p->unlock();
					inc_pending(1);
					ret = underlying->access(&req, 1);
					if (ret < 0) {
						perror("read");
						abort();
					}
					continue;
				}
				else {
					// At this moment, the page is being written back to the file
					// by another thread. We should queue the request to tht page,
					// so when the dirty page completes writing back, we can proceed
					// writing.
					p->lock();
					if (p->is_old_dirty()) {
						p->add_req(orig1);
						p->unlock();
						inc_pending(1);
						// the request has been added to the page, when the old dirty
						// data is written back to the file, the write request will be
						// reissued to the file.
						continue;
					}
					else {
						p->unlock();
						if (orig1 != orig)
							delete orig1;
					}
				}
			}

			/*
			 * Large access only makes sense for reading. As large writes
			 * essentially overwrite entire pages in the memory, so we may
			 * only need to read the first and the last pages.
			 * TODO It may be interesting to change our dirty page flushing
			 * algorithm so we can always write a large chunk of data in one
			 * go.
			 */
			if (orig->get_access_method() == WRITE) {
				/* We need to extract a page from the request. */
				io_request req;
				extract_pages(*orig, tmp_off, 1, req);
				req.set_orig(orig);
				req.set_partial(orig->get_size() > req.get_size());
				write(req, p);
			}
			else {
				pages[pg_idx++] = p;
				if (pg_idx == MAX_NUM_IOVECS) {
					io_request req;
					extract_pages(*orig, pages[0]->get_offset(), pg_idx, req);
					req.set_orig(orig);
					req.set_partial(orig->get_size() > req.get_size());
					read(req, pages, pg_idx);
					pg_idx = 0;
				}
			}
		}
		/*
		 * The only reason that pg_idx > 0 is that there is a large read request.
		 */
		if (pg_idx) {
			io_request req;
			extract_pages(*orig, pages[0]->get_offset(), pg_idx, req);
			req.set_orig(orig);
			req.set_partial(orig->get_size() > req.get_size());
			read(req, pages, pg_idx);
		}
	}
	return 0;
}

ssize_t global_cached_io::access(char *buf, off_t offset,
		ssize_t size, int access_method)
{
	assert(access_method == READ);

	ssize_t ret;
	off_t old_off = -1;
	thread_safe_page *p = (thread_safe_page *) (get_global_cache()
			->search(ROUND_PAGE(offset), old_off));

	/*
	 * the page isn't in the cache,
	 * so the cache evict a page and return it to us.
	 */
	if (old_off != ROUND_PAGE(offset) && old_off != -1) {
		/* 
		 * if the new page we get is dirty,
		 * we need to write its data back to the file
		 * before we can put data in the page. 
		 * Therefore, the data ready flag is definitely
		 * not set yet.
		 */
		if (p->is_dirty()) {
			unsigned long *l = (unsigned long *) p->get_data();
			unsigned long start = old_off / sizeof(long);
			if (*l != start)
				printf("start: %ld, l: %ld\n", start, *l);
			underlying->access((char *) p->get_data(),
					old_off, PAGE_SIZE, WRITE);
			p->set_dirty(false);
		}
	}

	if (!p->data_ready()) {
		/* if the page isn't io pending, set it io pending, and return
		 * original result. otherwise, just return the original value.
		 *
		 * This is an ugly hack, but with this atomic operation,
		 * I can avoid using locks.
		 */
		if(!p->test_and_set_io_pending()) {
			/* 
			 * Because of the atomic operation, it's guaranteed
			 * that only one thread can enter here.
			 * If other threads have reference to the page,
			 * they must be waiting for its data to be ready.
			 */
			/*
			 * It's possible that two threads go through here
			 * sequentially. For example, the second thread already
			 * sees the data isn't ready, but find io pending isn't
			 * set when the first thread resets io pending.
			 *
			 * However, in any case, when the second thread comes here,
			 * the data is already ready and the second thread should
			 * be able to see the data is ready when it comes here.
			 */
			if (!p->data_ready()) {
				/*
				 * No other threads set the page dirty at this moment,
				 * because if a thread can set the page dirty,
				 * it means the page isn't dirty and already has data
				 * ready at the first place.
				 */
				if (p->is_dirty())
					p->wait_cleaned();
				ret = underlying->access((char *) p->get_data(),
						ROUND_PAGE(offset), PAGE_SIZE, READ);
				if (ret < 0) {
					perror("read");
					exit(1);
				}
			}
			p->set_data_ready(true);
			p->set_io_pending(false);
		}
		else {
			num_waits++;
			// TODO this is a problem. It takes a long time to get data
			// from real storage devices. If we use busy waiting, 
			// we just waste computation time.
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

int global_cached_io::preload(off_t start, long size) {
	if (size > cache_size) {
		fprintf(stderr, "we can't preload data larger than the cache size\n");
		exit(1);
	}

	/* open the file. It's a hack, but it works for now. */
	underlying->init();

	assert(ROUND_PAGE(start) == start);
	for (long offset = start; offset < start + size; offset += PAGE_SIZE) {
		off_t old_off = -1;
		thread_safe_page *p = (thread_safe_page *) (get_global_cache()->search(ROUND_PAGE(offset), old_off));
		if (!p->data_ready()) {
			ssize_t ret = underlying->access((char *) p->get_data(),
					ROUND_PAGE(offset), PAGE_SIZE, READ);
			if (ret < 0) {
				p->dec_ref();
				perror("read");
				return ret;
			}
			p->set_io_pending(false);
			p->set_data_ready(true);
		}
		p->dec_ref();
	}
	/* close the file as it will be opened again in the real workload. */
	underlying->cleanup();
	return 0;
}

page_cache *global_cached_io::global_cache;
