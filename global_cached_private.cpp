#include "global_cached_private.h"
#include "container.cpp"

static void __complete_req(io_request *orig, thread_safe_page *p)
{
	int page_off = orig->get_offset() - ROUND_PAGE(orig->get_offset());

	p->lock();
	if (orig->get_access_method() == WRITE) {
		memcpy((char *) p->get_data() + page_off, orig->get_buf(),
				orig->get_size());
		p->set_dirty(true);
	}
	else 
		/* I assume the data I read never crosses the page boundary */
		memcpy(orig->get_buf(), (char *) p->get_data() + page_off,
				orig->get_size());
	p->unlock();
	p->dec_ref();
}

void read_complete(io_request *request)
{
	io_request *orig = (io_request *) request->get_priv();
	thread_safe_page *p = (thread_safe_page *) orig->get_priv();

	p->lock();
	p->set_data_ready(true);
	p->set_io_pending(false);
	io_request *req = p->get_io_req();
	assert(req);
	int access_method = req->get_access_method();
	assert(access_method == READ);
	// Now data is ready, no threads will add more read requests.
	io_request *old = p->reset_reqs();
	p->unlock();

	while (old) {
		io_request *next = old->get_next_req();
		__complete_req(old, p);
		io_interface *io = old->get_io();
		if (io->get_callback())
			io->get_callback()->invoke(old);
		// Now we can delete it.
		delete old;
		old = next;
	}
}

void reissue_read(io_request *request)
{
	io_request *orig = (io_request *) request->get_priv();
	thread_safe_page *p = (thread_safe_page *) orig->get_priv();

	p->lock();
	p->set_old_dirty(false);
	p->set_io_pending(false);
	io_request *req = p->get_io_req();
	assert(req);
	io_request *old = p->reset_reqs();
	p->unlock();
	global_cached_io *io = static_cast<global_cached_io *>(
			old->get_io());
	// We can't invoke write() here because it may block the thread.
	// Instead, we queue the request, so it will be issue to
	// the device by the user thread.
	io->queue_request(old);
	// These requests can't be deleted yet.
	// They will be deleted when these write requests are finally served.
}

void reissue_write(io_request *request)
{
	io_request *orig = (io_request *) request->get_priv();
	thread_safe_page *p = (thread_safe_page *) orig->get_priv();

	p->lock();
	// If we write data to part of a page, we need to first read
	// the entire page to memory first.
	if (request->get_access_method() == READ) {
		p->set_data_ready(true);
		p->set_io_pending(false);
	}
	// We just evict a page with dirty data and write the original
	// dirty data in the page to a file.
	else {
		p->set_old_dirty(false);
		p->set_io_pending(false);
	}
	io_request *req = p->get_io_req();
	assert(req);
	io_request *old = p->reset_reqs();
	p->unlock();

	if (request->get_access_method() == READ) {
		while (old) {
			io_request *next = old->get_next_req();
			// We are ready to write data to the page.
			__complete_req(old, p);
			io_interface *io = old->get_io();
			if (io->get_callback())
				io->get_callback()->invoke(old);

			delete old;
			old = next;
		}
	}
	else {
		global_cached_io *io = static_cast<global_cached_io *>(
				old->get_io());
		// We can't invoke write() here because it may block the thread.
		// Instead, we queue the request, so it will be issue to
		// the device by the user thread.
		io->queue_request(old);
		// These requests can't be deleted yet.
		// They will be deleted when these write requests are finally served.
	}
}

class access_page_callback: public callback
{
public:
	int invoke(io_request *request) {
		io_request *orig = (io_request *) request->get_priv();
		orig->get_cb()(request);
		return 0;
	}
};

global_cached_io::global_cached_io(io_interface *underlying): pending_requests(
		INIT_GCACHE_PENDING_SIZE)
{
	this->underlying = underlying;
	num_waits = 0;
	cache_size = 0;
	cb = NULL;
	cache_hits = 0;
}

global_cached_io::global_cached_io(io_interface *underlying, long cache_size,
		int cache_type): pending_requests(INIT_GCACHE_PENDING_SIZE)
{
	cb = NULL;
	cache_hits = 0;
	this->underlying = underlying;
	underlying->set_callback(new access_page_callback());
	num_waits = 0;
	this->cache_size = cache_size;
	if (global_cache == NULL) {
		global_cache = create_cache(cache_type, cache_size);
	}
}

ssize_t global_cached_io::__write(io_request *orig, thread_safe_page *p)
{
	ssize_t ret = 0;
	orig->set_priv((void *) p);
	orig->set_cb(reissue_write);
	p->lock();
	if (!p->data_ready()) {
		if(!p->is_io_pending()) {
			p->set_io_pending(true);
			assert(!p->is_dirty());

			// We are going to write to part of a page, therefore,
			// we need to first read the page.
			if (orig->get_size() < PAGE_SIZE) {
				io_request read_req((char *) p->get_data(),
						ROUND_PAGE(orig->get_offset()), PAGE_SIZE, READ,
						underlying, (char *) orig);
				p->add_req(orig);
				p->unlock();
				ret = underlying->access(&read_req, 1);
				if (ret < 0) {
					perror("read");
					exit(1);
				}
			}
			else {
				p->unlock();
				ret = PAGE_SIZE;
				__complete_req(orig, p);
				if (get_callback())
					get_callback()->invoke(orig);
				delete orig;
			}
		}
		else {
			// If there is an IO pending, it means a read request
			// has been issuded. It can't be a write request, otherwise,
			// the data in the page will be ready.
			p->add_req(orig);
			p->unlock();
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
		if (get_callback())
			get_callback()->invoke(orig);
		ret = orig->get_size();
		delete orig;
#ifdef STATISTICS
		cache_hits++;
#endif
	}
	return ret;
}

ssize_t global_cached_io::write(io_request &req, thread_safe_page *p)
{
	io_request *orig = new io_request(req);
	return __write(orig, p);
}

ssize_t global_cached_io::__read(io_request *orig, thread_safe_page *p)
{
	ssize_t ret = 0;
	orig->set_priv((void *) p);
	orig->set_cb(read_complete);
	p->lock();
	if (!p->data_ready()) {
		if(!p->is_io_pending()) {
			p->set_io_pending(true);
			assert(!p->is_dirty());

			p->add_req(orig);
			io_request req((char *) p->get_data(), p->get_offset(),
					/*
					 * it will notify the underlying IO,
					 * which then notifies global_cached_io.
					 */
					PAGE_SIZE, READ, underlying, (void *) orig);
			p->unlock();
			ret = underlying->access(&req, 1);
			if (ret < 0) {
				perror("read");
				exit(1);
			}
		}
		else {
			p->add_req(orig);
			p->unlock();
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

ssize_t global_cached_io::read(io_request &req, thread_safe_page *p)
{
	io_request *orig = new io_request(req);
	return __read(orig, p);
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
			thread_safe_page *p = static_cast<thread_safe_page *>(
					req->get_priv());
			assert(p);
			assert(!p->is_old_dirty());
			while (req) {
				io_request *next = req->get_next_req();
				assert(req->get_priv() == p);
				req->set_next_req(NULL);
				if (req->get_access_method() == WRITE)
					__write(req, p);
				else
					__read(req, p);
				req = next;
			}
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
		off_t old_off = -1;
		off_t offset = requests[i].get_offset();

		// TODO a request should be within a page.
		assert(offset - ROUND_PAGE(offset) + requests[i].get_size() <= PAGE_SIZE);
		thread_safe_page *p = (thread_safe_page *) (get_global_cache()
				->search(ROUND_PAGE(offset), old_off));
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
			/* The page is evicted in this thread */
			if (old_off != ROUND_PAGE(offset) && old_off != -1) {
				/*
				 * Only one thread can come here because only one thread
				 * can evict the dirty page and the thread gets its old
				 * offset, and only this thread can write back the old
				 * dirty page.
				 */
				io_request *orig = new io_request(requests[i]);
				orig->set_priv((void *) p);
				if (orig->get_access_method() == READ)
					orig->set_cb(reissue_read);
				else
					orig->set_cb(reissue_write);
				p->lock();
				p->add_req(orig);
				p->set_io_pending(true);
				io_request req((char *) p->get_data(),
						old_off, PAGE_SIZE, WRITE, underlying, (void *) orig);
				p->unlock();
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
				io_request *orig = new io_request(requests[i]);
				orig->set_priv((void *) p);
				orig->set_cb(reissue_write);
				p->lock();
				if (p->is_old_dirty()) {
					p->add_req(orig);
					p->unlock();
					// the request has been added to the page, when the old dirty
					// data is written back to the file, the write request will be
					// reissued to the file.
					continue;
				}
				else {
					p->unlock();
					delete orig;
				}
			}
		}

		if (requests[i].get_access_method() == WRITE)
			write(requests[i], p);
		else
			read(requests[i], p);
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
				perror("read");
				return ret;
			}
			p->set_io_pending(false);
			p->set_data_ready(true);
		}
	}
	/* close the file as it will be opened again in the real workload. */
	underlying->cleanup();
	return 0;
}

page_cache *global_cached_io::global_cache;
