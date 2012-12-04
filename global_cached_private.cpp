#include "global_cached_private.h"

static void __complete_req(io_request *orig, thread_safe_page *p)
{
	int page_off = orig->get_offset() - ROUND_PAGE(orig->get_offset());

	p->lock();
	if (orig->get_access_method() == WRITE) {
		unsigned long *l = (unsigned long *) ((char *) p->get_data()
				+ page_off);
		unsigned long start = orig->get_offset() / sizeof(long);
		if (*l != start)
			printf("write: start: %ld, l: %ld, offset: %ld\n",
					start, *l, p->get_offset() / sizeof(long));
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

class read_page_callback: public callback
{
public:
	int invoke(io_request *request) {
		io_request *orig = (io_request *) request->get_priv();
		thread_safe_page *p = (thread_safe_page *) orig->get_priv();
		std::vector<io_request> reqs;

		p->lock();
		p->set_data_ready(true);
		p->set_io_pending(false);
		io_request *req = p->get_io_req();
		assert(req);
		// TODO vector may allocate memory and block the thread.
		// I should be careful of that.
		while (req) {
			reqs.push_back(*req);
			req = req->get_next_req();
		}
		p->unlock();
		// At this point, data is ready, so any threads that want to add
		// requests to the page will fail.
		p->reset_reqs();

		for (unsigned i = 0; i < reqs.size(); i++) {
			__complete_req(&reqs[i], p);
			io_interface *io = reqs[i].get_io();
			if (io->get_callback())
				io->get_callback()->invoke(&reqs[i]);
		}
		return 0;
	}
};

global_cached_io::global_cached_io(io_interface *underlying,
		long cache_size, int cache_type) {
	cb = NULL;
	cache_hits = 0;
	this->underlying = underlying;
	underlying->set_callback(new read_page_callback());
	num_waits = 0;
	this->cache_size = cache_size;
	if (global_cache == NULL) {
		global_cache = create_cache(cache_type, cache_size);
	}
}

ssize_t global_cached_io::access(io_request *requests, int num)
{
	for (int i = 0; i < num; i++) {
		ssize_t ret;
		off_t old_off = -1;
		off_t offset = requests[i].get_offset();
		thread_safe_page *p = (thread_safe_page *) (get_global_cache()
				->search(ROUND_PAGE(offset), old_off));

		/*
		 * the page isn't in the cache,
		 * so the cache evict a page and return it to us.
		 */
		if (old_off != ROUND_PAGE(offset) && old_off != -1) {
			// TODO handle dirty pages later.
		}

		p->lock();
		if (!p->data_ready()) {
			if(!p->is_io_pending()) {
				p->set_io_pending(true);
				assert(!p->is_dirty());
				p->unlock();

				// If the thread can't add the request as the first IO request
				// of the page, it doesn't need to issue a IO request to the file.
				int status;
				io_request *orig = p->add_first_io_req(requests[i], status);
				if (status == ADD_REQ_SUCCESS) {
					io_request req((char *) p->get_data(),
							ROUND_PAGE(offset), PAGE_SIZE, READ,
							/*
							 * it will notify the underlying IO,
							 * which then notifies global_cached_io.
							 */
							underlying, (void *) orig);
					ret = underlying->access(&req, 1);
					if (ret < 0) {
						perror("read");
						exit(1);
					}
				}
				else {
#ifdef STATISTICS
					cache_hits++;
#endif
					// If data is ready, the request won't be added to the page.
					// so just serve the data.
					if (status == ADD_REQ_DATA_READY)
						goto serve_data;
				}
			}
			else {
				p->unlock();
#ifdef STATISTICS
				cache_hits++;
#endif
				// An IO request can be added only when the page is still
				// in IO pending state. Otherwise, it returns NULL.
				io_request *orig = p->add_io_req(requests[i]);
				// the page isn't in IO pending state any more.
				if (orig == NULL)
					goto serve_data;
			}
		}
		else {
			// If the data in the page is ready, we don't need to change any state
			// of the page and just read data.
			p->unlock();
#ifdef STATISTICS
			cache_hits++;
#endif
serve_data:
			__complete_req(&requests[i], p);
			if (get_callback())
				get_callback()->invoke(&requests[i]);
		}
	}
	return 0;
}

ssize_t global_cached_io::access(char *buf, off_t offset,
		ssize_t size, int access_method)
{
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
