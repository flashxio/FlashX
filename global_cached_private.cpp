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
	global_cached_io *io;
public:
	read_page_callback(global_cached_io *io) {
		this->io = io;
	}

	int invoke(io_request *request) {
		io_request *orig = (io_request *) request->get_priv();
		thread_safe_page *p = (thread_safe_page *) orig->get_priv();
		p->set_data_ready(true);
		p->set_io_pending(false);
		__complete_req(orig, p);
		if (io->get_cb())
			io->get_cb()->invoke(orig);
		// TODO use my own deallocator.
		delete orig;
		return 0;
	}
};

global_cached_io::global_cached_io(io_interface *underlying,
		long cache_size, int cache_type, memory_manager *manager) {
	cb = NULL;
	cache_hits = 0;
	this->underlying = underlying;
	underlying->set_callback(new read_page_callback(this));
	num_waits = 0;
	this->cache_size = cache_size;
	if (global_cache == NULL) {
		page::allocate_cache(cache_size);
		global_cache = create_cache(cache_type, cache_size, manager);
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

		if (!p->data_ready()) {
			/*
			 * TODO if two threads want to access the same page,
			 * both of them will come here and the page will be
			 * read twice from the file.
			 * I can't use the IO pending bit because a request
			 * can't be served immediately when it's issued.
			 * I just hope the case that two threads access 
			 * the same page doesn't happen very often.
			 */
			/*
			 * No other threads set the page dirty at this moment,
			 * because if a thread can set the page dirty,
			 * it means the page isn't dirty and already has data
			 * ready at the first place.
			 */
			if (p->is_dirty())
				p->wait_cleaned();

			// TODO I need to my own allocator for this.
			io_request *orig = new io_request(requests[i]);
			orig->set_priv((void *) p);

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
			__complete_req(&requests[i], p);
#ifdef STATISTICS
			cache_hits++;
#endif
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
