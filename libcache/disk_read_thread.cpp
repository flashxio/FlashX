#include "cache.h"
#include "disk_read_thread.h"
#include "parameters.h"
#include "aio_private.h"

disk_read_thread::disk_read_thread(const logical_file_partition &_partition,
		const std::tr1::unordered_map<int, aio_complete_thread *> &complete_threads,
		int node_id): queue(node_id, std::string("io-queue-") + itoa(node_id),
			IO_QUEUE_SIZE, INT_MAX, false), low_prio_queue(node_id,
				// TODO let's allow the low-priority queue to
				// be infinitely large for now.
				std::string("io-queue-low_prio-") + itoa(node_id),
				IO_QUEUE_SIZE, INT_MAX, false), partition(_partition)
{
	aio = new async_io(_partition, complete_threads, AIO_DEPTH_PER_FILE, node_id);
	this->node_id = node_id;
	num_accesses = 0;
	num_low_prio_accesses = 0;

	int ret = pthread_create(&id, NULL, process_requests, (void *) this);
	if (ret) {
		perror("pthread_create");
		exit(1);
	}
}

void ignore_flush(std::tr1::unordered_map<io_interface *, int> &ignored_flushes,
		const io_request &req)
{
	std::tr1::unordered_map<io_interface *, int>::iterator it
		= ignored_flushes.find(req.get_io());
	if (it == ignored_flushes.end()) {
		ignored_flushes.insert(std::pair<io_interface *, int>(req.get_io(), 0));
		it = ignored_flushes.find(req.get_io());
	}
	it->second++;
}

int disk_read_thread::process_low_prio_msg(message<io_request> &low_prio_msg,
		std::tr1::unordered_map<io_interface *, int> &ignored_flushes)
{
	int num_accesses = 0;
	int io_slots = aio->num_available_IO_slots();

	io_request req;
	while (low_prio_msg.has_next() && num_accesses < io_slots) {
		// We copy the request to the local stack.
		low_prio_msg.get_next(req);
		assert(req.get_num_bufs() == 1);
		// The request doesn't own the page, so the reference count
		// isn't increased while in the queue. Now we try to write
		// it back, we need to increase its reference. The only
		// safe way to do it is to use the search method of
		// the page cache.
		page_cache *cache = (page_cache *) req.get_priv();
		thread_safe_page *p = (thread_safe_page *) cache->search(
				req.get_offset());
		// The page has been evicted.
		if (p == NULL) {
			// The original page has been evicted, we should clear
			// the prepare-writeback flag on it.
			req.get_page(0)->set_prepare_writeback(false);
			num_ignored_low_prio_accesses++;
			ignore_flush(ignored_flushes, req);
			continue;
		}
		// If the original page has been evicted and the new page for
		// the offset has been added to the cache.
		if (p != req.get_page(0)) {
			p->dec_ref();
			// The original page has been evicted, we should clear
			// the prepare-writeback flag on it.
			req.get_page(0)->set_prepare_writeback(false);
			num_ignored_low_prio_accesses++;
			ignore_flush(ignored_flushes, req);
			continue;
		}
		// If we are here, it means the page is the one we are looking for.
		// We can be certain that the page won't be evicted because we have
		// a reference on it.

		// The object of page always exists, so we can always
		// lock a page.
		p->lock();
		// The page may have been written back by the applications.
		// But in either way, we need to reset the PREPARE_WRITEBACK
		// flag.
		p->set_prepare_writeback(false);
		// If the page is being written back or has been written back,
		// we can skip the request.
		if (p->is_io_pending() || !p->is_dirty()) {
			p->unlock();
			p->dec_ref();
			num_ignored_low_prio_accesses++;
			ignore_flush(ignored_flushes, req);
			continue;
		}
		assert(p == req.get_page(0));
		p->set_io_pending(true);
		p->unlock();
		num_accesses++;
		// The current private data points to the page cache.
		// Now the request owns the page, it's safe to point to
		// the page directly.
		req.set_priv(p);
		// This should block the thread.
		aio->access(&req, 1);
	}
	if (low_prio_msg.is_empty())
		low_prio_msg.clear();

	return num_accesses;
}

void disk_read_thread::run() {
	bind2node_id(node_id);
#ifdef DEBUG
	printf("disk read thread runs on node %d\n", node_id);
#endif
	aio->init();
	message<io_request> msg_buffer[LOCAL_BUF_SIZE];
	message<io_request> low_prio_msg;

	const int LOCAL_REQ_BUF_SIZE = IO_MSG_SIZE;
	io_request local_reqs[LOCAL_REQ_BUF_SIZE];
	std::tr1::unordered_map<io_interface *, int> ignored_flushes;
	while (true) {
		int num;
		num = queue.non_blocking_fetch(msg_buffer, LOCAL_BUF_SIZE);
		if (enable_debug)
			printf("I/O thread %d: queue size: %d, low-prio queue size: %d\n",
					get_node_id(), queue.get_num_entries(),
					low_prio_queue.get_num_entries());
		// The high-prio queue is empty.
		bool processed_low_prio = false;
		while (num == 0) {
			processed_low_prio = true;
			// we can process as many low-prio requests as possible,
			// but they should block the thread.
			if (!low_prio_queue.is_empty()
					&& aio->num_available_IO_slots() > 0) {
				if (low_prio_msg.is_empty()) {
					int num = low_prio_queue.fetch(&low_prio_msg, 1);
					assert(num == 1);
				}
				int ret = process_low_prio_msg(low_prio_msg, ignored_flushes);
				num_accesses += ret;
				num_low_prio_accesses += ret;
			}
			/* 
			 * this is the only thread that fetch requests from the queue.
			 * If there are no incoming requests and there are pending IOs,
			 * let's complete the pending IOs first.
			 */
			else if (aio->num_pending_ios() > 0) {
				aio->wait4complete(1);
			}
			// If there is no other work to do, let's wait for new requests.
			else
				break;

			// Let's try to fetch requests again.
			num = queue.non_blocking_fetch(msg_buffer, LOCAL_BUF_SIZE);
		}
		if (processed_low_prio) {
			// When we ignore flush requests, we also need to tell it.
			for (std::tr1::unordered_map<io_interface *, int>::iterator it
					= ignored_flushes.begin(); it != ignored_flushes.end(); it++) {
				io_interface *io = it->first;
				if (it->second > 0) {
					io->notify_completion(NULL, it->second);
					// We need to clear the counter for the next count.
					it->second = 0;
				}
			}
		}

		if (num == 0)
			num = queue.fetch(msg_buffer, LOCAL_BUF_SIZE, true, true);
		int num_flushes = flush_counter.get();
		if (num_flushes > 0) {
			// This thread is the only one that decreases the counter.
			flush_counter.dec(1);
			assert(flush_counter.get() >= 0);
			aio->flush_requests();
		}
		// We have been interrupted from waiting for IO requests.
		// Let's go back and try to process low-priority requests.
		if (num == 0)
			continue;
		for (int i = 0; i < num; i++) {
			int num_reqs = msg_buffer[i].get_num_objs();
			assert(num_reqs <= LOCAL_REQ_BUF_SIZE);
			msg_buffer[i].get_next_objs(local_reqs, num_reqs);
			aio->access(local_reqs, num_reqs);
			num_accesses += num_reqs;
			msg_buffer[i].clear();
		}
	}
	// TODO I need to call cleanup() of aio.
}

void *process_requests(void *arg)
{
	disk_read_thread *thread = (disk_read_thread *) arg;
	thread->run();
	return NULL;
}
