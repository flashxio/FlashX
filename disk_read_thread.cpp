#include "disk_read_thread.h"
#include "parameters.h"
#include "container.cpp"

/* just call the callback of the initiator. */
class initiator_callback: public callback
{
public:
	int invoke(io_request *rqs[], int num) {
		// TODO let's do it for now.
		for (int i = 0; i < num; i++) {
			io_request *rq = rqs[i];
			io_interface *io = rq->get_io();
			if (io->get_callback())
				io->get_callback()->invoke(&rq, 1);
		}
		return 0;
	}
};

disk_read_thread::disk_read_thread(const logical_file_partition &partition,
		const std::tr1::unordered_map<int, aio_complete_thread *> &complete_threads,
		long size, int node_id): queue(partition.get_file_name(0),
			IO_QUEUE_SIZE, IO_QUEUE_SIZE), low_prio_queue(
				// TODO let's allow the low-priority queue to
				// be infinitely large for now.
				partition.get_file_name(0) + "-low_prio", IO_QUEUE_SIZE, INT_MAX) {
	aio = new async_io(partition, complete_threads, size, AIO_DEPTH_PER_FILE, node_id);
	aio->set_callback(new initiator_callback());
	this->node_id = node_id;
	num_accesses = 0;
	num_low_prio_accesses = 0;

	int ret = pthread_create(&id, NULL, process_requests, (void *) this);
	if (ret) {
		perror("pthread_create");
		exit(1);
	}
}

void disk_read_thread::run() {
	bind2node_id(node_id);
	printf("disk read thread runs on node %d\n", node_id);
	aio->init();
	io_request reqs[AIO_DEPTH_PER_FILE];
	while (true) {
		int num;
		num = queue.non_blocking_fetch(reqs, AIO_DEPTH_PER_FILE);
		// The high-prio queue is empty.
		while (num == 0) {
			// we can process as many low-prio requests as possible,
			// but they should block the thread.
			if (!low_prio_queue.is_empty()
					&& aio->num_available_IO_slots() > 0) {
				int io_slots = aio->num_available_IO_slots();
				assert(io_slots <= AIO_DEPTH_PER_FILE);
				int num = low_prio_queue.fetch(reqs, io_slots);
				for (int i = 0; i < num; i++) {
					// The request doesn't own the page, so the reference count
					// isn't increased while in the queue. Now we try to write
					// it back, we need to increase its reference. The only
					// safe way to do it is to use the search method of
					// the page cache.
					page_cache *cache = (page_cache *) reqs[i].get_priv();
					thread_safe_page *p = (thread_safe_page *) cache->search(
							reqs[i].get_offset());
					if (p == NULL)
						continue;
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
						continue;
					}
					p->set_io_pending(true);
					p->unlock();
					num_accesses++;
					// The current private data points to the page cache.
					// Now the request owns the page, it's safe to point to
					// the page directly.
					reqs[i].set_priv(p);
					// This should block the thread.
					aio->access(&reqs[i], 1);
					num_low_prio_accesses++;
				}
			}
			/* 
			 * this is the only thread that fetch requests from the queue.
			 * If there are no incoming requests and there are pending IOs,
			 * let's complete the pending IOs first.
			 */
			else if (aio->num_pending_IOs() > 0) {
				aio->wait4complete();
			}
			// If there is no other work to do, let's wait for new requests.
			else
				break;

			// Let's try to fetch requests again.
			num = queue.non_blocking_fetch(reqs, AIO_DEPTH_PER_FILE);
		}

		if (num == 0)
			num = queue.fetch(reqs, AIO_DEPTH_PER_FILE, true, true);
		// We have been interrupted from waiting for IO requests.
		// Let's go back and try to process low-priority requests.
		if (num == 0)
			continue;
		num_accesses += num;
		aio->access(reqs, num);
	}
	// TODO I need to call cleanup() of aio.
}

void *process_requests(void *arg)
{
	disk_read_thread *thread = (disk_read_thread *) arg;
	printf("disk_read_thread: pid: %d, tid: %ld\n", getpid(), gettid());
	thread->run();
	return NULL;
}
