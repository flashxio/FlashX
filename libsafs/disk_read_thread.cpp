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

#include "cache.h"
#include "disk_read_thread.h"
#include "parameters.h"
#include "aio_private.h"
#include "debugger.h"

const int AIO_HIGH_PRIO_SLOTS = 7;
const int NUM_DIRTY_PAGES_TO_FETCH = 16 * 18;

disk_read_thread::disk_read_thread(const logical_file_partition &_partition,
		int node_id, page_cache *cache, int _disk_id): thread(
			std::string("io-thread-") + itoa(node_id), node_id), disk_id(
			_disk_id), queue(node_id, std::string("io-queue-") + itoa(node_id),
			IO_QUEUE_SIZE, INT_MAX, false), low_prio_queue(node_id,
				// TODO let's allow the low-priority queue to
				// be infinitely large for now.
				std::string("io-queue-low_prio-") + itoa(node_id),
				IO_QUEUE_SIZE, INT_MAX, false), partition(_partition), filter(
					open_files, _disk_id)
{
	this->cache = cache;
	aio = new async_io(_partition, AIO_DEPTH_PER_FILE, this);
#ifdef STATISTICS
	num_empty = 0;
	num_reads = 0;
	num_writes = 0;
	num_read_bytes = 0;
	num_write_bytes = 0;
	num_low_prio_accesses = 0;
	num_requested_flushes = 0;
	num_ignored_flushes_evicted = 0;
	num_ignored_flushes_cleaned = 0;
	num_ignored_flushes_old = 0;
	tot_flush_delay = 0;
	max_flush_delay = 0;
	min_flush_delay = LONG_MAX;
	num_msgs = 0;
#endif
	thread::start();
}

/**
 * Notify the IO issuer of the ignored flushes.
 * All flush requests must come from the same IO instance.
 */
void notify_ignored_flushes(io_request ignored_flushes[], int num_ignored)
{
	for (int i = 0; i < num_ignored; i++) {
		ignored_flushes[i].set_discarded(true);
		io_request *flush = &ignored_flushes[i];
		io_interface *io = flush->get_io();
		io->notify_completion(&flush, 1);
	}
}

int disk_read_thread::process_low_prio_msg(message<io_request> &low_prio_msg)
{
	int num_accesses = 0;

#ifdef STATISTICS
	struct timeval curr_time;
	gettimeofday(&curr_time, NULL);
#endif
	io_request req;
	stack_array<io_request> ignored_flushes(low_prio_msg.get_num_objs());
	int num_ignored = 0;
	while (low_prio_msg.has_next()
			&& aio->num_available_IO_slots() > AIO_HIGH_PRIO_SLOTS
			// We only submit requests to the disk when there aren't
			// high-prio requests.
			&& queue.is_empty()) {
		// We copy the request to the local stack.
		low_prio_msg.get_next(req);
#ifdef STATISTICS
		num_low_prio_accesses++;
#endif
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
#ifdef STATISTICS
			num_ignored_flushes_evicted++;
#endif
			ignored_flushes[num_ignored++] = req;
			continue;
		}
		// If the original page has been evicted and the new page for
		// the offset has been added to the cache.
		if (p != req.get_page(0)) {
			p->dec_ref();
			// The original page has been evicted, we should clear
			// the prepare-writeback flag on it.
			req.get_page(0)->set_prepare_writeback(false);
#ifdef STATISTICS
			num_ignored_flushes_evicted++;
#endif
			ignored_flushes[num_ignored++] = req;
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
		if (p->is_io_pending() || !p->is_dirty()
				|| p->get_flush_score() > DISCARD_FLUSH_THRESHOLD) {
			p->unlock();
			p->dec_ref();
#ifdef STATISTICS
			if (p->get_flush_score() > DISCARD_FLUSH_THRESHOLD)
				num_ignored_flushes_old++;
			else
				num_ignored_flushes_cleaned++;
#endif
			ignored_flushes[num_ignored++] = req;
			continue;
		}

#ifdef STATISTICS
		long delay = time_diff_us(req.get_timestamp(), curr_time);
		tot_flush_delay += delay;
		if (delay < min_flush_delay)
			min_flush_delay = delay;
		if (delay > max_flush_delay)
			max_flush_delay = delay;
		if (req.get_access_method() == READ) {
			num_reads++;
			num_read_bytes += req.get_size();
		}
		else {
			num_writes++;
			num_write_bytes += req.get_size();
		}
#endif

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

	if (num_ignored > 0)
		notify_ignored_flushes(ignored_flushes.data(), num_ignored);

	return num_accesses;
}

void disk_read_thread::run() {
	// First, check if we need to flush requests.
	int num_flushes = flush_counter.get();
	if (num_flushes > 0) {
		// This thread is the only one that decreases the counter.
		flush_counter.dec(1);
		assert(flush_counter.get() >= 0);
		aio->flush_requests();
	}

	message<io_request> msg_buffer[LOCAL_BUF_SIZE];
	message<io_request> low_prio_msg;

	const int LOCAL_REQ_BUF_SIZE = IO_MSG_SIZE;
	do {
		int num = queue.fetch(msg_buffer, LOCAL_BUF_SIZE);
		num_msgs += num;
		if (is_debug_enabled())
			printf("I/O thread %d: queue size: %d, low-prio queue size: %d\n",
					get_node_id(), queue.get_num_entries(),
					low_prio_queue.get_num_entries());
		// The high-prio queue is empty.
		while (num == 0) {
			// we can process as many low-prio requests as possible,
			// but they shouldn't block the thread.
			if (!low_prio_queue.is_empty()
					&& aio->num_available_IO_slots() > AIO_HIGH_PRIO_SLOTS) {
				if (low_prio_msg.is_empty()) {
					int num = low_prio_queue.fetch(&low_prio_msg, 1);
					num_msgs += num;
					assert(num == 1);
				}
				process_low_prio_msg(low_prio_msg);
			}
			/* 
			 * this is the only thread that fetch requests from the queue.
			 * If there are no incoming requests and there are pending IOs,
			 * let's complete the pending IOs first.
			 */
			else if (aio->num_pending_ios() > 0) {
				aio->wait4complete(1);
			}
			else if (cache) {
				int ret = cache->flush_dirty_pages(&filter, NUM_DIRTY_PAGES_TO_FETCH);
				if (ret == 0)
					break;
#ifdef STATISTICS
				num_requested_flushes += ret;
#endif
			}
			else
				break;

			// Let's try to fetch requests again.
			num = queue.fetch(msg_buffer, LOCAL_BUF_SIZE);
			num_msgs += num;
		}

		stack_array<io_request> local_reqs(LOCAL_REQ_BUF_SIZE);
		for (int i = 0; i < num; i++) {
			int num_reqs = msg_buffer[i].get_num_objs();
			assert(num_reqs <= LOCAL_REQ_BUF_SIZE);
			msg_buffer[i].get_next_objs(local_reqs.data(), num_reqs);
#ifdef STATISTICS
			for (int j = 0; j < num_reqs; j++) {
				if (local_reqs[j].get_access_method() == READ) {
					num_reads++;
					num_read_bytes += local_reqs[j].get_size();
				}
				else {
					num_writes++;
					num_write_bytes += local_reqs[j].get_size();
				}
			}
#endif
			aio->access(local_reqs.data(), num_reqs);
			msg_buffer[i].clear();
		}

		// We can't exit the loop if there are still pending AIO requests.
		// This thread is responsible for processing completed AIO requests.
	} while (aio->num_pending_ios() > 0);
}

int disk_read_thread::dirty_page_filter::filter(const thread_safe_page *pages[],
		int num, const thread_safe_page *returned_pages[])
{
	assert(mappers.size() == 1);
	int num_returned = 0;
	for (int i = 0; i < num; i++) {
		int id = mappers[0]->map2file(pages[i]->get_offset());
		if (this->disk_id == id)
			returned_pages[num_returned++] = pages[i];
	}
	return num_returned;
}

void disk_read_thread::print_state()
{
	printf("io thread %d has %d reqs and %d low-prio reqs in the queue\n",
			disk_id, queue.get_num_objs(), low_prio_queue.get_num_objs());
	aio->print_state();
}
