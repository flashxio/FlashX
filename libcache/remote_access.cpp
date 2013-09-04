#include "remote_access.h"
#include "parameters.h"
#include "slab_allocator.h"
#include "disk_read_thread.h"
#include "file_mapper.h"

void remote_disk_access::notify_completion(io_request *reqs[], int num)
{
	// There are a few cases for the incoming requests.
	//	the requests issued by the upper layer IO;
	//	the requests split by the current IO;
	//	the requests issued by an application.
	io_request *from_upper[num];
	io_request *from_app[num];
	io_interface *upper_io = NULL;
	int num_from_upper = 0;
	int num_from_app = 0;
	std::vector<io_request *> completes;
	for (int i = 0; i < num; i++) {
		// The requests issued by the upper layer IO.
		if (reqs[i]->get_io() != this) {
			if (upper_io == NULL)
				upper_io = reqs[i]->get_io();
			else
				// They should be from the same upper layer IO.
				assert(upper_io == reqs[i]->get_io());
			from_upper[num_from_upper++] = reqs[i];
			continue;
		}
		if (reqs[i]->get_io() == this && !reqs[i]->is_extended_req()) {
			from_app[num_from_app++] = reqs[i];
			continue;
		}

		io_request *orig = reqs[i]->get_orig();
		io_request *req = reqs[i];
		orig->inc_complete_count();
		if (orig->complete_size(req->get_size()))
			completes.push_back(orig);
		else
			orig->dec_complete_count();
	}
	if (num_from_upper > 0)
		upper_io->notify_completion(from_upper, num_from_upper);
	if (num_from_app > 0 && this->get_callback())
		this->get_callback()->invoke(from_app, num_from_app);

	for (unsigned i = 0; i < completes.size(); i++) {
		io_request *orig = completes[i];
		assert(orig->is_extended_req());
		io_interface *io = orig->get_io();
		// It's from an application.
		if (io == this) {
			if (io->get_callback())
				io->get_callback()->invoke(&orig, 1);
		}
		else
			io->notify_completion(&orig, 1);
		orig->dec_complete_count();
		orig->wait4unref();
		// Now we can delete it.
		delete orig;
	}
}

remote_disk_access::remote_disk_access(const std::vector<disk_read_thread *> &remotes,
		aio_complete_thread *complete_thread, file_mapper *mapper, int node_id,
		int max_reqs): io_interface(node_id), max_disk_cached_reqs(max_reqs)
{
	this->io_threads = remotes;
	if (complete_thread == NULL)
		this->complete_queue = NULL;
	else
		this->complete_queue = complete_thread->get_queue();
	// TODO I need to deallocate it later.
	msg_allocator = new slab_allocator(IO_MSG_SIZE, IO_MSG_SIZE * 1024,
			INT_MAX, node_id);
	senders.resize(remotes.size());
	low_prio_senders.resize(remotes.size());
	// create a msg sender for each disk read thread.
	for (unsigned i = 0; i < remotes.size(); i++) {
		senders[i] = request_sender::create(node_id, msg_allocator,
				remotes[i]->get_queue());
		low_prio_senders[i] = request_sender::create(node_id, msg_allocator,
				remotes[i]->get_low_prio_queue());
	}
	cb = NULL;
	this->block_mapper = mapper;
	num_completed_reqs = 0;
}

remote_disk_access::~remote_disk_access()
{
	assert(senders.size() == low_prio_senders.size());
	int num_senders = senders.size();
	for (int i = 0; i < num_senders; i++) {
		request_sender::destroy(senders[i]);
		request_sender::destroy(low_prio_senders[i]);
	}
}

io_interface *remote_disk_access::clone() const
{
	remote_disk_access *copy = new remote_disk_access(this->get_node_id(),
			this->max_disk_cached_reqs);
	copy->io_threads = this->io_threads;
	copy->senders.resize(this->senders.size());
	copy->low_prio_senders.resize(this->low_prio_senders.size());
	assert(copy->senders.size() == copy->low_prio_senders.size());
	for (unsigned i = 0; i < copy->senders.size(); i++) {
		copy->senders[i] = request_sender::create(this->get_node_id(),
				msg_allocator, this->senders[i]->get_queue());
		copy->low_prio_senders[i] = request_sender::create(this->get_node_id(),
				msg_allocator, this->low_prio_senders[i]->get_queue());
	}
	copy->cb = this->cb;
	copy->block_mapper = this->block_mapper;
	copy->complete_queue = this->complete_queue;
	copy->msg_allocator = this->msg_allocator;
	return copy;
}

void remote_disk_access::cleanup()
{
	for (unsigned i = 0; i < senders.size(); i++) {
		senders[i]->flush_all();
		low_prio_senders[i]->flush_all();
	}
	int num;
	do {
		num = 0;
		assert(senders.size() == low_prio_senders.size());
		for (unsigned i = 0; i < senders.size(); i++) {
			num += senders[i]->get_queue()->get_num_entries();
			num += low_prio_senders[i]->get_queue()->get_num_entries();
		}
		/* 
		 * if there are still messages in the queue, wait.
		 * this might be the best I can do right now
		 * unless the queues can notify me when they are
		 * empty.
		 */
		if (num > 0) {
			// Let's wake up all IO threads if there are still
			// some low-priority requests that need to be processed.
			for (unsigned i = 0; i < senders.size(); i++) {
				senders[i]->get_queue()->wakeup();
			}
			usleep(100000);
		}
	} while (num > 0);

	for (unsigned i = 0; i < io_threads.size(); i++)
		io_threads[i]->flush_requests();
}

void remote_disk_access::access(io_request *requests, int num,
		io_status *status)
{
	// It marks whether a low-priority sender gets a request.
	bool has_msgs[low_prio_senders.size()];
	memset(has_msgs, 0, sizeof(has_msgs[0]) * low_prio_senders.size());

	bool syncd = false;
	for (int i = 0; i < num; i++) {
		assert(requests[i].get_size() >= MIN_BLOCK_SIZE);

		if (requests[i].is_flush()) {
			syncd = true;
			continue;
		}
		else if (requests[i].is_sync()) {
			syncd = true;
		}

		// If the request accesses one RAID block, it's simple.
		if (inside_RAID_block(requests[i])) {
			off_t pg_off = requests[i].get_offset() / PAGE_SIZE;
			int idx = block_mapper->map2file(pg_off);
			// The cache inside a sender is extensible, so it can absorb
			// all requests.
			int ret;
			if (requests[i].is_high_prio())
				ret = senders[idx]->send_cached(&requests[i]);
			else {
				has_msgs[idx] = true;
				ret = low_prio_senders[idx]->send_cached(&requests[i]);
			}
			assert(ret == 1);
		}
		else {
			// If the request accesses multiple RAID blocks, we have to
			// split the request.
			// I still use the default memory allocator, but since it is used
			// when the request size is large, it should normally be OK.
			// TODO I can use slab allocators later.
			io_request *orig = new io_request(true);
			// global_cached_io doesn't issue requests across a block boundary.
			// It can only be application issued requst, so it shouldn't have
			// extension.
			assert(!requests[i].is_extended_req());
			orig->init(requests[i]);
			off_t end = orig->get_offset() + orig->get_size();
			const off_t RAID_block_size = params.get_RAID_block_size() * PAGE_SIZE;
			for (off_t begin = orig->get_offset(); begin < end;
					begin = ROUND(begin + RAID_block_size, RAID_block_size)) {
				io_request req(true);
				int size = ROUND(begin + RAID_block_size, RAID_block_size) - begin;
				size = min(size, end - begin);
				// It only supports to extract a specified request from
				// a single-buffer request.
				extract_pages(*orig, begin, size / PAGE_SIZE, req);
				req.set_orig(orig);
				req.set_io(this);
				assert(inside_RAID_block(req));

				// Send a request.
				off_t pg_off = req.get_offset() / PAGE_SIZE;
				int idx = block_mapper->map2file(pg_off);
				// The cache inside a sender is extensible, so it can absorb
				// all requests.
				int ret;
				if (req.is_high_prio())
					ret = senders[idx]->send_cached(&req);
				else {
					has_msgs[idx] = true;
					ret = low_prio_senders[idx]->send_cached(&req);
				}
				assert(ret == 1);
			}
		}
	}

	int num_remaining = 0;
	assert(senders.size() == low_prio_senders.size());
	for (unsigned i = 0; i < senders.size(); i++) {
		num_remaining += senders[i]->get_num_remaining();
		num_remaining += low_prio_senders[i]->get_num_remaining();
	}
	if (num_remaining > MAX_DISK_CACHED_REQS) {
		flush_requests(MAX_DISK_CACHED_REQS);
	}
	if (syncd)
		flush_requests();

	if (status)
		for (int i = 0; i < num; i++)
			status[i] = IO_PENDING;

	// The IO threads are never blocked by the low-priority queues,
	// but they may be blocked by the regular message queues.
	// If so, we need to explicitly wake up the IO threads.
	for (unsigned i = 0; i < senders.size(); i++) {
		if (has_msgs[i])
			senders[i]->get_queue()->wakeup();
	}
}

void remote_disk_access::flush_requests()
{
	flush_requests(0);
}

void remote_disk_access::flush_requests(int max_cached)
{
	if (complete_queue)
		num_completed_reqs += complete_queue->process(1000, false);
	int num_high_prio_remaining = 0;
	int num_low_prio_remaining = 0;
	// Now let's flush requests to the queues, but we first try to
	// flush requests non-blockingly.
	assert(senders.size() == low_prio_senders.size());
	int num_senders = senders.size();
	for (int i = 0; i < num_senders; i++) {
		senders[i]->flush(false);
		low_prio_senders[i]->flush(false);
		num_high_prio_remaining += senders[i]->get_num_remaining();
		num_low_prio_remaining += low_prio_senders[i]->get_num_remaining();
	}
	// If all requests have been flushed successfully, return immediately.
	if (num_high_prio_remaining + num_low_prio_remaining == 0)
		return;

	int base_idx;
	if (num_senders == 1)
		base_idx = 0;
	else
		base_idx = random() % num_senders;
	int i = 0;
	// We only allow cache that many requests. If we have more than
	// we want, continue flushing, but try harder this time.
	while (num_high_prio_remaining + num_low_prio_remaining > max_cached
			&& num_high_prio_remaining > 0) {
		int idx = (base_idx + i) % num_senders;
		int orig_remaining = senders[idx]->get_num_remaining();
		senders[idx]->flush(true);
		int num_flushed = orig_remaining - senders[idx]->get_num_remaining();
		num_high_prio_remaining -= num_flushed;
		i++;
	}
	// When we reach this point, it means either that the total number of
	// remaining requests is lower than max_cached or there aren't high-
	// priority requests left.
	// In this way, we can make sure high-priority requests have been moved
	// to the IO threads before low-priority requests are moved.
	i = 0;
	while (num_low_prio_remaining > max_cached) {
		int idx = (base_idx + i) % num_senders;
		int orig_remaining = low_prio_senders[idx]->get_num_remaining();
		low_prio_senders[idx]->flush(true);
		int num_flushed = orig_remaining
			- low_prio_senders[idx]->get_num_remaining();
		num_low_prio_remaining -= num_flushed;
		i++;
	}
}

int remote_disk_access::get_file_id() const
{
	return block_mapper->get_file_id();
}
