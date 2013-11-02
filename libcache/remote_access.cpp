#include "remote_access.h"
#include "parameters.h"
#include "slab_allocator.h"
#include "disk_read_thread.h"
#include "file_mapper.h"

const int NUM_PROCESS_COMPLETED_REQS = 8;

/**
 * The maximal number of pending IOs is decided by the maximal pending IOs
 * allowed on each I/O thread divided by the number of remote_disk_access.
 * Thus, the performance isn't decided by the number of remote_disk_access.
 */
int remote_disk_access::get_max_num_pending_ios() const
{
	return io_interface::get_max_num_pending_ios() *
		io_threads.size() / num_ios.get();
}

void remote_disk_access::notify_completion(io_request *reqs[], int num)
{
	stack_array<io_request> req_copies(num);
	for (int i = 0; i < num; i++) {
		req_copies[i] = *reqs[i];
		assert(req_copies[i].get_io());
	}

	int ret = complete_queue.add(req_copies.data(), num);
	assert(ret == num);
	get_thread()->activate();
}

remote_disk_access::remote_disk_access(const std::vector<disk_read_thread *> &remotes,
		file_mapper *mapper, thread *t, int max_reqs): io_interface(
			// TODO I hope the queue size is large enough.
			t), max_disk_cached_reqs(max_reqs), complete_queue(t->get_node_id(),
				COMPLETE_QUEUE_SIZE)
{
	int node_id = t->get_node_id();
	num_ios.inc(1);
	this->io_threads = remotes;
	// TODO I need to deallocate it later.
	msg_allocator = new slab_allocator(std::string("disk_msg_allocator-")
			+ itoa(node_id), IO_MSG_SIZE * sizeof(io_request),
			IO_MSG_SIZE * sizeof(io_request) * 1024, INT_MAX, node_id);
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

io_interface *remote_disk_access::clone(thread *t) const
{
	// An IO may not be associated to any threads.
	ASSERT_TRUE(t);
	num_ios.inc(1);
	remote_disk_access *copy = new remote_disk_access(io_threads,
			block_mapper, t, this->max_disk_cached_reqs);
	copy->cb = this->cb;
	return copy;
}

void remote_disk_access::cleanup()
{
	process_completed_requests(complete_queue.get_num_entries());

	num_ios.dec(1);
	for (unsigned i = 0; i < senders.size(); i++) {
		senders[i]->flush();
		low_prio_senders[i]->flush();
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
			for (unsigned i = 0; i < io_threads.size(); i++) {
				io_threads[i]->activate();
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
	ASSERT_EQ(get_thread(), thread::get_curr_thread());
	num_issued_reqs.inc(num);

	bool syncd = false;
	for (int i = 0; i < num; i++) {
		ASSERT_LTEQ(requests[i].get_size(), MIN_BLOCK_SIZE);

		if (requests[i].is_flush()) {
			syncd = true;
			num_completed_reqs.inc(1);
			continue;
		}
		else if (requests[i].is_sync()) {
			syncd = true;
		}

		// If the request accesses one RAID block, it's simple.
		if (requests[i].inside_RAID_block()) {
			off_t pg_off = requests[i].get_offset() / PAGE_SIZE;
			int idx = block_mapper->map2file(pg_off);
			// The cache inside a sender is extensible, so it can absorb
			// all requests.
			int ret;
			if (requests[i].is_high_prio())
				ret = senders[idx]->send_cached(&requests[i]);
			else {
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
			io_req_extension *ext = new io_req_extension();
			io_request *orig = new io_request(ext, 0, 0, NULL, 0);
			// global_cached_io doesn't issue requests across a block boundary.
			// It can only be application issued requst, so it shouldn't have
			// extension.
			assert(!requests[i].is_extended_req());
			orig->init(requests[i]);
			off_t end = orig->get_offset() + orig->get_size();
			const off_t RAID_block_size = params.get_RAID_block_size() * PAGE_SIZE;
			for (off_t begin = orig->get_offset(); begin < end;
					begin = ROUND(begin + RAID_block_size, RAID_block_size)) {
				io_req_extension *ext = new io_req_extension();
				ext->set_orig(orig);
				io_request req(ext, 0, 0, NULL, 0);
				int size = ROUND(begin + RAID_block_size, RAID_block_size) - begin;
				size = min(size, end - begin);
				// It only supports to extract a specified request from
				// a single-buffer request.
				orig->extract(begin, size, req);
				req.set_io(this);
				assert(req.inside_RAID_block());

				// Send a request.
				off_t pg_off = req.get_offset() / PAGE_SIZE;
				int idx = block_mapper->map2file(pg_off);
				// The cache inside a sender is extensible, so it can absorb
				// all requests.
				int ret;
				if (req.is_high_prio())
					ret = senders[idx]->send_cached(&req);
				else {
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
}

void remote_disk_access::flush_requests()
{
	flush_requests(0);
}

int remote_disk_access::process_completed_requests(int num)
{
	if (num > 0) {
		stack_array<io_request> reqs(num);
		int ret = complete_queue.fetch(reqs.data(), num);
		process_completed_requests(reqs.data(), ret);
		return ret;
	}
	else
		return 0;
}

int remote_disk_access::process_completed_requests(io_request reqs[], int num)
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
	int num_part_reqs = 0;
	std::vector<io_request *> completes;
	for (int i = 0; i < num; i++) {
		assert(reqs[i].get_io());
		// The requests issued by the upper layer IO.
		if (reqs[i].get_io() != this) {
			if (upper_io == NULL)
				upper_io = reqs[i].get_io();
			else
				// They should be from the same upper layer IO.
				assert(upper_io == reqs[i].get_io());
			from_upper[num_from_upper++] = &reqs[i];
			continue;
		}
		if (reqs[i].get_io() == this && !reqs[i].is_extended_req()) {
			from_app[num_from_app++] = &reqs[i];
			continue;
		}

		io_request *orig = reqs[i].get_orig();
		io_request *req = &reqs[i];
		orig->inc_complete_count();
		if (orig->complete_size(req->get_size()))
			completes.push_back(orig);
		else {
			orig->dec_complete_count();
			num_part_reqs++;
		}
		delete req->get_extension();
	}
	if (num_from_upper > 0) {
		assert(upper_io);
		upper_io->notify_completion(from_upper, num_from_upper);
	}
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
		delete orig->get_extension();
		delete orig;
	}

	num_completed_reqs.inc(num - num_part_reqs);
	return num - num_part_reqs;
}

void remote_disk_access::flush_requests(int max_cached)
{
	// Now let's flush requests to the queues, but we first try to
	// flush requests non-blockingly.
	assert(senders.size() == low_prio_senders.size());
	int num_senders = senders.size();
	for (int i = 0; i < num_senders; i++) {
		senders[i]->flush();
		low_prio_senders[i]->flush();
		assert(senders[i]->get_num_remaining() == 0);
		assert(low_prio_senders[i]->get_num_remaining() == 0);
	}
	for (unsigned i = 0; i < io_threads.size(); i++) {
		io_threads[i]->activate();
	}
}

int remote_disk_access::get_file_id() const
{
	return block_mapper->get_file_id();
}

/**
 * We wait for at least the specified number of requests to complete.
 */
int remote_disk_access::wait4complete(int num_to_complete)
{
	flush_requests();
	int pending = num_pending_ios();
	num_to_complete = min(pending, num_to_complete);

	process_completed_requests(complete_queue.get_num_entries());
	int iters = 0;
	while (pending - num_pending_ios() < num_to_complete) {
		iters++;
		get_thread()->wait();
		process_completed_requests(complete_queue.get_num_entries());
	}
	return pending - num_pending_ios();
}

void remote_disk_access::print_state()
{
	printf("remote_io %d has %d pending reqs, %d completed reqs\n",
			get_io_idx(), num_pending_ios(), complete_queue.get_num_entries());
	for (unsigned i = 0; i < senders.size(); i++)
		printf("\tsender %d: remain %d reqs\n", i,
				senders[i]->get_num_remaining());
	for (unsigned i = 0; i < low_prio_senders.size(); i++)
		printf("\tlow-prio sender %d: remain %d reqs\n", i,
				low_prio_senders[i]->get_num_remaining());
}

atomic_integer remote_disk_access::num_ios;
