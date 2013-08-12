#include "remote_access.h"
#include "parameters.h"

const int INIT_DISK_QUEUE_SIZE = 32;

class request_assemble_callback: public callback
{
public:
	int invoke(io_request *requests[], int num);
};

int request_assemble_callback::invoke(io_request *requests[], int num)
{
	std::vector<io_request *> completes;
	for (int i = 0; i < num; i++) {
		io_request *orig = requests[i]->get_orig();
		io_request *req = requests[i];
		orig->inc_complete_count();
		if (orig->complete_size(req->get_size()))
			completes.push_back(orig);
		else
			orig->dec_complete_count();
	}

	for (unsigned i = 0; i < completes.size(); i++) {
		io_request *orig = completes[i];
		if (orig->get_io()->get_callback())
			orig->get_io()->get_callback()->invoke(&orig, 1);
		orig->dec_complete_count();
		orig->wait4unref();
		// Now we can delete it.
		delete orig;
	}
	return 0;
}

class request_intercepter: public io_interface
{
	callback *cb;
public:
	request_intercepter(): io_interface(-1) {
		cb = new request_assemble_callback();
	}

	virtual callback *get_callback() {
		return cb;
	}
};

remote_disk_access::remote_disk_access(disk_read_thread **remotes,
		aio_complete_thread *complete_thread, int num_remotes,
		file_mapper *mapper, int node_id): io_interface(node_id)
{
	if (complete_thread == NULL)
		this->complete_queue = NULL;
	else
		this->complete_queue = complete_thread->get_queue();
	assert(num_remotes == mapper->get_num_files());
	senders = new request_sender *[num_remotes];
	low_prio_senders = new request_sender *[num_remotes];
	num_senders = num_remotes;
	// create a msg sender for each disk read thread.
	for (int i = 0; i < num_remotes; i++) {
		blocking_FIFO_queue<io_request> *queue = remotes[i]->get_queue();
		senders[i] = new request_sender(queue, INIT_DISK_QUEUE_SIZE);
		low_prio_senders[i] = new request_sender(remotes[i]->get_low_prio_queue(),
				INIT_DISK_QUEUE_SIZE);
	}
	cb = NULL;
	this->block_mapper = mapper;
	num_completed_reqs = 0;
	this->req_intercepter = new request_intercepter();
}

remote_disk_access::~remote_disk_access()
{
	for (int i = 0; i < num_senders; i++) {
		delete senders[i];
		delete low_prio_senders[i];
	}
	delete [] senders;
	delete [] low_prio_senders;
	delete req_intercepter;
}

io_interface *remote_disk_access::clone() const
{
	remote_disk_access *copy = new remote_disk_access(this->get_node_id());
	copy->num_senders = this->num_senders;
	copy->senders = new request_sender *[this->num_senders];
	copy->low_prio_senders = new request_sender *[this->num_senders];
	for (int i = 0; i < copy->num_senders; i++) {
		copy->senders[i] = new request_sender(this->senders[i]->get_queue(),
				INIT_DISK_QUEUE_SIZE);
		copy->low_prio_senders[i] = new request_sender(
				this->low_prio_senders[i]->get_queue(), INIT_DISK_QUEUE_SIZE);
	}
	copy->cb = this->cb;
	copy->block_mapper = this->block_mapper;
	copy->complete_queue = this->complete_queue;
	return copy;
}

void remote_disk_access::cleanup()
{
	for (int i = 0; i < num_senders; i++) {
		senders[i]->flush_all();
		low_prio_senders[i]->flush_all();
	}
	int num;
	do {
		num = 0;
		for (int i = 0; i < num_senders; i++) {
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
			for (int i = 0; i < num_senders; i++) {
				senders[i]->get_queue()->wakeup();
			}
			usleep(100000);
		}
	} while (num > 0);
}

void remote_disk_access::access(io_request *requests, int num,
		io_status *status)
{
	// It marks whether a low-priority sender gets a request.
	bool has_msgs[num_senders];
	memset(has_msgs, 0, sizeof(has_msgs[0]) * num_senders);

	bool syncd = false;
	for (int i = 0; i < num; i++) {
		requests[i].set_file_id(block_mapper->get_file_id());
		assert(requests[i].get_size() > 0);

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
			io_request *orig = new io_request();
			*orig = requests[i];
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
				req.set_file_id(orig->get_file_id());
				req.set_orig(orig);
				req.set_io(req_intercepter);
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
	for (int i = 0; i < num_senders; i++) {
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
	for (int i = 0; i < num_senders; i++) {
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
