#include "remote_access.h"
#include "parameters.h"

const int INIT_DISK_QUEUE_SIZE = 32;

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
	total_size = 0;
	local_size = 0;
	// create a msg sender for each disk read thread.
	for (int i = 0; i < num_remotes; i++) {
		total_size += remotes[i]->get_size();
		if (remotes[i]->get_node_id() == node_id)
			local_size += remotes[i]->get_size();
		blocking_FIFO_queue<io_request> *queue = remotes[i]->get_queue();
		senders[i] = new request_sender(queue, INIT_DISK_QUEUE_SIZE);
		low_prio_senders[i] = new request_sender(remotes[i]->get_low_prio_queue(),
				INIT_DISK_QUEUE_SIZE);
	}
	cb = NULL;
	this->block_mapper = mapper;
	num_completed_reqs = 0;
}

remote_disk_access::~remote_disk_access()
{
	for (int i = 0; i < num_senders; i++) {
		delete senders[i];
		delete low_prio_senders[i];
	}
	delete [] senders;
	delete [] low_prio_senders;
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
	copy->total_size = this->total_size;
	copy->local_size = this->local_size;
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

	for (int i = 0; i < num; i++) {
		assert(requests[i].get_size() > 0);
		// TODO data is striped on disks, we have to make sure the data
		// can be accessed from the remote disk access.
		// and I assume the request size is aligned with the strip size.
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

	int num_remaining = 0;
	for (int i = 0; i < num_senders; i++) {
		num_remaining += senders[i]->get_num_remaining();
		num_remaining += low_prio_senders[i]->get_num_remaining();
	}
	if (num_remaining > MAX_DISK_CACHED_REQS) {
		flush_requests(MAX_DISK_CACHED_REQS);
	}
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
