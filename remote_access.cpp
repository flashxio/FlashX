#include "remote_access.h"
#include "parameters.h"

const int INIT_DISK_QUEUE_SIZE = 32;
const int MAX_DISK_CACHED_REQS = 32;

class request_sender
{
	fifo_queue<io_request> buf;
	blocking_FIFO_queue<io_request> *queue;
public:
	/**
	 * buf_size: the number of messages that can be buffered in the sender.
	 */
	request_sender(blocking_FIFO_queue<io_request> *queue): buf(
			INIT_DISK_QUEUE_SIZE, true) {
		this->queue = queue;
	}

	int flush(bool blocking) {
		if (buf.is_empty()) {
			return 0;
		}
		if (blocking)
			return queue->add(&buf);
		else
			return queue->non_blocking_add(&buf);
	}

	void flush_all() {
		while (!buf.is_empty())
			queue->add(&buf);
	}

	int get_num_remaining() {
		return buf.get_num_entries();
	}

	int send_cached(io_request *msg) {
		if (buf.is_full())
			buf.expand_queue(buf.get_size() * 2);
		buf.push_back(*msg);
		return 1;
	}

	blocking_FIFO_queue<io_request> *get_queue() const {
		return queue;
	}
};

remote_disk_access::remote_disk_access(disk_read_thread **remotes,
		int num_remotes, file_mapper *mapper,
		int node_id): io_interface(node_id)
{
	assert(num_remotes == mapper->get_num_files());
	senders = new request_sender *[num_remotes];
	num_senders = num_remotes;
	total_size = 0;
	local_size = 0;
	// create a msg sender for each disk read thread.
	for (int i = 0; i < num_remotes; i++) {
		total_size += remotes[i]->get_size();
		if (remotes[i]->get_node_id() == node_id)
			local_size += remotes[i]->get_size();
		blocking_FIFO_queue<io_request> *queue = remotes[i]->get_queue();
		senders[i] = new request_sender(queue);
	}
	cb = NULL;
	this->block_mapper = mapper;
}

remote_disk_access::~remote_disk_access()
{
	for (int i = 0; i < num_senders; i++)
		delete senders[i];
	delete [] senders;
}

io_interface *remote_disk_access::clone() const
{
	remote_disk_access *copy = new remote_disk_access(this->get_node_id());
	copy->num_senders = this->num_senders;
	copy->senders = new request_sender *[this->num_senders];
	for (int i = 0; i < copy->num_senders; i++) {
		copy->senders[i] = new request_sender(this->senders[i]->get_queue());
	}
	copy->cb = this->cb;
	copy->total_size = this->total_size;
	copy->local_size = this->local_size;
	copy->block_mapper = this->block_mapper;
	return copy;
}

void remote_disk_access::cleanup()
{
	for (int i = 0; i < num_senders; i++) {
		senders[i]->flush_all();
	}
	int num;
	do {
		num = 0;
		for (int i = 0; i < num_senders; i++) {
			num += senders[i]->get_queue()->get_num_entries();
		}
		/* 
		 * if there are still messages in the queue, wait.
		 * this might be the best I can do right now
		 * unless the queues can notify me when they are
		 * empty.
		 */
		if (num > 0)
			usleep(100000);
	} while (num > 0);
}

ssize_t remote_disk_access::access(io_request *requests, int num)
{
	int num_remaining = 0;
	for (int i = 0; i < num; i++) {
		assert(requests[i].get_size() > 0);
		// TODO data is striped on disks, we have to make sure the data
		// can be accessed from the remote disk access.
		// and I assume the request size is aligned with the strip size.
		off_t pg_off = requests[i].get_offset() / PAGE_SIZE;
		int idx = block_mapper->map2file(pg_off);
		// The cache inside a sender is extensible, so it can absorb
		// all requests.
		int ret = senders[idx]->send_cached(&requests[i]);
		num_remaining += senders[idx]->get_num_remaining();
		assert(ret == 1);
	}

	if (num_remaining > MAX_DISK_CACHED_REQS) {
		num_remaining = 0;
		// Now let's flush requests to the queues, but we first try to
		// flush requests non-blockingly.
		for (int i = 0; i < num_senders; i++) {
			senders[i]->flush(false);
			num_remaining += senders[i]->get_num_remaining();
		}

		int base_idx;
		if (num_senders == 1)
			base_idx = 0;
		else
			base_idx = random() % num_senders;
		int i = 0;
		// We only allow cache that many requests. If we have more than
		// we want, continue flushing, but try harder this time.
		while (num_remaining > MAX_DISK_CACHED_REQS) {
			int idx = (base_idx + i) % num_senders;
			int orig_remaining = senders[idx]->get_num_remaining();
			senders[idx]->flush(true);
			assert(senders[idx]->get_num_remaining() == 0);
			num_remaining -= orig_remaining;
			i++;
		}
	}
	return 0;
}
