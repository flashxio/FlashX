#include "remote_access.h"
#include "parameters.h"

remote_disk_access::remote_disk_access(disk_read_thread **remotes,
			int num_remotes, int node_id): io_interface(node_id)
{
	senders = new msg_sender<io_request> *[num_remotes];
	queues = new thread_safe_FIFO_queue<io_request> *[num_remotes];
	num_senders = num_remotes;
	total_size = 0;
	local_size = 0;
	// create a msg sender for each disk read thread.
	for (int i = 0; i < num_remotes; i++) {
		total_size += remotes[i]->get_size();
		if (remotes[i]->get_node_id() == node_id)
			local_size += remotes[i]->get_size();
		thread_safe_FIFO_queue<io_request> *queue = remotes[i]->get_queue();
		senders[i] = new msg_sender<io_request>(MSG_SEND_BUF_SIZE,
				&queue, 1);
		queues[i] = queue;
	}
	cb = NULL;
}

remote_disk_access::~remote_disk_access()
{
	for (int i = 0; i < num_senders; i++)
		delete senders[i];
	delete [] senders;
	delete [] queues;
}

io_interface *remote_disk_access::clone() const
{
	remote_disk_access *copy = new remote_disk_access(this->get_node_id());
	copy->num_senders = this->num_senders;
	copy->senders = new msg_sender<io_request> *[this->num_senders];
	copy->queues = new thread_safe_FIFO_queue<io_request> *[this->num_senders];
	for (int i = 0; i < copy->num_senders; i++) {
		copy->queues[i] = this->queues[i];
		copy->senders[i] = new msg_sender<io_request>(MSG_SEND_BUF_SIZE,
				&copy->queues[i], 1);
	}
	copy->cb = this->cb;
	copy->total_size = this->total_size;
	copy->local_size = this->local_size;
	return copy;
}

void remote_disk_access::cleanup()
{
	for (int i = 0; i < num_senders; i++) {
		msg_sender<io_request> *sender = senders[i];
		sender->flush_all();
	}
	int num;
	do {
		num = 0;
		for (int i = 0; i < num_senders; i++) {
			num += queues[i]->get_num_entries();
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
	for (int i = 0; i < num; i++) {
		assert(requests[i].get_size() > 0);
		// TODO data is striped on disks, we have to make sure the data
		// can be accessed from the remote disk access.
		// and I assume the request size is aligned with the strip size.
		off_t pg_off = requests[i].get_offset() / PAGE_SIZE;
		int idx = file_hash(pg_off, num_senders);
		int ret = senders[idx]->send_cached(&requests[i]);
		/* send should always succeed. */
		assert(ret > 0);
	}
	return 0;
}
