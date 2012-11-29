#include "remote_access.h"
#include "parameters.h"

remote_disk_access::remote_disk_access(disk_read_thread **remotes,
			int num_remotes)
{
	senders = new msg_sender<io_request> *[num_remotes];
	queues = new thread_safe_FIFO_queue<io_request> *[num_remotes];
	num_senders = num_remotes;
	// create a msg sender for each disk read thread.
	for (int i = 0; i < num_remotes; i++) {
		thread_safe_FIFO_queue<io_request> *queue = remotes[i]->get_queue();
		senders[i] = new msg_sender<io_request>(MSG_SEND_BUF_SIZE,
				&queue, 1);
		queues[i] = queue;
	}
}

remote_disk_access::~remote_disk_access()
{
	for (int i = 0; i < num_senders; i++)
		delete senders[i];
	delete [] senders;
	delete [] queues;
}

void remote_disk_access::cleanup()
{
	for (int i = 0; i < num_senders; i++) {
		msg_sender<io_request> *sender = senders[i];
		while (sender->num_msg())
			sender->flush();
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
		// TODO I assume data is striped on disks.
		// and I assume the request size is aligned with the strip size.
		int idx = requests[i].get_offset() / PAGE_SIZE % num_senders;
		int ret = senders[idx]->send_cached(&requests[i]);
		/* send should always succeed. */
		assert(ret > 0);
	}
	return 0;
}
