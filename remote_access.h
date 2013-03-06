#ifndef __REMOTE_ACCESS_H__
#define __REMOTE_ACCESS_H__

#include "disk_read_thread.h"
#include "messaging.h"
#include "container.h"

/**
 * This class is to help the local thread send IO requests to remote threads
 * dedicated to accessing SSDs. Each SSD has such a thread.
 * However, the helper class isn't thread safe, so each local thread has to
 * reference its own helper object.
 */
class remote_disk_access: public io_interface
{
	msg_sender<io_request> **senders;
	thread_safe_FIFO_queue<io_request> **queues;
	int num_senders;
	callback *cb;

	// The total size of data accessible with this IO interface.
	long total_size;
	// The size of data on the local node.
	long local_size;

	remote_disk_access(int node_id): io_interface(node_id) {
		senders = NULL;
		queues = NULL;
		num_senders = 0;
		cb = NULL;
	}
public:
	remote_disk_access(disk_read_thread **remotes,
			int num_remotes, int node_id);

	~remote_disk_access();

	virtual int thread_init() {
		return 0;
	}

	virtual bool support_aio() {
		return true;
	}

	virtual void cleanup();

	virtual bool set_callback(callback *cb) {
		this->cb = cb;
		return true;
	}

	virtual callback *get_callback() {
		return cb;
	}

	long get_size() const {
		return total_size;
	}

	long get_local_size() const {
		return local_size;
	}

	virtual ssize_t access(io_request *requests, int num);
	virtual io_interface *clone() const;
};

#endif
