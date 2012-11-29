#ifndef __REMOTE_ACCESS_H__
#define __REMOTE_ACCESS_H__

#include "disk_read_thread.h"
#include "messaging.h"
#include "container.h"

class remote_disk_access: public io_interface
{
	msg_sender<io_request> **senders;
	thread_safe_FIFO_queue<io_request> **queues;
	int num_senders;
	callback *cb;
public:
	remote_disk_access(disk_read_thread **remotes,
			int num_remotes);

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

	virtual ssize_t access(io_request *requests, int num);
};

#endif
