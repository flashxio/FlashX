#ifndef __REMOTE_ACCESS_H__
#define __REMOTE_ACCESS_H__

#include "disk_read_thread.h"
#include "messaging.h"

class remote_disk_access: public thread_private
{
	msg_sender<io_request> *sender;
public:
	remote_disk_access(disk_read_thread **remotes,
			int num_remotes): thread_private(0, 0) {
		bulk_queue<io_request> *queues[num_remotes];
		for (int i = 0; i < num_remotes; i++) {
			queues[i] = remotes[i]->get_queue();
		}
		sender = new msg_sender<io_request>(BUF_SIZE, queues, num_remotes);
	}

	~remote_disk_access() {
		delete sender;
	}

	virtual int thread_init() {
		return 0;
	}

	virtual bool support_bulk() {
		return true;
	}

	virtual ssize_t access(io_request *requests, int num, int access_method) {
		for (int i = 0; i < num; i++) {
			int ret = sender->send_cached(&requests[i]);
			/* send should always succeed. */
			assert(ret > 0);
		}
		return 0;
	}
};

#endif
