#ifndef __REMOTE_ACCESS_H__
#define __REMOTE_ACCESS_H__

#include "disk_read_thread.h"
#include "messaging.h"

class remote_disk_access: public io_interface
{
	msg_sender<io_request> *sender;
	bulk_queue<io_request> **queues;
	int num_queues;
	callback *cb;
public:
	remote_disk_access(disk_read_thread **remotes,
			int num_remotes) {
		queues = new bulk_queue<io_request> *[num_remotes];
		num_queues = num_remotes;
		for (int i = 0; i < num_remotes; i++) {
			queues[i] = remotes[i]->get_queue();
		}
		sender = new msg_sender<io_request>(BUF_SIZE, queues, num_remotes);
	}

	~remote_disk_access() {
		delete sender;
		delete queues;
	}

	virtual int thread_init() {
		return 0;
	}

	virtual bool support_aio() {
		return true;
	}

	virtual void cleanup() {
		while (sender->num_msg())
			sender->flush();
		int num;
		do {
			num = 0;
			for (int i = 0; i < num_queues; i++) {
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

	virtual bool set_callback(callback *cb) {
		this->cb = cb;
		return true;
	}

	virtual callback *get_callback() {
		return cb;
	}

	virtual ssize_t access(io_request *requests, int num) {
		for (int i = 0; i < num; i++) {
			int ret = sender->send_cached(&requests[i]);
			/* send should always succeed. */
			assert(ret > 0);
		}
		return 0;
	}
};

#endif
