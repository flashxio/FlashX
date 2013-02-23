#ifndef __DISK_READ_THREAD_H__
#define __DISK_READ_THREAD_H__

#include <string>

#include "messaging.h"
#include "aio_private.h"
#include "container.h"


void *process_requests(void *arg);

class disk_read_thread
{
	blocking_FIFO_queue<io_request> queue;
	pthread_t id;
	async_io *aio;
	int node_id;

public:
	disk_read_thread(const char *name, long size, int node_id);

	blocking_FIFO_queue<io_request> *get_queue() {
		return &queue;
	}

	~disk_read_thread() {
		delete aio;
	}

	void run();
};

#endif
