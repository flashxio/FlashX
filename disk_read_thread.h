#ifndef __DISK_READ_THREAD_H__
#define __DISK_READ_THREAD_H__

#include <string>

#include "messaging.h"
#include "aio_private.h"
#include "container.h"
#include "file_partition.h"

void *process_requests(void *arg);

class disk_read_thread
{
	blocking_FIFO_queue<io_request> queue;
	pthread_t id;
	async_io *aio;
	int node_id;
	int num_accesses;

public:
	disk_read_thread(const logical_file_partition &partition, long size,
			int node_id);

	blocking_FIFO_queue<io_request> *get_queue() {
		return &queue;
	}

	long get_size() const {
		return aio->get_size();
	}

	int get_node_id() const {
		return node_id;
	}

	int get_num_accesses() const {
		return num_accesses;
	}

	int get_num_iowait() const {
		return aio->get_num_iowait();
	}

	~disk_read_thread() {
		delete aio;
	}

	void run();
};

#endif
