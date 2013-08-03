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
	blocking_FIFO_queue<io_request> low_prio_queue;
	logical_file_partition partition;
	std::vector<file_mapper *> open_files;

	pthread_t id;
	async_io *aio;
	int node_id;
	int num_accesses;
	int num_low_prio_accesses;

public:
	disk_read_thread(const logical_file_partition &partition,
			const std::tr1::unordered_map<int, aio_complete_thread *> &complete_threads,
			int node_id);

	blocking_FIFO_queue<io_request> *get_queue() {
		return &queue;
	}

	blocking_FIFO_queue<io_request> *get_low_prio_queue() {
		return &low_prio_queue;
	}

	int get_node_id() const {
		return node_id;
	}

	int get_num_accesses() const {
		return num_accesses;
	}

	int get_num_low_prio_accesses() const {
		return num_low_prio_accesses;
	}

	int get_num_iowait() const {
		return aio->get_num_iowait();
	}

	int get_num_completed_reqs() const {
		return aio->get_num_completed_reqs();
	}

	int get_num_local_alloc() const {
		return aio->get_num_local_alloc();
	}

	const std::string get_file_name() const {
		logical_file_partition *part = partition.create_file_partition(open_files[0]);
		std::string name = part->get_file_name(0);
		delete part;
		return name;
	}

	// It open a new file. The mapping is still the same.
	int open_file(file_mapper *mapper) {
		open_files.push_back(mapper);
		logical_file_partition *part = partition.create_file_partition(mapper);
		int ret = aio->open_file(*part);
		delete part;
		return ret;
	}

	~disk_read_thread() {
		delete aio;
	}

	void run();
};

#endif
