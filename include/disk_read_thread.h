#ifndef __DISK_READ_THREAD_H__
#define __DISK_READ_THREAD_H__

#include <unistd.h>

#include <string>
#include <tr1/unordered_map>

#include "aio_private.h"
#include "io_request.h"
#include "container.h"
#include "file_partition.h"
#include "messaging.h"
#include "cache.h"

void *process_requests(void *arg);

class async_io;

class disk_read_thread
{
	static const int LOCAL_BUF_SIZE = 16;

	const int disk_id;

	msg_queue<io_request> queue;
	msg_queue<io_request> low_prio_queue;
	logical_file_partition partition;
	std::vector<file_mapper *> open_files;

	pthread_t id;
	async_io *aio;
	int node_id;
#ifdef STATISTICS
	long num_reads;
	long num_writes;
	long num_read_bytes;
	long num_write_bytes;
	long num_low_prio_accesses;
	long num_requested_flushes;
	long num_ignored_flushes_evicted;
	long num_ignored_flushes_cleaned;
	long num_ignored_flushes_old;
	long tot_flush_delay;	// in us
	long max_flush_delay;
	long min_flush_delay;
#endif

	volatile bool running;

	atomic_integer flush_counter;

	class dirty_page_filter: public page_filter {
		const std::vector<file_mapper *> &mappers;
		int disk_id;
	public:
		dirty_page_filter(const std::vector<file_mapper *> &_mappers,
				int disk_id): mappers(_mappers) {
			this->disk_id = disk_id;
		}

		int filter(const thread_safe_page *pages[], int num,
				const thread_safe_page *returned_pages[]);
	};

	page_cache *cache;
	dirty_page_filter filter;

	int process_low_prio_msg(message<io_request> &low_prio_msg);

	int get_num_high_prio_reqs() {
		return queue.get_num_objs();
	}

	int get_num_low_prio_reqs() {
		return low_prio_queue.get_num_objs();
	}

public:
	disk_read_thread(const logical_file_partition &partition, int node_id,
			page_cache *cache, int disk_id);

	msg_queue<io_request> *get_queue() {
		return &queue;
	}

	msg_queue<io_request> *get_low_prio_queue() {
		return &low_prio_queue;
	}

	int get_node_id() const {
		return node_id;
	}

	const std::string get_file_name() const {
		if (open_files.empty())
			return "";

		logical_file_partition *part = partition.create_file_partition(open_files[0]);
		std::string name = part->get_file_name(0);
		delete part;
		return name;
	}

	/**
	 * Flush threads asynchronously.
	 * The invoker of this function shouldn't be the I/O thread.
	 * So we need to wake up the I/O thread and notify the I/O thread
	 * to flush requests.
	 */
	void flush_requests() {
		flush_counter.inc(1);
		// If the I/O thread is blocked by the request queue, we should
		// wake the thread up.
		queue.wakeup();
	}

	// It open a new file. The mapping is still the same.
	int open_file(file_mapper *mapper) {
		open_files.push_back(mapper);
		logical_file_partition *part = partition.create_file_partition(mapper);
		int ret = aio->open_file(*part);
		delete part;
		return ret;
	}

	void register_cache(page_cache *cache) {
		this->cache = cache;
	}

	~disk_read_thread() {
		delete aio;
	}

	void run();

	void stop() {
		running = false;
	}

	void print_stat() {
#ifdef STATISTICS
		printf("queue on file %s wait for requests for %d times, is full for %d times,\n",
				get_file_name().c_str(), get_queue()->get_num_empty(),
				get_queue()->get_num_full());
		printf("\t%ld reads (%ld bytes), %ld writes (%ld bytes) and %d io waits, complete %d reqs and %ld low-prio reqs,\n",
				num_reads, num_read_bytes, num_writes, num_write_bytes, aio->get_num_iowait(), aio->get_num_completed_reqs(),
				num_low_prio_accesses);
		printf("\trequest %ld flushes, ignore flushes: %ld evicted, %ld cleaned, %ld out-of-date\n",
				num_requested_flushes, num_ignored_flushes_evicted,
				num_ignored_flushes_cleaned, num_ignored_flushes_old);
		if (num_low_prio_accesses > 0)
			printf("\tavg flush delay: %ldus, max flush delay: %ldus, min flush delay: %ldus\n",
					tot_flush_delay / num_low_prio_accesses, max_flush_delay,
					min_flush_delay);
		printf("\tremain %d high-prio requests, %d low-prio requests\n",
				get_num_high_prio_reqs(), get_num_low_prio_reqs());
#endif
	}
};

#endif
