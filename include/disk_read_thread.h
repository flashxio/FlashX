#ifndef __DISK_READ_THREAD_H__
#define __DISK_READ_THREAD_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <unistd.h>

#include <string>
#include <tr1/unordered_map>

#include "aio_private.h"
#include "io_request.h"
#include "container.h"
#include "file_partition.h"
#include "messaging.h"
#include "cache.h"
#include "thread.h"

void *process_requests(void *arg);

class async_io;

class disk_io_thread: public thread
{
	static const int LOCAL_BUF_SIZE = 16;

	/**
	 * This is a remote command that is to be executed in the I/O thread.
	 * It is mainly used to open and close physical files in the I/O thread.
	 */
	class remote_comm
	{
		pthread_mutex_t mutex;
		pthread_cond_t cond;
		bool is_complete;
		int status;
	public:
		remote_comm() {
			status = 0;
			is_complete = false;
			pthread_mutex_init(&mutex, NULL);
			pthread_cond_init(&cond, NULL);
		}

		virtual ~remote_comm() {
		}

		virtual void run() = 0;

		void complete() {
			pthread_mutex_lock(&mutex);
			is_complete = true;
			pthread_mutex_unlock(&mutex);
			pthread_cond_signal(&cond);
		}

		void wait4complete() {
			pthread_mutex_lock(&mutex);
			while (!is_complete) {
				pthread_cond_wait(&cond, &mutex);
			}
			pthread_mutex_unlock(&mutex);
		}

		void set_status(int status) {
			this->status = status;
		}

		int get_status() const {
			return status;
		}
	};

	/**
	 * The command opens a file in the I/O thread.
	 */
	class open_comm: public remote_comm
	{
		file_mapper *mapper;
		const logical_file_partition *partition;
		async_io *aio;
	public:
		open_comm(async_io *aio, file_mapper *mapper,
				const logical_file_partition *partition) {
			this->aio = aio;
			this->mapper = mapper;
			this->partition = partition;
		}

		void run() {
			logical_file_partition *part = partition->create_file_partition(
					mapper);
			int ret = aio->open_file(*part);
			delete part;
			set_status(ret);
		}
	};

	class close_comm: public remote_comm
	{
		async_io *aio;
		int file_id;
	public:
		close_comm(async_io *aio, int file_id) {
			this->aio = aio;
			this->file_id = file_id;
		}

		void run() {
			int ret = aio->close_file(file_id);
			set_status(ret);
		}
	};

	const int disk_id;

	msg_queue<io_request> queue;
	msg_queue<io_request> low_prio_queue;
	thread_safe_FIFO_queue<remote_comm *> comm_queue;
	logical_file_partition partition;

	async_io *aio;
#ifdef STATISTICS
	long num_empty;
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
	long num_msgs;
#endif

	atomic_integer flush_counter;

	class dirty_page_filter: public page_filter {
		const file_mapper *mapper;
		int disk_id;
	public:
		dirty_page_filter(const file_mapper *_mapper, int disk_id) {
			this->disk_id = disk_id;
			this->mapper = _mapper;
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

	void run_commands(thread_safe_FIFO_queue<remote_comm *> &);

	int execute_remote_comm(remote_comm *comm) {
		comm_queue.add(&comm, 1);
		// We need to wake up the I/O thread to run the command.
		this->activate();
		// And then wait for the command to be executed.
		comm->wait4complete();
		int ret = comm->get_status();
		delete comm;
		return ret;
	}
public:
	disk_io_thread(const logical_file_partition &partition, int node_id,
			page_cache *cache, int disk_id);

	msg_queue<io_request> *get_queue() {
		return &queue;
	}

	msg_queue<io_request> *get_low_prio_queue() {
		return &low_prio_queue;
	}

	/**
	 * Flush threads asynchronously.
	 * The invoker of this function shouldn't be the I/O thread.
	 * So we need to wake up the I/O thread and notify the I/O thread
	 * to flush requests.
	 */
	void flush_requests() {
		flush_counter.inc(1);
		activate();
	}

	// It open a new file. The mapping is still the same.
	int open_file(file_mapper *mapper) {
		remote_comm *comm = new open_comm(aio, mapper, &partition);
		return execute_remote_comm(comm);
	}

	int close_file(file_mapper *mapper) {
		remote_comm *comm = new close_comm(aio, mapper->get_file_id());
		return execute_remote_comm(comm);
	}

	void register_cache(page_cache *cache) {
		this->cache = cache;
	}

	~disk_io_thread() {
		delete aio;
	}

	void run();
	void init() {
		aio->init();
	}
	virtual void cleanup() {
		aio->cleanup();
	}

	void print_stat() {
#ifdef STATISTICS
		printf("queue on disk %d wait for requests for %ld times,\n",
				disk_id, num_empty);
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
		printf("\tremain %d high-prio requests, %d low-prio requests, %ld messages in total\n",
				get_num_high_prio_reqs(), get_num_low_prio_reqs(), num_msgs);
		aio->print_stat(1);
#endif
	}

	void print_state();
};

#endif
