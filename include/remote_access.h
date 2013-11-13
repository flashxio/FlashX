#ifndef __REMOTE_ACCESS_H__
#define __REMOTE_ACCESS_H__

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

#include "io_interface.h"
#include "container.h"

class request_sender;
class disk_io_thread;
class slab_allocator;
class file_mapper;

const int COMPLETE_QUEUE_SIZE = 10240;

/**
 * This class is to help the local thread send IO requests to remote threads
 * dedicated to accessing SSDs. Each SSD has such a thread.
 * However, the helper class isn't thread safe, so each local thread has to
 * reference its own helper object.
 */
class remote_io: public io_interface
{
	static atomic_integer num_ios;
	const int max_disk_cached_reqs;
	// They work as buffers for requests and are only used to
	// send high-priority requests.
	std::vector<request_sender *> senders;
	// They are used to send low-priority requests.
	std::vector<request_sender *> low_prio_senders;
	std::vector<disk_io_thread *> io_threads;
	callback *cb;
	file_mapper *block_mapper;
	thread_safe_FIFO_queue<io_request> complete_queue;
	slab_allocator *msg_allocator;

	atomic_integer num_completed_reqs;
	atomic_integer num_issued_reqs;
public:
	remote_io(const std::vector<disk_io_thread *> &remotes,
			file_mapper *mapper, thread *t,
			int max_reqs = MAX_DISK_CACHED_REQS);

	~remote_io();

	virtual int process_completed_requests(io_request reqs[], int num);
	int process_completed_requests(int num);

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

	virtual int get_file_id() const;

	virtual void access(io_request *requests, int num,
			io_status *status = NULL);
	virtual void notify_completion(io_request *reqs[], int num);
	virtual int wait4complete(int num_to_complete);
	virtual int num_pending_ios() const {
		return num_issued_reqs.get() - num_completed_reqs.get();
	}
	virtual int get_max_num_pending_ios() const;
	virtual io_interface *clone(thread *t) const;
	void flush_requests(int max_cached);
	virtual void flush_requests();
	virtual void print_stat(int nthreads) {
		printf("remote_io: %d reqs, %d completed reqs\n",
				num_issued_reqs.get(), num_completed_reqs.get());
	}
	virtual void print_state();
};

#endif
