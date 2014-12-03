#ifndef __REMOTE_ACCESS_H__
#define __REMOTE_ACCESS_H__

/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "slab_allocator.h"
#include "io_interface.h"
#include "container.h"

namespace safs
{

class request_sender;
class disk_io_thread;
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
	std::vector<std::shared_ptr<disk_io_thread> > io_threads;
	callback::ptr cb;
	file_mapper *block_mapper;
	thread_safe_FIFO_queue<io_request> complete_queue;
	slab_allocator &msg_allocator;

	atomic_integer num_completed_reqs;
	atomic_integer num_issued_reqs;
public:
	remote_io(const std::vector<std::shared_ptr<disk_io_thread> > &remotes,
			slab_allocator &msg_allocator, file_mapper *mapper, thread *t,
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

	virtual bool set_callback(callback::ptr cb) {
		this->cb = cb;
		return true;
	}

	virtual bool have_callback() const {
		return cb != NULL;
	}

	virtual callback &get_callback() {
		return *cb;
	}

	virtual int get_file_id() const;

	virtual void access(io_request *requests, int num,
			io_status *status = NULL);
	virtual void notify_completion(io_request *reqs[], int num);
	virtual int wait4complete(int num_to_complete);
	virtual int num_pending_ios() const {
		return num_issued_reqs.get() - num_completed_reqs.get();
	}
	virtual io_interface *clone(thread *t) const;
	void flush_requests(int max_cached);
	virtual void flush_requests();
	virtual void print_state();

	size_t get_num_reqs() const {
		return num_issued_reqs.get();
	}
};

}

#endif
