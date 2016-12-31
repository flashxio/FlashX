#ifndef __COMP_ACCESS_H__
#define __COMP_ACCESS_H__

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

#include "container.h"

#include "io_interface.h"

namespace safs
{

class remote_io;
class comp_io_scheduler;
class direct_byte_array_allocator;

/**
 * This class supports asynchronous user-task I/O requests.
 * Unlike global_cached_io, it doesn't support page cache.
 */
class direct_comp_io: public io_interface
{
	size_t num_disk_bytes;
	size_t num_req_bytes;
	size_t num_issued_areqs;
	size_t num_completed_areqs;

	// The memory that have been allocated for the I/O requests that have
	// been issued to the underlying I/O.
	size_t alloc_mem_size;
	fifo_queue<io_request> req_buf;
	std::shared_ptr<remote_io> underlying;
	std::shared_ptr<comp_io_scheduler> comp_sched;
	std::unique_ptr<direct_byte_array_allocator> arr_alloc;

	void process_buf_reqs();
	void process_incomplete_computes();
public:
	direct_comp_io(std::shared_ptr<remote_io> io);
	~direct_comp_io();

	/*
	 * This method is called when an I/O request is completed by
	 * the underlying I/O.
	 */
	void complete_req(const io_request &req);

	virtual int get_file_id() const;
	virtual void cleanup();
	virtual bool support_aio() {
		return true;
	}

	virtual void access(io_request *requests, int num,
			io_status *status = NULL);
	virtual void flush_requests();
	virtual int wait4complete(int num);
	virtual int num_pending_ios() const;

	size_t get_num_disk_bytes() {
		return num_disk_bytes;
	}

	size_t get_num_req_bytes() {
		return num_req_bytes;
	}

	size_t get_num_reqs() {
		return num_issued_areqs;
	}

	size_t get_num_completed_reqs() {
		return num_completed_areqs;
	}
};

}

#endif
