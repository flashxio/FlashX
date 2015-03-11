#ifndef __IN_MEM_IO_H__
#define __IN_MEM_IO_H__

/*
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

#include "io_interface.h"
#include "comp_io_scheduler.h"
#include "cache.h"

namespace safs
{

/*
 * This class provides a single I/O interface for accessing data in memory.
 * The main reason of having this class is to allow us to reuse a lot of code
 * when data is stored in memory.
 */
class in_mem_io: public io_interface
{
	// The I/O interface doesn't own the byte array.
	// TODO we should use a smart pointer here.
	char *data;
	int file_id;
	fifo_queue<io_request> req_buf;
	comp_io_scheduler::ptr comp_io_sched;
	std::unique_ptr<byte_array_allocator> array_allocator;

	callback::ptr cb;

	void process_req(const io_request &req);
	void process_computes();
public:
	in_mem_io(char *data, int file_id, thread *t);

	virtual int get_file_id() const {
		return file_id;
	}

	virtual bool support_aio() {
		return true;
	}

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

	virtual void flush_requests() { }

	virtual int num_pending_ios() const {
		return 0;
	}

	virtual io_status access(char *buf, off_t off, ssize_t size,
			int access_method) {
		assert(access_method == READ);
		memcpy(buf, data + off, size);
		return IO_OK;
	}

	virtual void access(io_request *requests, int num, io_status *status);
	virtual int wait4complete(int) {
		return 0;
	}
};

class in_mem_io_factory: public file_io_factory
{
	// The I/O interface doesn't own the byte array.
	// TODO we should use a smart pointer here.
	char *data;
	int file_id;
public:
	in_mem_io_factory(char *data, int file_id,
			const std::string file_name): file_io_factory(file_name) {
		this->data = data;
		this->file_id = file_id;
	}

	virtual int get_file_id() const {
		return file_id;
	}

	virtual io_interface::ptr create_io(thread *t) {
		return io_interface::ptr(new in_mem_io(data, file_id, t));
	}

	virtual void destroy_io(io_interface &) {
	}
};

}

#endif
