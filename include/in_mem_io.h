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

#include "NUMA_mapper.h"

#include "io_interface.h"
#include "comp_io_scheduler.h"
#include "cache.h"

namespace safs
{

class NUMA_buffer
{
	std::vector<std::shared_ptr<char> > bufs;
	// This has the length of individual buffers.
	// The sum of the physical buffer lengths >= the total length;
	std::vector<size_t> buf_lens;
	// This is the total length of the buffer.
	size_t length;
	NUMA_mapper mapper;

	struct data_loc_info {
		int node_id;
		off_t local_off;
		size_t local_size;
		data_loc_info(int node_id, off_t local_off, size_t local_size) {
			this->node_id = node_id;
			this->local_off = local_off;
			this->local_size = local_size;
		}
	};
	data_loc_info get_data_loc(off_t off, size_t size) const;

	NUMA_buffer(size_t length, const NUMA_mapper &mapper);
public:
	typedef std::pair<const char *, size_t> cdata_info;
	typedef std::pair<char *, size_t> data_info;
	typedef std::shared_ptr<NUMA_buffer> ptr;

	/*
	 * Load data in a file to the buffer.
	 */
	static ptr load(const std::string &file, const NUMA_mapper &mapper);
	static ptr load_safs(const std::string &file, const NUMA_mapper &mapper);

	static ptr create(size_t length, const NUMA_mapper &mapper) {
		return ptr(new NUMA_buffer(length, mapper));
	}

	size_t get_length() const {
		return length;
	}

	/*
	 * Get the data in the specified location.
	 * Since the data in the buffer isn't stored contiguously, the size of
	 * the returned data may be smaller than the specified size.
	 */
	cdata_info get_data(off_t off, size_t size) const;
	data_info get_data(off_t off, size_t size);

	/*
	 * Write the data in the given buffer to the specified location.
	 */
	void copy_from(const char *buf, size_t size, off_t off);
	/*
	 * Copy the data in the specified location to the given buffer.
	 */
	void copy_to(char *buf, size_t size, off_t off) const;

	void dump(const std::string &file);
};

/*
 * This class provides a single I/O interface for accessing data in memory.
 * The main reason of having this class is to allow us to reuse a lot of code
 * when data is stored in memory.
 */
class in_mem_io: public io_interface
{
	NUMA_buffer::ptr data;
	int file_id;
	fifo_queue<io_request> req_buf;
	comp_io_scheduler::ptr comp_io_sched;
	std::unique_ptr<byte_array_allocator> array_allocator;

	callback::ptr cb;

	void process_req(const io_request &req);
	void process_computes();
public:
	in_mem_io(NUMA_buffer::ptr data, int file_id, thread *t);

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
			int access_method);
	virtual void access(io_request *requests, int num, io_status *status);
	virtual int wait4complete(int) {
		return 0;
	}
};

class in_mem_io_factory: public file_io_factory
{
	// The I/O interface doesn't own the byte array.
	NUMA_buffer::ptr data;
	int file_id;
public:
	in_mem_io_factory(NUMA_buffer::ptr data, int file_id,
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
