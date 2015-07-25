#ifndef __AIO_PRIVATE_H__
#define __AIO_PRIVATE_H__

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

#include <deque>
#include <tr1/unordered_map>

#include "wpaio.h"
#include "io_interface.h"
#include "thread.h"
#include "container.h"
#include "io_request.h"

namespace safs
{

void aio_callback(io_context_t, struct iocb*, void *, long, long);

struct thread_callback_s;

class buffered_io;
class logical_file_partition;
class callback_allocator;
class virt_data_impl;

class async_io: public io_interface
{
	int buf_idx;
	aio_ctx *ctx;
	callback::ptr cb;
	const int AIO_DEPTH;
	callback_allocator *cb_allocator;
	int open_flags;

	int num_iowait;
	int num_completed_reqs;

	class io_ref
	{
		std::shared_ptr<buffered_io> io;
		int count;
	public:
		io_ref() {
			this->count = 0;
		}

		io_ref(buffered_io *io);

		void inc_ref() {
			count++;
		}

		void dec_ref();

		int get_count() const {
			return count;
		}

		const buffered_io &get_io() const {
			return *io;
		}

		buffered_io &get_io() {
			return *io;
		}

		bool is_valid() const {
			return io != NULL;
		}
	};
	// file id <-> buffered io
	std::tr1::unordered_map<int, io_ref> open_files;
	io_ref default_io;

#if 0
	virt_data_impl *data;
#endif

	struct iocb *construct_req(io_request &io_req, callback_t cb_func);
public:
	/**
	 * @aio_depth_per_file
	 * @node_id: the NUMA node where the disks to be read are connected to.
	 */
	async_io(const logical_file_partition &partition,
			int aio_depth_per_file, thread *t, const safs_header &header,
			int flags);

	virtual ~async_io();

	virtual io_status access(char *, off_t, ssize_t, int) {
		return IO_UNSUPPORTED;
	}
	virtual void access(io_request *requests, int num, io_status *status = NULL);

	bool set_callback(callback::ptr cb) {
		this->cb = cb;
		return true;
	}

	bool have_callback() const {
		return this->cb != NULL;
	}

	callback &get_callback() {
		return *cb;
	}

	bool support_aio() {
		return true;
	}

	int get_file_id() const;

	virtual void cleanup();

	void return_cb(thread_callback_s *tcbs[], int num);

	int num_available_IO_slots() const {
		return ctx->max_io_slot();
	}

	virtual int num_pending_ios() const {
		return AIO_DEPTH - ctx->max_io_slot();
	}

	virtual void notify_completion(io_request *reqs[], int num);
	int wait4complete(int num) {
		return ctx->io_wait(NULL, num);
	}
	virtual int get_max_num_pending_ios() const {
		return AIO_DEPTH;
	}
	virtual void set_max_num_pending_ios(int max) {
	}

	int get_num_iowait() const {
		return num_iowait;
	}

	int get_num_completed_reqs() const {
		return num_completed_reqs;
	}

	virtual void flush_requests();

	// These two interfaces allow users to open and close more files.
	
	/*
	 * It opens a virtual file.
	 * Actually, it opens physical files on the underlying filesystems
	 * within the partition of the virtual file, managed by the IO interface.
	 */
	int open_file(const logical_file_partition &partition);
	int close_file(int file_id);

	virtual void print_state() {
		printf("aio %d has %ld open files, %d pending reqs\n",
				get_io_id(), open_files.size(), num_pending_ios());
	}
};

void init_aio(std::vector<int> node_ids);
void destroy_aio();

}

#endif
