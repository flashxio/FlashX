#ifndef __IO_INTERFACE_H__
#define __IO_INTERFACE_H__

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

#include <stdlib.h>

#include <vector>
#include <memory>

#include "exception.h"
#include "common.h"
#include "concurrency.h"
#include "thread.h"
#include "io_request.h"
#include "comp_io_scheduler.h"

class io_request;

class callback
{
public:
	virtual ~callback() {
	}

	virtual int invoke(io_request *reqs[], int num) = 0;
};

class io_status
{
	long status: 8;
	long priv_data: 56;
public:
	io_status() {
		status = 0;
		priv_data = 0;
	}

	io_status(int status) {
		this->status = status;
		priv_data = 0;
	}

	void set_priv_data(long data) {
		priv_data = data;
	}

	long get_priv_data() const {
		return priv_data;
	}

	io_status &operator=(int status) {
		this->status = status;
		return *this;
	}

	bool operator==(int status) {
		return this->status == status;
	}
};

enum
{
	IO_OK,
	IO_PENDING = -1,
	IO_FAIL = -2,
	IO_UNSUPPORTED = -3,
};

enum {
	READ_ACCESS,
	DIRECT_ACCESS,
	AIO_ACCESS,
	REMOTE_ACCESS,
	GLOBAL_CACHE_ACCESS,
	PART_GLOBAL_ACCESS,
};

/**
 * The interface for all IO classes.
 */
class io_interface
{
	thread *curr;		// the thread where the IO instance runs.

	// This is an index for locating this IO object in a global table.
	int io_idx;
	int max_num_pending_ios;
	static atomic_integer io_counter;

public:
	io_interface(thread *t) {
		this->curr = t;
		this->io_idx = io_counter.inc(1) - 1;
		max_num_pending_ios = 32;
	}

	virtual ~io_interface() { }

	thread *get_thread() const {
		assert(curr);
		return curr;
	}

	int get_node_id() const {
		assert(curr);
		return curr->get_node_id();
	}

	int get_io_idx() const {
		return io_idx;
	}

	/**
	 * The number of requests are allowed to be sent to the IO instance.
	 */
	int get_remaining_io_slots() const {
		return get_max_num_pending_ios() - num_pending_ios();
	}

	/* When a thread begins, this method will be called. */
	virtual int init() {
		return 0;
	}

	/**
	 * Each IO instance can only access one file. This returns the ID of
	 * the file being accessed by the IO instance.
	 */
	virtual int get_file_id() const = 0;

	/**
	 * This method needs to be called before destroying the IO instance.
	 */
	virtual void cleanup() {
	}

	/**
	 * This method prints the final statistics of the IO instance.
	 */
	virtual void print_stat(int nthreads) {
	}

	/**
	 * This method prints the current state of the IO instance.
	 * For the sake of performance, the method may not be thread-safe.
	 * It should be used with caution.
	 */
	virtual void print_state() {
	}

	virtual io_interface *clone(thread *t) const {
		return NULL;
	}

	/**
	 * Indicate whether it supports asynchronous IO interface.
	 */
	virtual bool support_aio() {
		return false;
	}

	/**
	 * The asynchronous IO should implement some of the following methods.
	 */

	/**
	 * The main interface to send requests.
	 */
	virtual void access(io_request *requests, int num, io_status *status = NULL) {
		throw unsupported_exception();
	}
	/**
	 * When requests are passed to the access method, an IO layer may buffer
	 * the requests. This method guarantees that all requests are flushed to
	 * the underlying devices.
	 */
	virtual void flush_requests() {
		throw unsupported_exception();
	}
	/**
	 * This method waits for at least the specified number of requests currently
	 * being sent by the access method to complete.
	 */
	virtual int wait4complete(int num) {
		throw unsupported_exception();
	}
	virtual int num_pending_ios() const {
		throw unsupported_exception();
	}
	virtual int get_max_num_pending_ios() const {
		return max_num_pending_ios;
	}
	virtual void set_max_num_pending_ios(int max) {
		this->max_num_pending_ios = max;
	}

	/**
	 * This method gives the underlying layer an interface to notify
	 * the current IO of the completed requests.
	 * The method may be called by multiple threads, so it has to be made
	 * thread-safe.
	 * The requests should be issued by this IO.
	 */
	virtual void notify_completion(io_request *reqs[], int num) {
		if (get_callback())
			get_callback()->invoke(reqs, num);
	}

	/**
	 * set the callback if the class supports the asynchronous fashion.
	 * If the class doesn't support async IO, return false.
	 */
	virtual bool set_callback(callback *cb) {
		return false;
	}

	virtual callback *get_callback() {
		return NULL;
	}

	/**
	 * The synchronous IO interface.
	 */
	virtual io_status access(char *, off_t, ssize_t, int) {
		return IO_UNSUPPORTED;
	}
};

class comp_io_sched_creater
{
public:
	virtual comp_io_scheduler *create(int node_id) const = 0;
};

/**
 * The interface of creating IOs to access a file.
 */
class file_io_factory
{
	comp_io_sched_creater *creater;
	// The name of the file.
	const std::string name;
public:
	typedef std::shared_ptr<file_io_factory> shared_ptr;

	file_io_factory(const std::string _name): name(_name) {
		creater = NULL;
	}

	virtual ~file_io_factory() {
	}

	void set_sched_creater(comp_io_sched_creater *creater) {
		this->creater = creater;
	}

	comp_io_sched_creater *get_sched_creater() const {
		return creater;
	}

	const std::string &get_name() const {
		return name;
	}

	virtual int get_file_id() const = 0;

	virtual io_interface *create_io(thread *t) = 0;
	virtual void destroy_io(io_interface *) = 0;
	virtual void print_state() {
	}

	ssize_t get_file_size() const;
};

class cache_config;
class RAID_config;

file_io_factory::shared_ptr create_io_factory(const std::string &file_name,
		const int access_option);

void init_io_system(const config_map &map);
void destroy_io_system();
const RAID_config &get_sys_RAID_conf();

// This interface is used for debugging.
void print_io_thread_stat();

#endif
