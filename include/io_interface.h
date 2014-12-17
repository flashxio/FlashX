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

#include "config_map.h"
#include "exception.h"
#include "common.h"
#include "concurrency.h"
#include "thread.h"
#include "io_request.h"
#include "comp_io_scheduler.h"

namespace safs
{

class io_request;

/**
 * The callback interface to notify the completion of I/O requests.
 */
class callback
{
public:
	typedef std::shared_ptr<callback> ptr;

	virtual ~callback() {
	}

	/**
	 * The user-defined code is invoked on the completed I/O requests.
	 * \param reqs an array of completed I/O requests.
	 * \param num the number of I/O requests in the array.
	 */
	virtual int invoke(io_request *reqs[], int num) = 0;
};

/**
 * The I/O status for an I/O request.
 * Right now, it's not used.
 */
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

/**
 * This defines the method of accessing a SAFS file. There are six options.
 */
enum {
	/*
	 * The three methods below are used internally or for testing only.
	 */
	READ_ACCESS,
	DIRECT_ACCESS,
	AIO_ACCESS,

	/**
	 * This method is equivalent to Linux direct I/O. It accesses data
	 * in a SAFS file without caching. It only supports asynchronous I/O.
	 */
	REMOTE_ACCESS,

	/**
	 * This method is equivalent to Linux buffered I/O. It accesses data
	 * in a SAFS file with caching. It supports both synchronous I/O and
	 * asynchronous I/O.
	 */
	GLOBAL_CACHE_ACCESS,

	/**
	 * This method also accesses a SAFS file with caching. It localizes
	 * data access in the page cache. The implementation of this method
	 * is incomplete right now, so it shouldn't be used.
	 */
	PART_GLOBAL_ACCESS,
};

/**
 * This class defines the interface of accessing a SAFS file.
 * An I/O instance is not thread-safe, so we need to create an instance
 * for each thread. In the case of asynchronous I/O, the user-defined
 * callback function is guaranteed to be invoked in the same threads
 * as the I/O request was issued.
 */
class io_interface
{
	thread *curr;		// the thread where the IO instance runs.

	// This is an index for locating this IO object in a global table.
	int io_idx;
	int max_num_pending_ios;
	static atomic_integer io_counter;

protected:
	io_interface(thread *t) {
		this->curr = t;
		this->io_idx = io_counter.inc(1) - 1;
		max_num_pending_ios = params.get_max_num_pending_ios();
	}

public:
	typedef std::shared_ptr<io_interface> ptr;

	virtual ~io_interface() { }

	/**
	 * This method get the thread that the I/O instance is associated with.
	 * \return the thread.
	 */
	thread *get_thread() const {
		assert(curr);
		return curr;
	}

	/**
	 * This method gets the NUMA node that the I/O instance is used.
	 * \return the NUMA node id.
	 */
	int get_node_id() const {
		assert(curr);
		return curr->get_node_id();
	}

	/**
	 * This method gets the ID of the I/O instance.
	 * \return the ID of the I/O instance.
	 */
	int get_io_id() const {
		return io_idx;
	}

	/**
	 * This method gets the number of requests still allowed to be sent to
	 * the IO instance.
	 * \return the number of requests.
	 */
	int get_remaining_io_slots() const {
		return get_max_num_pending_ios() - num_pending_ios();
	}

	/**
	 * This method returns the ID of the file being accessed by the IO
	 * instance.
	 * \return the file ID.
	 */
	virtual int get_file_id() const = 0;

	/**
	 * This method is called when the I/O instance is destroyed.
	 */
	virtual void cleanup() {
	}

	/*
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
	 * This method indicates whether it supports asynchronous IO interface.
	 * \return boolean.
	 */
	virtual bool support_aio() {
		return false;
	}

	/**
	 * The asynchronous IO should implement some of the following methods.
	 */

	/**
	 * The main interface of sending asynchronous I/O requests.
	 * \param requests an array of I/O requests to be issued.
	 * \param num the number of I/O requests to be issued.
	 * \param status an array of I/O status, one for each I/O request.
	 */
	virtual void access(io_request *requests, int num, io_status *status = NULL) {
		throw unsupported_exception();
	}
	/**
	 * This method flushes I/O requests.
	 * When requests are passed to the access method, an IO layer may buffer
	 * the requests. This method guarantees that all requests are flushed to
	 * the underlying devices.
	 */
	virtual void flush_requests() {
		throw unsupported_exception();
	}
	/**
	 * This method waits for at least the specified number of requests
	 * issued by the access method to complete.
	 * \param the number of requests that need to be completed.
	 */
	virtual int wait4complete(int num) {
		throw unsupported_exception();
	}
	/**
	 * This method gets the number of I/O requests pending in the I/O instance.
	 * \return the number of pending I/O requests.
	 */
	virtual int num_pending_ios() const {
		throw unsupported_exception();
	}
	/**
	 * This method gets the maximal number of I/O requests allowed to
	 * be pending in the I/O instance.
	 * \return the maximal number of pending I/O requests.
	 */
	virtual int get_max_num_pending_ios() const {
		return max_num_pending_ios;
	}
	/**
	 * This method sets the maximal number of I/O requests allowed to be
	 * pending in the I/O instance.
	 * \param the maximal number of I/O requests.
	 */
	virtual void set_max_num_pending_ios(int max) {
		this->max_num_pending_ios = max;
	}

	/**
	 * \internal
	 * This method gives the underlying layer an interface to notify
	 * the current IO of the completed requests.
	 * The method may be called by multiple threads, so it has to be made
	 * thread-safe.
	 * The requests should be issued by this IO.
	 */
	virtual void notify_completion(io_request *reqs[], int num) {
		if (have_callback())
			get_callback().invoke(reqs, num);
	}

	/**
	 * This method sets the callback if the class supports the asynchronous
	 * I/O.
	 * \param cb the user-defined callback.
	 * \return false if the class doesn't support async I/O.
	 */
	virtual bool set_callback(callback::ptr cb) {
		throw unsupported_exception();
	}

	virtual bool have_callback() const {
		return false;
	}

	/**
	 * This method gets the user-defined callback.
	 * \return the user-defined callback.
	 */
	virtual callback &get_callback() {
		throw unsupported_exception();
	}

	/**
	 * The synchronous IO interface.
	 * \param buf the data buffer.
	 * \param off the location in the file where data should be read
	 * from or written to.
	 * \param size the size of data.
	 * \param access_method whether to read or write.
	 * \return I/O status for the request.
	 */
	virtual io_status access(char *buf, off_t off, ssize_t size,
			int access_method) {
		return IO_UNSUPPORTED;
	}
};

/**
 * This class creates an I/O scheduler used in the page cache.
 */
class comp_io_sched_creator
{
public:
	typedef std::shared_ptr<comp_io_sched_creator> ptr;
	/**
	 * This method creates a I/O scheduler on the specifies NUMA node.
	 * \param node_id the NUMA node ID.
	 * \return the I/O scheduler.
	 */
	virtual comp_io_scheduler::ptr create(int node_id) const = 0;
};

/**
 * This class defines the interface of creating I/O instances of accessing
 * a file.
 */
class file_io_factory
{
	comp_io_sched_creator::ptr creator;
	// The name of the file.
	const std::string name;
public:
	typedef std::shared_ptr<file_io_factory> shared_ptr;

	file_io_factory(const std::string _name): name(_name) {
	}

	virtual ~file_io_factory() {
	}

	/**
	 * This method sets an I/O scheduler creator in the I/O factory.
	 * An I/O scheduler only works in the page cache, so it has no effect
	 * if the I/O instance doesn't have a page cache.
	 * \param creater the I/O scheduler creator.
	 */
	void set_sched_creator(comp_io_sched_creator::ptr creator) {
		this->creator = creator;
	}

	/**
	 * This method gets an I/O scheduler creator in the I/O factory.
	 * \return the I/O scheduler creator.
	 */
	comp_io_sched_creator::ptr get_sched_creator() const {
		return creator;
	}

	/**
	 * This method gets the name of the SAFS file that the I/O instances
	 * in the I/O factory access.
	 * \return the file name.
	 */
	const std::string &get_name() const {
		return name;
	}

	/**
	 * This method gets the Id of the file associated with the I/O factory.
	 * \return the file ID.
	 */
	virtual int get_file_id() const = 0;

	/**
	 * This method creates an I/O instance for the specified thread.
	 * \param t the thread where the I/O instance can be used.
	 * \return the I/O instance.
	 */
	virtual io_interface::ptr create_io(thread *t) = 0;

	/*
	 * This method destroys an I/O instance created in the I/O factory.
	 * \param io the I/O instance to be destroyed.
	 */
	virtual void destroy_io(io_interface *io) = 0;
	virtual void print_state() {
	}

	virtual void collect_stat(io_interface &) {
	}

	virtual void print_statistics() const {
	}

	/**
	 * This method gets the size of the file accessed by the I/O factory.
	 * \return the file size.
	 */
	ssize_t get_file_size() const;
};

class cache_config;
class RAID_config;

/**
 * This function creates an I/O factory of the specified I/O method.
 * \param file_name the SAFS file accessed by the I/O factory.
 * \param access_option the I/O method of accessing the SAFS file.
 * The I/O method can be one of REMOTE_ACCESS, GLOBAL_CACHE_ACCESS and
 * PART_GLOBAL_ACCESS.
 */
file_io_factory::shared_ptr create_io_factory(const std::string &file_name,
		const int access_option);

/**
 * This function initializes SAFS. It should be called at the beginning
 * of a program. It can be invoked multiple times. If it is executed multiple
 * time successfully, destroy_io_system() needs to be invoked the same number
 * of times to complete clean up the I/O system.
 * \param map the SAFS configuration.
 * \param with_cache determine whether the I/O system is initialized with
 * page cache.
 */
void init_io_system(config_map::ptr map, bool with_cache = true);

/**
 * This function destroys SAFS. It should be called at the end of a program.
 */
void destroy_io_system();

/**
 * This function tells whether SAFS is initialized successfully.
 */
bool is_safs_init();

/**
 * This function gets the RAID configuration of SAFS.
 */
const RAID_config &get_sys_RAID_conf();

/**
 * \internal
 * This method print the I/O statistic information. It's used for debugging.
 */
void print_io_thread_stat();

/**
 * The users can set the weight of a file. The file weight is used by
 * the page cache. The file with a higher weight can have its data in
 * the page cache longer.
 * This function isn't thread-safe and should be used before I/O instances
 * are created.
 * \param file_name The file name.
 * \param weight The new weight of the file.
 */
void set_file_weight(const std::string &file_name, int weight);

}

#endif
