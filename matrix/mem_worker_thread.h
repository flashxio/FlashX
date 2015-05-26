#ifndef __MEM_WORKER_THREAD_H__
#define __MEM_WORKER_THREAD_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
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

#include <memory>

#include "thread.h"

#include "local_vec_store.h"

namespace safs
{
class file_io_factory;
class io_interface;
class io_request;
}

namespace fm
{

namespace detail
{

class pool_task_thread: public task_thread
{
	int pool_thread_id;
public:
	pool_task_thread(int pool_thread_id, const std::string &name,
			int node): task_thread(name, node) {
		this->pool_thread_id = pool_thread_id;
	}

	int get_pool_thread_id() const {
		return pool_thread_id;
	}
};

/*
 * This is designed to replace openmp for parallelization while respecting
 * NUMA data locality.
 */
class mem_thread_pool
{
	std::vector<size_t> ntasks_per_node;
	std::vector<std::vector<std::shared_ptr<pool_task_thread> > > threads;

	mem_thread_pool(int num_nodes, int nthreads_per_node);
public:
	typedef std::shared_ptr<mem_thread_pool> ptr;

	static ptr get_global_mem_threads();
	static void init_global_mem_threads(int num_nodes, int nthreads_per_node);

	static ptr create(int num_nodes, int nthreads_per_node) {
		return ptr(new mem_thread_pool(num_nodes, nthreads_per_node));
	}

	size_t get_num_nodes() const {
		return ntasks_per_node.size();
	}
	size_t get_num_threads() const {
		assert(!threads.empty());
		return threads.size() * threads.front().size();
	}

	void process_task(int node_id, thread_task *task);

	void wait4complete();
};

class io_worker_task;

/*
 * This runs on the portion of the data in a data container when the portion
 * of data is available in memory.
 */
class portion_compute
{
	// This object doesn't own the worker task. We have to make sure
	// the worker task is alive when this object is alive.
	io_worker_task *worker_task;
public:
	typedef std::shared_ptr<portion_compute> ptr;

	virtual ~portion_compute() {
	}

	io_worker_task &get_worker_task() {
		return *worker_task;
	}

	void set_worker(io_worker_task *task) {
		worker_task = task;
	}

	bool register_portion_compute(const safs::io_request &req,
			portion_compute::ptr compute);

	virtual void run(char *buf, size_t size) = 0;
};

/*
 * When we write data to disks, we need to have something to hold the buffer.
 * This holds the local buffer until the write completes.
 */
class portion_write_complete: public portion_compute
{
	local_buf_vec_store::ptr store;
public:
	portion_write_complete(local_buf_vec_store::ptr store) {
		this->store = store;
	}

	virtual void run(char *buf, size_t size) {
		assert(store->get_raw_arr() == buf);
		assert(store->get_length() * store->get_entry_size() == size);
	}
};

/*
 * This is an I/O task that defines computation on the portion of data
 * in an external-memory data container.
 */
class portion_io_task
{
	// This object doesn't own the worker task. We have to make sure
	// the worker task is alive when this object is alive.
	io_worker_task *worker_task;
public:
	typedef std::shared_ptr<portion_io_task> ptr;

	virtual ~portion_io_task() {
	}

	io_worker_task &get_worker_task() {
		return *worker_task;
	}

	void set_worker(io_worker_task *task) {
		worker_task = task;
	}

	bool register_portion_compute(const safs::io_request &req,
			portion_compute::ptr compute);

	virtual void run(std::shared_ptr<safs::io_interface> read_io,
			std::shared_ptr<safs::io_interface> write_io) = 0;
};

/*
 * This defines a set of I/O tasks that process an entire data container.
 */
class task_dispatcher
{
public:
	typedef std::shared_ptr<task_dispatcher> ptr;

	virtual ~task_dispatcher() {
	}
	/*
	 * Return a task to the invoker.
	 * This method must be thread-safe.
	 */
	virtual portion_io_task::ptr get_task() = 0;
};

class portion_callback;

class io_worker_task: public thread_task
{
	std::shared_ptr<safs::file_io_factory> read_factory;
	std::shared_ptr<safs::file_io_factory> write_factory;
	std::shared_ptr<portion_callback> cb;
	task_dispatcher::ptr dispatch;
	int max_pending_ios;
public:
	/*
	 * If the read and write I/O factories are the same, only one I/O instance
	 * is created and is used for both read and write. If only the read
	 * I/O factory is created, only a read I/O instance is created.
	 * If only the write I/O factory is created, only a write I/O instance
	 * is created.
	 */
	io_worker_task(std::shared_ptr<safs::file_io_factory> read_factory,
			std::shared_ptr<safs::file_io_factory> write_factory,
			task_dispatcher::ptr dispatch, int max_pending_ios = 16) {
		this->read_factory = read_factory;
		this->write_factory = write_factory;
		this->dispatch = dispatch;
		this->max_pending_ios = max_pending_ios;
	}

	bool register_portion_compute(const safs::io_request &req,
			portion_compute::ptr compute);

	virtual void run();
};

}

}

#endif
