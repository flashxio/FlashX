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
#include <set>
#include <unordered_map>

#include "thread.h"
#include "io_interface.h"

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
	size_t tot_num_tasks;
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

	size_t get_num_pending() const;

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
public:
	typedef std::shared_ptr<portion_compute> ptr;

	virtual ~portion_compute() {
	}

	virtual void run(char *buf, size_t size) = 0;
};

class portion_callback: public safs::callback
{
	std::unordered_map<char *, portion_compute::ptr> computes;
public:
	typedef std::shared_ptr<portion_callback> ptr;

	virtual ~portion_callback() {
		assert(computes.empty());
	}

	void add(const safs::io_request &req, portion_compute::ptr compute) {
		auto ret = computes.insert(std::pair<char *, portion_compute::ptr>(
					req.get_buf(), compute));
		assert(ret.second);
	}

	virtual int invoke(safs::io_request *reqs[], int num) {
		for (int i = 0; i < num; i++) {
			auto it = computes.find(reqs[i]->get_buf());
			assert(it != computes.end());
			portion_compute::ptr compute = it->second;
			computes.erase(it);
			compute->run(reqs[i]->get_buf(), reqs[i]->get_size());
		}
		return 0;
	}
};

/*
 * When we write data to disks, we need to have something to hold the buffer.
 * This holds the local buffer until the write completes.
 */
class portion_write_complete: public portion_compute
{
	local_buf_vec_store::const_ptr store;
public:
	portion_write_complete(local_buf_vec_store::const_ptr store) {
		this->store = store;
	}

	virtual void run(char *buf, size_t size) {
		assert(store->get_raw_arr() == buf);
		assert(store->get_length() * store->get_entry_size() == size);
	}
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
	 * Issue a task.
	 * This method must be thread-safe.
	 */
	virtual bool issue_task() = 0;
};

class EM_object;

class io_worker_task: public thread_task
{
	pthread_spinlock_t lock;
	std::set<EM_object *> EM_objs;

	std::shared_ptr<portion_callback> cb;
	task_dispatcher::ptr dispatch;
	int max_pending_ios;
public:
	io_worker_task(task_dispatcher::ptr dispatch, int max_pending_ios = 16) {
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		this->dispatch = dispatch;
		this->max_pending_ios = max_pending_ios;
	}

	~io_worker_task() {
		pthread_spin_destroy(&lock);
	}

	void register_EM_obj(EM_object *obj) {
		pthread_spin_lock(&lock);
		auto it = EM_objs.find(obj);
		assert(it == EM_objs.end());
		EM_objs.insert(obj);
		pthread_spin_unlock(&lock);
	}

	virtual void run();
};

}

}

#endif
