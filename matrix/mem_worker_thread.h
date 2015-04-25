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

}

}

#endif
