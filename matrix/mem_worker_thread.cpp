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

#include <unordered_map>

#include "io_interface.h"

#include "matrix_config.h"
#include "mem_worker_thread.h"
#include "EM_object.h"

namespace fm
{

namespace detail
{

mem_thread_pool::mem_thread_pool(int num_nodes, int nthreads_per_node)
{
	tot_num_tasks = 0;
	threads.resize(num_nodes);
	for (int i = 0; i < num_nodes; i++) {
		threads[i].resize(nthreads_per_node);
		for (int j = 0; j < nthreads_per_node; j++) {
			std::string name
				= std::string("mem-worker-") + itoa(i) + "-" + itoa(j);
			threads[i][j] = std::shared_ptr<pool_task_thread>(
					new pool_task_thread(i * nthreads_per_node + j, name, i));
			threads[i][j]->start();
		}
	}
	ntasks_per_node.resize(num_nodes);
}

/*
 * This method dispatches tasks in a round-robin fashion.
 */
void mem_thread_pool::process_task(int node_id, thread_task *task)
{
	if (node_id < 0)
		node_id = tot_num_tasks % get_num_nodes();

	assert((size_t) node_id < threads.size());
	size_t idx = ntasks_per_node[node_id] % threads[node_id].size();
	threads[node_id][idx]->add_task(task);
	ntasks_per_node[node_id]++;
	tot_num_tasks++;
}

void mem_thread_pool::wait4complete()
{
	for (size_t i = 0; i < threads.size(); i++) {
		size_t nthreads = threads[i].size();
		for (size_t j = 0; j < nthreads; j++)
			threads[i][j]->wait4complete();
	}
}

size_t mem_thread_pool::get_num_pending() const
{
	size_t ret = 0;
	for (size_t i = 0; i < threads.size(); i++) {
		size_t nthreads = threads[i].size();
		for (size_t j = 0; j < nthreads; j++)
			ret += threads[i][j]->get_num_pending();
	}
	return ret;
}

static mem_thread_pool::ptr global_threads;

mem_thread_pool::ptr mem_thread_pool::get_global_mem_threads()
{
	if (global_threads == NULL) {
		int nthreads_per_node
			= matrix_conf.get_num_threads() / matrix_conf.get_num_nodes();
		assert(nthreads_per_node > 0);
		global_threads = mem_thread_pool::create(matrix_conf.get_num_nodes(),
				nthreads_per_node);
	}
	return global_threads;
}

void mem_thread_pool::init_global_mem_threads(int num_nodes,
		int nthreads_per_node)
{
	if (global_threads == NULL)
		global_threads = mem_thread_pool::create(num_nodes, nthreads_per_node);
}

void io_worker_task::run()
{
	std::vector<safs::io_interface::ptr> ios;
	pthread_spin_lock(&lock);
	for (auto it = EM_objs.begin(); it != EM_objs.end(); it++) {
		std::vector<safs::io_interface::ptr> tmp = (*it)->create_ios();
		ios.insert(ios.end(), tmp.begin(), tmp.end());
	}
	pthread_spin_unlock(&lock);

	// The task runs until there are no tasks left in the queue.
	while (dispatch->issue_task()) {
		// TODO we need to make this a parameter.
		for (size_t i = 0; i < ios.size(); i++) {
			ios[i]->wait4complete(0);
			while (ios[i]->num_pending_ios() > max_pending_ios)
				ios[i]->wait4complete(1);
		}
	}
	// Test if all I/O instances have processed all requests.
	bool complete;
	do {
		complete = true;
		for (size_t i = 0; i < ios.size(); i++) {
			ios[i]->wait4complete(0);
			// If there is still an I/O instance has pending requests,
			// we need to start over and test all I/O instances again.
			if (ios[i]->num_pending_ios() > 0) {
				complete = false;
				ios[i]->wait4complete(1);
			}
		}
		// If all I/O instances have no pending I/O requests left.
	} while (!complete);

	pthread_spin_lock(&lock);
	EM_objs.clear();
	pthread_spin_unlock(&lock);
}

}

}
