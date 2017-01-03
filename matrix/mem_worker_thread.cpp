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
#include "local_mem_buffer.h"

namespace fm
{

namespace detail
{

#ifdef USE_HWLOC
std::vector<int> get_cpus(int node_id)
{
	std::vector<int> io_cpus = safs::get_io_cpus();
	std::vector<int> logical_units = cpus.get_node(node_id).get_logical_units();
	std::set<int> cpu_set(logical_units.begin(), logical_units.end());
	// Remove the logical units where I/O threads run.
	for (size_t j = 0; j < io_cpus.size(); j++)
		cpu_set.erase(io_cpus[j]);
	return std::vector<int>(cpu_set.begin(), cpu_set.end());
}
#endif

mem_thread_pool::mem_thread_pool(int num_nodes, int nthreads_per_node)
{
	tot_num_tasks = 0;
	threads.resize(num_nodes);
	for (int i = 0; i < num_nodes; i++) {
		// Get the CPU cores that are in node i.
#ifdef USE_HWLOC
		std::vector<int> cpus = get_cpus(i);
#endif
		threads[i].resize(nthreads_per_node);
		for (int j = 0; j < nthreads_per_node; j++) {
			std::string name
				= std::string("mem-worker-") + itoa(i) + "-" + itoa(j);
#ifdef USE_HWLOC
			if (safs::get_io_cpus().empty())
				threads[i][j] = std::shared_ptr<pool_task_thread>(
						new pool_task_thread(i * nthreads_per_node + j, name, i));
			else
				threads[i][j] = std::shared_ptr<pool_task_thread>(
						new pool_task_thread(i * nthreads_per_node + j, name,
							cpus, i));
#else
			threads[i][j] = std::shared_ptr<pool_task_thread>(
					new pool_task_thread(i * nthreads_per_node + j, name, i));
#endif
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

	// After all workers complete, we should try to clear the memory buffers
	// in each worker thread to reduce memory consumption.
	// We might want to keep the memory buffer for I/O on dense matrices.
	if (matrix_conf.is_keep_mem_buf())
		detail::local_mem_buffer::clear_bufs(
				detail::local_mem_buffer::MAT_PORTION);
	else
		detail::local_mem_buffer::clear_bufs();
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
			= matrix_conf.get_num_DM_threads() / matrix_conf.get_num_nodes();
		assert(nthreads_per_node > 0);
		global_threads = mem_thread_pool::create(matrix_conf.get_num_nodes(),
				nthreads_per_node);
	}
	return global_threads;
}

size_t mem_thread_pool::get_global_num_threads()
{
	// We also count the main thread.
	return get_global_mem_threads()->get_num_threads() + 1;
}

int mem_thread_pool::get_curr_thread_id()
{
	// It return 0 for the main thread. The worker thread Id starts with 1.
	detail::pool_task_thread *curr
		= dynamic_cast<detail::pool_task_thread *>(thread::get_curr_thread());
	if (curr)
		return curr->get_pool_thread_id() + 1;
	// If this isn't a pool thread, it must be the main thread.
	else
		return 0;
}

void mem_thread_pool::init_global_mem_threads(int num_nodes,
		int nthreads_per_node)
{
	if (global_threads == NULL)
		global_threads = mem_thread_pool::create(num_nodes, nthreads_per_node);
}

void mem_thread_pool::destroy()
{
	global_threads = NULL;
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
	safs::io_select::ptr select = safs::create_io_select(ios);

	// The task runs until there are no tasks left in the queue.
	while (dispatch->issue_task())
		safs::wait4ios(select, max_pending_ios);
	// Test if all I/O instances have processed all requests.
	size_t num_pending = safs::wait4ios(select, 0);
	assert(num_pending == 0);

	pthread_spin_lock(&lock);
	EM_objs.clear();
	pthread_spin_unlock(&lock);

	for (size_t i = 0; i < ios.size(); i++) {
		portion_callback &cb = static_cast<portion_callback &>(
				ios[i]->get_callback());
		assert(!cb.has_callback());
	}
}

global_counter::global_counter()
{
	counts.resize(mem_thread_pool::get_global_num_threads());
	reset();
}

void global_counter::inc(size_t val)
{
	int id = mem_thread_pool::get_curr_thread_id();
	counts[id].count += val;
}

void global_counter::reset()
{
	for (size_t i = 0; i < counts.size(); i++)
		counts[i].count = 0;
}

size_t global_counter::get() const
{
	size_t tot = 0;
	for (size_t i = 0; i < counts.size(); i++)
		tot += counts[i].count;
	return tot;
}

}

}
