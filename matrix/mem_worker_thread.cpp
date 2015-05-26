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

namespace fm
{

namespace detail
{

mem_thread_pool::mem_thread_pool(int num_nodes, int nthreads_per_node)
{
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
	assert((size_t) node_id < threads.size());
	size_t idx = ntasks_per_node[node_id] % threads[node_id].size();
	threads[node_id][idx]->add_task(task);
	ntasks_per_node[node_id]++;
}

void mem_thread_pool::wait4complete()
{
	for (size_t i = 0; i < threads.size(); i++) {
		size_t nthreads = threads[i].size();
		for (size_t j = 0; j < nthreads; j++)
			threads[i][j]->wait4complete();
	}
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
	global_threads = mem_thread_pool::create(num_nodes, nthreads_per_node);
}

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

bool io_worker_task::register_portion_compute(const safs::io_request &req,
		portion_compute::ptr compute)
{
	if (cb) {
		cb->add(req, compute);
		return true;
	}
	else
		return false;
}

void io_worker_task::run()
{
	safs::io_interface::ptr read_io, write_io;
	if (read_factory == write_factory)
		read_io = write_io = create_io(read_factory, thread::get_curr_thread());
	else if (read_factory == NULL)
		write_io = create_io(write_factory, thread::get_curr_thread());
	else if (write_factory == NULL)
		read_io = create_io(read_factory, thread::get_curr_thread());
	else {
		write_io = create_io(write_factory, thread::get_curr_thread());
		read_io = create_io(read_factory, thread::get_curr_thread());
	}
	cb = portion_callback::ptr(new portion_callback());
	if (read_io)
		read_io->set_callback(cb);
	if (write_io)
		write_io->set_callback(cb);

	portion_io_task::ptr task;
	// The task runs until there are no tasks left in the queue.
	while ((task = dispatch->get_task()) != NULL) {
		task->set_worker(this);
		task->run(read_io, write_io);
		// TODO we need to make this a parameter.
		while (read_io && read_io->num_pending_ios() > max_pending_ios)
			read_io->wait4complete(1);
		while (write_io && write_io->num_pending_ios() > max_pending_ios)
			write_io->wait4complete(1);
	}
	if (read_io)
		while (read_io->num_pending_ios() > 0)
			read_io->wait4complete(read_io->num_pending_ios());
	if (write_io)
		while (write_io->num_pending_ios() > 0)
			write_io->wait4complete(write_io->num_pending_ios());
}

bool portion_compute::register_portion_compute(const safs::io_request &req,
		portion_compute::ptr compute)
{
	compute->set_worker(worker_task);
	return worker_task->register_portion_compute(req, compute);
}

bool portion_io_task::register_portion_compute(const safs::io_request &req,
		portion_compute::ptr compute)
{
	compute->set_worker(worker_task);
	return worker_task->register_portion_compute(req, compute);
}

}

}
