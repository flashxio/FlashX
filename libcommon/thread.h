#ifndef __MY_THREAD_H__
#define __MY_THREAD_H__

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

#include <pthread.h>
#ifdef USE_HWLOC
#include <hwloc.h>
#endif

#include <atomic>
#include <string>
#include <set>

#include "concurrency.h"
#include "common.h"
#include "container.h"

class thread
{
	static pthread_key_t thread_key;
	static atomic_integer num_threads;

	volatile pid_t tid;
	int thread_idx;
	int node_id;
	// Indicate which CPU cores the thread is bound to.
	std::vector<int> cpu_affinity;
	pthread_t id;
	bool blocking;
	std::string name;

	volatile void *user_data;

	volatile bool _is_sleeping;
	volatile bool _is_running;
	volatile bool _has_exit;
	volatile bool _is_activated;

	pthread_mutex_t mutex;
	pthread_cond_t cond;

	void construct_init();

	friend void init_thread_class();
	friend void *thread_run(void *arg);
public:
	thread(std::string name, int node_id, bool blocking = true);
	thread(std::string name, const std::vector<int> &cpu_affinity,
			bool blocking = true);

	void set_user_data(void *user_data) {
		assert(this->user_data == NULL);
		this->user_data = user_data;
	}

	void *get_user_data() const {
		return (void *) user_data;
	}

	virtual ~thread() {
		stop();
		if (thread_idx >= 0)
			join();
	}

	void activate() {
		bool sleeping;
		pthread_mutex_lock(&mutex);
		_is_activated = true;
		sleeping = _is_sleeping;
		pthread_mutex_unlock(&mutex);
		if (sleeping)
			pthread_cond_signal(&cond);
	}

	void wait() {
		pthread_mutex_lock(&mutex);
		if (!_is_activated) {
			while (!_is_activated && _is_running) {
				_is_sleeping = true;
				int ret = pthread_cond_wait(&cond, &mutex);
				_is_sleeping = false;
				if (ret)
					perror("pthread_cond_wait");
			}
		}
		_is_activated = false;
		pthread_mutex_unlock(&mutex);
	}

	bool is_running() const {
		return _is_running;
	}

	bool is_blocking() const {
		return blocking;
	}

	int get_node_id() const {
		return node_id;
	}

	const std::vector<int> get_cpu_affinity() const {
		return cpu_affinity;
	}

	void stop() {
		pthread_mutex_lock(&mutex);
		// This is the only place where `is_running' is changed.
		_is_running = false;
		pthread_cond_signal(&cond);
		pthread_mutex_unlock(&mutex);
	}

	void exit() {
		_has_exit = true;
	}

	bool has_exit() const {
		return _has_exit;
	}

	void join() {
		pthread_join(id, NULL);
#ifdef DEBUG
		printf("stop thread %s\n", name.c_str());
#endif
		thread_idx = -1;
	}

	int get_id() const {
		return thread_idx;
	}
	int get_tid() const {
		if (get_curr_thread() == NULL || get_curr_thread() != this) {
			// TODO I need a better way to wait for tid to be initialized.
			while (tid < 0) { }
			return tid;
		}
		else
			return gettid();
	}

	void start();
	virtual void run() = 0;
	virtual void init() {
	}
	virtual void cleanup() {
	}

	/*
	 * This is to initialize the thread class instead of a single thread.
	 */
	static void thread_class_init();

	static thread *get_curr_thread();

	/*
	 * This creates a thread instance to represent the current thread context.
	 * It is used when the current thread isn't created by the thread class.
	 */
	static thread *represent_thread(int node_id);

	const std::string &get_thread_name() const {
		return name;
	}
};

class thread_task
{
public:
	virtual ~thread_task() {
	}
	virtual void run() = 0;
};

class task_thread: public thread
{
	fifo_queue<thread_task *> tasks;
	std::atomic<size_t> num_pending;
	bool all_complete;
	pthread_mutex_t mutex;
	pthread_cond_t cond;
public:
	task_thread(const std::string &name, int node): thread(name,
			node), tasks(node, 1024, true) {
		pthread_mutex_init(&mutex, NULL);
		pthread_cond_init(&cond, NULL);
		all_complete = false;
		num_pending = 0;
	}

	task_thread(const std::string &name, const std::vector<int> &cpus,
			int node): thread(name, cpus), tasks(node, 1024, true) {
		pthread_mutex_init(&mutex, NULL);
		pthread_cond_init(&cond, NULL);
		all_complete = false;
		num_pending = 0;
	}

	void add_task(thread_task *t);
	void run();
	void wait4complete();

	size_t get_num_pending() const {
		return num_pending;
	}
};

#ifdef USE_HWLOC

class CPU_core
{
	std::vector<int> logical_units;
public:
	CPU_core(hwloc_obj_t core);

	const std::vector<int> get_units() const {
		return logical_units;
	}

	off_t get_logical_unit(size_t idx) const {
		return logical_units[idx];
	}

	size_t get_num_units() const {
		return logical_units.size();
	}
};

class NUMA_node
{
	std::vector<CPU_core> cores;
	std::set<int> lus;
public:
	/* This constructor works for the machine without NUMA nodes. */
	NUMA_node(hwloc_topology_t topology);
	/* This constructor works for the machine with NUMA nodes. */
	NUMA_node(hwloc_obj_t node);

	bool contain_lu(int unit) const {
		return lus.find(unit) != lus.end();
	}

	std::vector<int> get_logical_units() const;

	const CPU_core &get_core(size_t idx) const {
		return cores[idx];
	}

	size_t get_num_cores() const {
		return cores.size();
	}

	size_t get_num_logical_units() const {
		return cores.size() * cores.front().get_num_units();
	}
};

class CPU_hierarchy
{
	std::vector<NUMA_node> nodes;
public:
	CPU_hierarchy();

	const NUMA_node &get_node(size_t idx) const {
		return nodes[idx];
	}

	size_t get_num_nodes() const {
		return nodes.size();
	}

	size_t get_num_cores() const {
		return nodes.size() * nodes.front().get_num_cores();
	}

	size_t get_num_logical_units() const {
		return nodes.size() * nodes.front().get_num_logical_units();
	}

	std::vector<int> lus2node(const std::vector<int> &lus) const;
};

extern CPU_hierarchy cpus;

#endif

#endif
