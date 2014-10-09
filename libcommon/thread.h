#ifndef __MY_THREAD_H__
#define __MY_THREAD_H__

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

#include <pthread.h>

#include <string>

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

	friend void init_thread_class();
	friend void *thread_run(void *arg);
public:
	thread(std::string name, int node_id, bool blocking = true) {
		thread_class_init();

		tid = -1;
		thread_idx = num_threads.inc(1);
		this->name = name + "-" + itoa(thread_idx);
		this->node_id = node_id;
		this->blocking = blocking;
		this->id = 0;

		_is_activated = false;
		_has_exit = false;
		_is_running = true;
		_is_sleeping = true;

		pthread_mutex_init(&mutex, NULL);
		pthread_cond_init(&cond, NULL);
		user_data = NULL;
	}

	void set_user_data(void *user_data) {
		assert(this->user_data == NULL);
		this->user_data = user_data;
	}

	void *get_user_data() const {
		return (void *) user_data;
	}

	virtual ~thread() {
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

	/**
	 * This is to initialize the thread class instead of a single thread.
	 */
	static void thread_class_init();

	static thread *get_curr_thread();

	/**
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
	bool all_complete;
	pthread_mutex_t mutex;
	pthread_cond_t cond;
public:
	task_thread(const std::string &name, int node): thread(name,
			node), tasks(node, 1024, true) {
		pthread_mutex_init(&mutex, NULL);
		pthread_cond_init(&cond, NULL);
		all_complete = false;
	}

	void add_task(thread_task *t) {
		pthread_mutex_lock(&mutex);
		all_complete = false;
		if (tasks.is_full())
			tasks.expand_queue(tasks.get_size() * 2);
		tasks.push_back(t);
		pthread_mutex_unlock(&mutex);
		activate();
	}

	void run() {
		const int TASK_BUF_SIZE = 128;
		thread_task *local_tasks[TASK_BUF_SIZE];
		pthread_mutex_lock(&mutex);
		while (!tasks.is_empty()) {
			int num_tasks = tasks.fetch(local_tasks, TASK_BUF_SIZE);
			pthread_mutex_unlock(&mutex);
			for (int i = 0; i < num_tasks; i++) {
				local_tasks[i]->run();
				delete local_tasks[i];
			}
			pthread_mutex_lock(&mutex);
		}
		all_complete = true;
		pthread_mutex_unlock(&mutex);
		pthread_cond_signal(&cond);
	}

	void wait4complete() {
		pthread_mutex_lock(&mutex);
		while (!all_complete) {
			pthread_cond_wait(&cond, &mutex);
		}
		pthread_mutex_unlock(&mutex);
	}
};

#endif
