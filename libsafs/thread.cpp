/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <numa.h>
#include <pthread.h>

#include "thread.h"
#include "common.h"

void *thread_run(void *arg)
{
	thread *t = (thread *) arg;
	pthread_setspecific(thread::thread_key, t);
	t->tid = gettid();

	int node_id = t->get_node_id();
	if (node_id >= 0)
		bind2node_id(node_id);
	t->init();
	while (t->is_running()) {
		t->run();
		if (t->is_blocking())
			t->wait();
	}
	t->cleanup();
	t->exit();
	return NULL;
}

void thread::start()
{
	assert(id == 0);
	int ret = pthread_create(&id, NULL, thread_run, (void *) this);
	if (ret) {
		perror("pthread_create");
		::exit(1);
	}
}

static pthread_once_t once_control = PTHREAD_ONCE_INIT;

void init_thread_class()
{
	pthread_key_create(&thread::thread_key, NULL);
}

void thread::thread_class_init()
{
	pthread_once(&once_control, init_thread_class);
}

thread *thread::get_curr_thread()
{
	thread *curr = (thread *) pthread_getspecific(thread_key);
	return curr;
}

pthread_key_t thread::thread_key;

class thread_representer: public thread
{
	void run() {
	}
public:
	thread_representer(int node_id): thread(std::string(
				"representer-thread-node") + itoa(node_id), node_id) {
	}
};

thread *thread::represent_thread(int node_id)
{
	if (node_id >= 0)
		bind2node_id(node_id);

	thread *curr = thread::get_curr_thread();
	assert(curr == NULL);
	curr = new thread_representer(node_id);
	pthread_setspecific(thread::thread_key, curr);
	curr->_is_sleeping = false;
	curr->tid = gettid();
	return curr;
}

atomic_integer thread::num_threads;
