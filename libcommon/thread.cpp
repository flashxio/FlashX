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

#include <stdio.h>
#include <stdlib.h>
#include <numa.h>
#include <pthread.h>

#include "thread.h"
#include "common.h"

static void bind2node_id(int node_id)
{
	struct bitmask *bmp = numa_allocate_nodemask();
	numa_bitmask_setbit(bmp, node_id);
	numa_bind(bmp);
	numa_free_nodemask(bmp);
}

#if 0
static int get_numa_run_node()
{
	struct bitmask *bmp = numa_get_run_node_mask();
	int nbytes = numa_bitmask_nbytes(bmp);
	int num_nodes = 0;
	int node_id = -1;
	int i;
	for (i = 0; i < nbytes * 8; i++)
		if (numa_bitmask_isbitset(bmp, i)) {
			num_nodes++;
			printf("bind to node %d\n", i);
			node_id = i;
		}
	return node_id;
}

static int numa_get_mem_node()
{
	struct bitmask *bmp = numa_get_membind();
	int nbytes = numa_bitmask_nbytes(bmp);
	int num_nodes = 0;
	int node_id = -1;
	int i;
	for (i = 0; i < nbytes * 8; i++)
		if (numa_bitmask_isbitset(bmp, i)) {
			num_nodes++;
			node_id = i;
		}
	assert(num_nodes == 1);
	return node_id;
}
#endif

void *thread_run(void *arg)
{
	thread *t = (thread *) arg;
	pthread_setspecific(thread::thread_key, t);
	t->tid = gettid();

	int node_id = t->get_node_id();
	if (node_id >= 0) {
		bind2node_id(node_id);
		numa_set_bind_policy(1);
	}
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
