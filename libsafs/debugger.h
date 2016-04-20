#ifndef __DEBUGGER_H__
#define __DEBUGGER_H__

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

#include <map>

#include "concurrency.h"

namespace safs
{

class debug_task
{
public:
	virtual ~debug_task() {
	}

	virtual void run() = 0;
};

class debugger
{
	// We have to use map instead of unordered_map here because
	// we want to preserve the order of the tasks.
	std::map<int, debug_task *> tasks;
	pthread_spinlock_t lock;
	atomic_integer task_id_gen;
public:
	debugger();

	~debugger();

	int register_task(debug_task *task) {
		int id = task_id_gen.inc(1);
		pthread_spin_lock(&lock);
		tasks.insert(std::pair<int, debug_task *>(id, task));
		pthread_spin_unlock(&lock);
		return id;
	}

	void remove_task(int task_id) {
		pthread_spin_lock(&lock);
		tasks.erase(task_id);
		pthread_spin_unlock(&lock);
	}

	void run();
};

bool is_debug_enabled();

}

#endif
