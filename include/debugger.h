#ifndef __DEBUGGER_H__
#define __DEBUGGER_H__

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

#include <pthread.h>

#include <map>

#include "concurrency.h"

class debug_task
{
public:
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

extern debugger debug;

bool is_debug_enabled();

#endif
