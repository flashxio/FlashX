#ifndef __MY_TIMER_H__
#define __MY_TIMER_H__

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

#include <sys/time.h>

#include <map>

#include "concurrency.h"

static inline int64_t get_curr_time_us()
{
	struct timeval curr;
	gettimeofday(&curr, NULL);
	int64_t curr_us = curr.tv_sec;
	curr_us = curr_us * 1000000 + curr.tv_usec;
	return curr_us;
}

class timer_task
{
	static atomic_number<long> task_count;
	long task_id;
	int64_t abs_timeout;
	long timeout;
public:
	// timeout: in microsecond
	timer_task(long timeout) {
		task_id = task_count.inc(1);
		int64_t curr = get_curr_time_us();
		abs_timeout = curr + timeout;
		this->timeout = timeout;
	}

	virtual ~timer_task() {
	}

	long get_id() const {
		return task_id;
	}

	int64_t get_timeout() const {
		return timeout;
	}

	int64_t get_abs_timeout() const {
		return abs_timeout;
	}

	void renew() {
		int64_t curr = get_curr_time_us();
		abs_timeout = curr + timeout;
	}

	virtual void run() = 0;
};

class thread;

class timer
{
	static atomic_integer timer_count;

	// I treat this as a sorted list.
	// key is the absolute timeout.
	std::multimap<int64_t, timer_task *> task_list;
	spin_lock lock;
	timer_t timerid;
	thread *t;
	int timer_id;

	void set_timeout(int64_t timeout);
	void init_timer();
	int64_t _add_task(timer_task *task);
public:
	timer();

	timer(thread *t);

	bool add_task(timer_task *task);

	void run_tasks();

	void print_state();

	int64_t get_next_timeout();

	int get_id() const {
		return timer_id;
	}
};

class periodic_timer
{
	static atomic_integer timer_count;
	int timer_id;
	timer_t timerid;
	class timer_task *task;
	thread *t;

	void set_timeout(int64_t timeout);
public:
	periodic_timer(thread *t, timer_task *task);
	void run_task();
	int get_id() const {
		return timer_id;
	}

	bool is_enabled() const;
	bool start();
};

#endif

