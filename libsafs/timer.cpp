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

#include <signal.h>
#include <time.h>

#include "timer.h"
#include "thread.h"

static void handler(int sig, siginfo_t *si, void *uc)
{
	timer *t = (timer *) si->si_value.sival_ptr;
	t->run_tasks();
}

static void disable_timer()
{
	sigset_t mask;

	sigemptyset(&mask);
	sigaddset(&mask, SIGRTMIN);
	if (sigprocmask(SIG_SETMASK, &mask, NULL) == -1) {
		perror("sigprocmask");
		exit(1);
	}
}

static void enable_timer()
{
	sigset_t mask;

	sigemptyset(&mask);
	sigaddset(&mask, SIGRTMIN);
	if (sigprocmask(SIG_UNBLOCK, &mask, NULL) == -1) {
		perror("sigprocmask");
		exit(1);
	}
}

class timer_thread: public thread
{
public:
	timer_thread(): thread("timer_thread", 0) {
	}

	void run() {
	}
};

void timer::init_timer()
{
	/**
	 * The code here is copied from the example code in the manual of
	 * timer_create.
	 */
	struct sigevent sev;
	struct sigaction sa;

	/* Establish handler for timer signal */

	sa.sa_flags = SA_SIGINFO;
	sa.sa_sigaction = handler;
	sigemptyset(&sa.sa_mask);
	if (sigaction(SIGRTMIN, &sa, NULL) == -1) {
		perror("sigaction");
		exit(1);
	}

	// This is a Linux feature, but we really want the signal sent to a dedicated
	// thread, so whatever the timer task does, there won't be deadlock.
	sev.sigev_notify = SIGEV_THREAD_ID;
	// It supposes to be sigev_notify_thread_id
	sev._sigev_un._tid = t->get_tid();
	sev.sigev_signo = SIGRTMIN;
	sev.sigev_value.sival_ptr = this;
	if (timer_create(CLOCK_REALTIME, &sev, &timerid) == -1) {
		perror("timer_create");
		exit(1);
	}
}

timer::timer(thread *t)
{
	timer_id = timer_count.inc(1);
	this->t = t;

	init_timer();
}

timer::timer()
{
	timer_id = timer_count.inc(1);
	t = new timer_thread();
	t->start();

	init_timer();
}

/**
 * Add the new timer task to the task list.
 * If the new task is the first to timeout, return timeout in microseconds,
 * otherwise, return -1.
 */
int64_t timer::_add_task(timer_task *task)
{
	// relative timeout in microseconds.
	int64_t timeout = -1;
	lock.lock();
	if (!task_list.empty()) {
		// We want to have the new task to be after the first timer task.
		// It can save us some trouble.
		assert(task_list.begin()->second->get_abs_timeout()
				<= task->get_abs_timeout());
		// If the new task is the first task in the list.
		if (task_list.begin()->second == task)
			timeout = task->get_abs_timeout() - get_curr_time_us();
	}
	else {
		timeout = task->get_abs_timeout() - get_curr_time_us();
		assert(timeout >= 0);
	}
	task_list.insert(std::pair<int64_t, timer_task *>(
				task->get_abs_timeout(), task));
	lock.unlock();
	return timeout;
}

/**
 * This method also adds a timer task to the list.
 * It should be called by the user.
 */
bool timer::add_task(timer_task *task)
{
	int64_t next_timeout = _add_task(task);
	if (next_timeout > 0) {
		set_timeout(next_timeout);
	}
	return true;
}

void timer::set_timeout(int64_t timeout)
{
	assert(timeout > 0);

	struct itimerspec its;

	disable_timer();

	/* Start the timer */

	its.it_value.tv_sec = 0;
	its.it_value.tv_nsec = timeout * 1000;
	its.it_interval.tv_sec = 0;
	its.it_interval.tv_nsec = 0;

	struct itimerspec old_value;
	if (timer_settime(timerid, 0, &its, &old_value) == -1) {
		perror("timer_settime");
		exit(1);
	}
	assert(old_value.it_interval.tv_sec == 0
			&& old_value.it_interval.tv_nsec == 0
			&& old_value.it_value.tv_sec == 0
			&& old_value.it_value.tv_nsec == 0);

	enable_timer();
}

void timer::run_tasks()
{
	int64_t curr = get_curr_time_us();
	std::vector<timer_task *> tasks;
	int64_t next_timeout = -1;

	// The current timer has expired. Since the timer is set to timeout
	// once, each timer has its own thread to get its timeout signal,
	// it's guaranteed that the thread won't get any new timeout signal,
	// so it's safe to use lock in the procedure.

	// Get all timeout tasks.
	lock.lock();
	while (!task_list.empty()) {
		timer_task *task = task_list.begin()->second;
		if (task->get_abs_timeout() <= curr) {
			tasks.push_back(task);
			task_list.erase(task_list.begin());
		}
		else
			break;
	}
	lock.unlock();

	// Execute all timeout tasks.
	for (unsigned i = 0; i < tasks.size(); i++) {
		tasks[i]->run();
		delete tasks[i];
	}

	lock.lock();
	if (!task_list.empty())
		next_timeout = task_list.begin()->second->get_abs_timeout() - curr;
	lock.unlock();

	// TODO There might be a race condition here.
	// If a timer task is added to the list and the new task need to timeout
	// earlier than any existing tasks. This invocation of set_timeout() may
	// override the timeout set by the new timer task.
	if (next_timeout > 0)
		set_timeout(next_timeout);
	assert(get_next_timeout() > 0);
}

int64_t timer::get_next_timeout()
{
	struct itimerspec curr_value;
	int ret = timer_gettime(timerid, &curr_value);
	if (ret < 0) {
		perror("timer_gettime");
		assert(0);
	}
	return ((int64_t) curr_value.it_value.tv_sec) * 1000000
		+ curr_value.it_value.tv_nsec / 1000;
}

void timer::print_state()
{
	printf("there are %ld timer tasks and next timeout is in %ldus\n",
			task_list.size(), get_next_timeout());
}

static void periodic_timer_handler(int sig, siginfo_t *si, void *uc)
{
	periodic_timer *t = (periodic_timer *) si->si_value.sival_ptr;
	t->run_task();
}

periodic_timer::periodic_timer(thread *t, timer_task *task)
{
	timer_id = timer_count.inc(1);
	this->t = t;
	this->task = task;

	/**
	 * The code here is copied from the example code in the manual of
	 * timer_create.
	 */
	struct sigevent sev;
	struct sigaction sa;

	/* Establish handler for timer signal */

	sa.sa_flags = SA_SIGINFO;
	sa.sa_sigaction = periodic_timer_handler;
	sigemptyset(&sa.sa_mask);
	if (sigaction(SIGRTMIN, &sa, NULL) == -1) {
		perror("sigaction");
		exit(1);
	}

	memset(&sev, 0, sizeof(sev));
	// This is a Linux feature, but we really want the signal sent to a dedicated
	// thread, so whatever the timer task does, there won't be deadlock.
	sev.sigev_notify = SIGEV_THREAD_ID;
	// It supposes to be sigev_notify_thread_id
	sev._sigev_un._tid = t->get_tid();
	sev.sigev_signo = SIGRTMIN;
	sev.sigev_value.sival_ptr = this;
	if (timer_create(CLOCK_REALTIME, &sev, &timerid) == -1) {
		perror("timer_create");
		exit(1);
	}
}

void periodic_timer::set_timeout(int64_t timeout)
{
	assert(timeout > 0);

	struct itimerspec its;

	disable_timer();

	/* Start the timer */

	its.it_value.tv_sec = 0;
	its.it_value.tv_nsec = timeout * 1000;
	its.it_interval.tv_sec = its.it_value.tv_sec;
	its.it_interval.tv_nsec = its.it_value.tv_nsec;

	struct itimerspec old_value;
	if (timer_settime(timerid, 0, &its, &old_value) == -1) {
		perror("timer_settime");
		exit(1);
	}

	enable_timer();
}

bool periodic_timer::start()
{
	set_timeout(task->get_timeout());

	struct itimerspec curr_value;
	int ret = timer_gettime(timerid, &curr_value);
	if (ret < 0) {
		perror("timer_gettime");
		assert(0);
	}

	return true;
}

bool periodic_timer::is_enabled() const
{
	struct itimerspec curr_value;
	int ret = timer_gettime(timerid, &curr_value);
	if (ret < 0) {
		perror("timer_gettime");
		assert(0);
	}
	return ((int64_t) curr_value.it_interval.tv_sec) * 1000000
		+ curr_value.it_interval.tv_nsec / 1000 > 0;
}

void periodic_timer::run_task()
{
	task->run();
}

atomic_number<long> timer_task::task_count;
atomic_integer timer::timer_count;
atomic_integer periodic_timer::timer_count;
