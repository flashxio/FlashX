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

#include <signal.h>
#include <time.h>

#include "timer.h"
#include "thread.h"

static void handler(int sig, siginfo_t *si, void *uc)
{
	// We can't call the flush requests to wake up I/O threads.
	// There might be a deadlock if the lock has been held by the thread
	// interrupted by the signal.

	timer *t = (timer *) thread::get_curr_thread()->get_user_data();
	t->run_tasks();
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
	sev.sigev_value.sival_ptr = &timerid;
	if (timer_create(CLOCK_REALTIME, &sev, &timerid) == -1) {
		perror("timer_create");
		exit(1);
	}
}

timer::timer(thread *t)
{
	this->t = t;
	// TODO I should change this. There is no guarantee that the current
	// thread doesn't have any user data.
	assert(t->get_user_data() == NULL);
	t->set_user_data(this);

	init_timer();
}

timer::timer()
{
	t = new timer_thread();
	t->start();
	t->set_user_data(this);

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

void timer::disable_timer()
{
	sigset_t mask;

	sigemptyset(&mask);
	sigaddset(&mask, SIGRTMIN);
	if (sigprocmask(SIG_SETMASK, &mask, NULL) == -1) {
		perror("sigprocmask");
		exit(1);
	}
}

void timer::enable_timer()
{
	sigset_t mask;

	sigemptyset(&mask);
	sigaddset(&mask, SIGRTMIN);
	if (sigprocmask(SIG_UNBLOCK, &mask, NULL) == -1) {
		perror("sigprocmask");
		exit(1);
	}
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

	if (timer_settime(timerid, 0, &its, NULL) == -1) {
		perror("timer_settime");
		exit(1);
	}

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
		if (tasks[i]->repeat()) {
			tasks[i]->renew();
			_add_task(tasks[i]);
		}
		else
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

atomic_number<long> timer_task::task_count;
