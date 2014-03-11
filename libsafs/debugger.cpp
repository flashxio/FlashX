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
#include <stdio.h>

#include <boost/foreach.hpp>

#include "debugger.h"

static bool enable_debug = false;

bool is_debug_enabled()
{
	return enable_debug;
}

static void enable_debug_handler(int sig, siginfo_t *si, void *uc)
{
	printf("debug mode is enabled\n");
	debug.run();
}

static void set_enable_debug_signal()
{
	struct sigaction sa;

	/* Establish handler for timer signal */

	sa.sa_flags = SA_SIGINFO;
	sa.sa_sigaction = enable_debug_handler;
	sigemptyset(&sa.sa_mask);
	if (sigaction(SIGUSR1, &sa, NULL) == -1) {
		perror("sigaction");
		exit(1);
	}
}

debugger::debugger()
{
	pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	set_enable_debug_signal();
}

debugger::~debugger()
{
	typedef std::map<int, debug_task *> task_map_t;
	BOOST_FOREACH(task_map_t::value_type &p, tasks) {
		delete p.second;
	}
}

void debugger::run()
{
	pthread_spin_lock(&lock);
	std::map<int, debug_task *> task_copies = tasks;
	pthread_spin_unlock(&lock);
	for (std::map<int, debug_task *>::const_iterator it
			= task_copies.begin(); it != task_copies.end(); it++) {
		it->second->run();
	}
}

debugger debug;
