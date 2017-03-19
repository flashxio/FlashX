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
#include <stdio.h>

#include <boost/foreach.hpp>
#include <system_error>

#include "debugger.h"

namespace safs
{

static bool enable_debug = false;

bool is_debug_enabled()
{
	return enable_debug;
}

static void enable_debug_handler(int sig, siginfo_t *si, void *uc)
{
#if 0
	printf("debug mode is enabled\n");
	debug.run();
#endif
}

static void set_enable_debug_signal()
{
	struct sigaction sa;

	/* Establish handler for timer signal */

	sa.sa_flags = SA_SIGINFO;
	sa.sa_sigaction = enable_debug_handler;
	sigemptyset(&sa.sa_mask);
	if (sigaction(SIGUSR1, &sa, NULL) == -1)
		throw std::system_error(std::make_error_code((std::errc) errno),
				"sigaction error");
}

debugger::debugger()
{
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
	lock.lock();
	std::map<int, debug_task *> task_copies = tasks;
	lock.unlock();
	for (std::map<int, debug_task *>::const_iterator it
			= task_copies.begin(); it != task_copies.end(); it++) {
		it->second->run();
	}
}

}
