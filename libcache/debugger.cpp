#include <signal.h>
#include <stdio.h>

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
