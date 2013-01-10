#include <stdio.h>
#include <stdlib.h>

#include "thread.h"

static void *thread_run(void *arg)
{
	thread *t = (thread *) arg;
	while (t->is_running()) {
		t->wait();
		t->run();
	}
	return NULL;
}

void thread::start()
{
	int ret = pthread_create(&id, NULL, thread_run, (void *) this);
	if (ret) {
		perror("pthread_create");
		exit(1);
	}
}
