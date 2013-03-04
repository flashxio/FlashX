#include <stdio.h>
#include <stdlib.h>
#include <numa.h>

#include "thread.h"

static void *thread_run(void *arg)
{
	thread *t = (thread *) arg;
	int node_id = t->get_node_id();
	if (node_id >= 0)
		numa_run_on_node(node_id);
	while (t->is_running()) {
		if (t->is_blocking())
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
