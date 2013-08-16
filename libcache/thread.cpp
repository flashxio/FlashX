#include <stdio.h>
#include <stdlib.h>
#include <numa.h>

#include "thread.h"
#include "common.h"

static void *thread_run(void *arg)
{
	thread *t = (thread *) arg;
	int node_id = t->get_node_id();
	if (node_id >= 0)
		bind2node_id(node_id);
	while (t->is_running()) {
		t->run();
		if (t->is_blocking())
			t->wait();
	}
	t->exit();
	return NULL;
}

void thread::start()
{
	int ret = pthread_create(&id, NULL, thread_run, (void *) this);
	if (ret) {
		perror("pthread_create");
		::exit(1);
	}
}
