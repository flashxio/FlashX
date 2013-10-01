#include <stdio.h>
#include <stdlib.h>
#include <numa.h>
#include <pthread.h>

#include "thread.h"
#include "common.h"

static void *thread_run(void *arg)
{
	thread *t = (thread *) arg;
	int node_id = t->get_node_id();
	if (node_id >= 0)
		bind2node_id(node_id);
	t->init();
	while (t->is_running()) {
		t->run();
		if (t->is_blocking())
			t->wait();
	}
	t->cleanup();
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

static pthread_once_t once_control = PTHREAD_ONCE_INIT;

void init_thread_class()
{
	printf("init thread key\n");
	pthread_key_create(&thread::thread_key, NULL);
}

void thread::thread_class_init()
{
	pthread_once(&once_control, init_thread_class);
}

thread *thread::get_curr_thread()
{
	thread *curr = (thread *) pthread_getspecific(thread_key);
	return curr;
}

pthread_key_t thread::thread_key;
