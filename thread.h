#ifndef __MY_THREAD_H__
#define __MY_THREAD_H__

#include <pthread.h>

class thread
{
	int node_id;
	pthread_t id;

	volatile bool _is_running;

	bool is_activate;
	pthread_mutex_t mutex;
	pthread_cond_t cond;
public:
	thread(int node_id = -1) {
		this->node_id = node_id;
		is_activate = false;
		_is_running = true;
		pthread_mutex_init(&mutex, NULL);
		pthread_cond_init(&cond, NULL);
	}

	void activate() {
		pthread_mutex_lock(&mutex);
		is_activate = true;
		pthread_cond_signal(&cond);
		pthread_mutex_unlock(&mutex);
	}

	void wait() {
		pthread_mutex_lock(&mutex);
		is_activate = false;
		while (!is_activate) {
			pthread_cond_wait(&cond, &mutex);
		}
		pthread_mutex_unlock(&mutex);
	}

	bool is_running() const {
		return _is_running;
	}

	int get_node_id() const {
		return node_id;
	}

	void stop() {
		// This is the only place where `is_running' is changed.
		_is_running = false;
	}
	void start();
	virtual void run() = 0;
};

#endif
