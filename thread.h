#ifndef __MY_THREAD_H__
#define __MY_THREAD_H__

#include <pthread.h>

#include <string>

class thread
{
	int node_id;
	pthread_t id;
	bool blocking;
	std::string name;

	volatile bool _is_running;

	bool is_activate;
	pthread_mutex_t mutex;
	pthread_cond_t cond;
public:
	thread(std::string name, int node_id = -1, bool blocking = true) {
		this->name = name;
		this->node_id = node_id;
		this->blocking = blocking;
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
		while (!is_activate && _is_running) {
			pthread_cond_wait(&cond, &mutex);
		}
		pthread_mutex_unlock(&mutex);
	}

	bool is_running() const {
		return _is_running;
	}

	bool is_blocking() const {
		return blocking;
	}

	int get_node_id() const {
		return node_id;
	}

	void stop() {
		pthread_mutex_lock(&mutex);
		// This is the only place where `is_running' is changed.
		_is_running = false;
		pthread_cond_signal(&cond);
		pthread_mutex_unlock(&mutex);
	}

	void join() {
		pthread_join(id, NULL);
		printf("stop thread %s\n", name.c_str());
	}

	void start();
	virtual void run() = 0;
};

#endif
