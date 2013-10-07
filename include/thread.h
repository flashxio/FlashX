#ifndef __MY_THREAD_H__
#define __MY_THREAD_H__

#include <pthread.h>

#include <string>

class thread
{
	static pthread_key_t thread_key;

	int node_id;
	pthread_t id;
	bool blocking;
	std::string name;

	volatile bool _is_sleeping;
	volatile bool _is_running;
	volatile bool _has_exit;
	volatile bool _is_activated;

	pthread_mutex_t mutex;
	pthread_cond_t cond;

	friend void init_thread_class();
	friend void *thread_run(void *arg);
public:
	thread(std::string name, int node_id = -1, bool blocking = true) {
		thread_class_init();

		this->name = name;
		this->node_id = node_id;
		this->blocking = blocking;
		this->id = 0;

		_is_activated = false;
		_has_exit = false;
		_is_running = true;
		_is_sleeping = true;

		pthread_mutex_init(&mutex, NULL);
		pthread_cond_init(&cond, NULL);
	}

	virtual ~thread() {
	}

	void activate() {
		_is_activated = true;
		if (_is_sleeping) {
			pthread_cond_signal(&cond);
		}
	}

	void wait() {
		if (!_is_activated) {
			pthread_mutex_lock(&mutex);
			_is_sleeping = true;
			while (!_is_activated && _is_running) {
				int ret = pthread_cond_wait(&cond, &mutex);
				if (ret)
					perror("pthread_cond_wait");
			}
			_is_sleeping = false;
			pthread_mutex_unlock(&mutex);
		}
		_is_activated = false;
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

	void exit() {
		_has_exit = true;
	}

	bool has_exit() const {
		return _has_exit;
	}

	void join() {
		pthread_join(id, NULL);
#ifdef DEBUG
		printf("stop thread %s\n", name.c_str());
#endif
	}

	void start();
	virtual void run() = 0;
	virtual void init() {
	}
	virtual void cleanup() {
	}

	/**
	 * This is to initialize the thread class instead of a single thread.
	 */
	static void thread_class_init();

	static thread *get_curr_thread();

	/**
	 * This creates a thread instance to represent the current thread context.
	 * It is used when the current thread isn't created by the thread class.
	 */
	static thread *represent_thread(int node_id);
};

#endif
