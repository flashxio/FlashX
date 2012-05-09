#ifndef __DISK_READ_THREAD_H__
#define __DISK_READ_THREAD_H__

#include "messaging.h"
#include "aio_private.h"


void *process_requests(void *arg);

template<class T>
class io_queue: public bulk_queue<T>
{
	/* when the queue becomes empty */
	pthread_cond_t empty_cond;
	pthread_mutex_t empty_mutex;

	/* when the queue becomes full */
	pthread_cond_t full_cond;
	pthread_mutex_t full_mutex;
public:
	io_queue(int size): bulk_queue<T>(size) {
		pthread_mutex_init(&empty_mutex, NULL);
		pthread_cond_init(&empty_cond, NULL);
		pthread_mutex_init(&full_mutex, NULL);
		pthread_cond_init(&full_cond, NULL);
	}

	virtual int fetch(T *entries, int num);

	virtual int add(T *entries, int num);
};

class disk_read_thread
{
	io_queue<io_request> queue;
	pthread_t id;
	aio_private *aio;

public:
	disk_read_thread(const char *name, long size);

	io_queue<io_request> *get_queue() {
		return &queue;
	}

	~disk_read_thread() {
		delete aio;
	}

	void run();
};

#endif
