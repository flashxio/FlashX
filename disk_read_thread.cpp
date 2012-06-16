#include "disk_read_thread.h"

const int MAX_FETCH_REQS = 100;

template<class T>
int io_queue<T>::fetch(T *entries, int num) {
	/* we have to wait for coming requests. */
	pthread_mutex_lock(&empty_mutex);
	while(this->is_empty()) {
		pthread_cond_wait(&empty_cond, &empty_mutex);
	}
	pthread_mutex_unlock(&empty_mutex);

	int ret = bulk_queue<T>::fetch(entries, num);

	/* wake up all threads to send more requests */
	pthread_cond_broadcast(&full_cond);

	return ret;
}

/**
 * This is a blocking version.
 * It adds all entries to the queue. If the queue is full,
 * wait until it can add all entries.
 */
template<class T>
int io_queue<T>::add(T *entries, int num) {
	int orig_num = num;

	while (num > 0) {
		int ret = bulk_queue<T>::add(entries, num);
		entries += ret;
		num -= ret;
		/* signal the thread of reading disk to wake up. */
		pthread_cond_signal(&empty_cond);

		pthread_mutex_lock(&full_mutex);
		while (this->is_full()) {
			pthread_cond_wait(&full_cond, &full_mutex);
		}
		pthread_mutex_unlock(&full_mutex);
	}
	return orig_num;
}

/* just call the callback of the initiator. */
class initiator_callback: public callback
{
public:

	int invoke(io_request *rq) {
		io_interface *io = rq->get_io();
		/* 
		 * after a request is processed,
		 * we need to notify the initiator thread.
		 * It's possible the initiator thread is itself,
		 * we need to stop the infinite loop.
		 */
		if (io->get_callback())
			io->get_callback()->invoke(rq);
		return 0;
	}
};

disk_read_thread::disk_read_thread(const char *name,
		long size): queue(1024) {
	aio = new async_io(&name, 1, size);
	aio->set_callback(new initiator_callback());

	int ret = pthread_create(&id, NULL, process_requests, (void *) this);
	if (ret) {
		perror("pthread_create");
		exit(1);
	}
}

void disk_read_thread::run() {
	aio->init();
	io_request reqs[MAX_FETCH_REQS];
	while (true) {

		/* 
		 * this is the only thread that fetch requests
		 * from the queue.
		 * get all requests in the queue
		 */
		int num = queue.fetch(reqs, MAX_FETCH_REQS);

		aio->access(reqs, num);
	}
	// TODO I need to call cleanup() of aio.
}

void *process_requests(void *arg)
{
	disk_read_thread *thread = (disk_read_thread *) arg;
	thread->run();
	return NULL;
}
