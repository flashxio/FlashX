#ifndef __MESSAGING_H__
#define __MESSAGING_H__

#include <pthread.h>
#include <numa.h>
#include <assert.h>

#include "common.h"

class thread_private;
class io_request
{
	char *buf;
	off_t offset;
	ssize_t size: 32;
	int access_method: 1;
	thread_private *thread;
	void *priv;
public:
	io_request() {
		init(NULL, 0, 0, READ, NULL);
	}

	io_request(char *buf, off_t off, ssize_t size,
			int access_method, thread_private *t, void *priv = NULL) {
		init(buf, off, size, access_method, t, priv);
	}

	void init(char *buf, off_t off, ssize_t size,
			int access_method, thread_private *t, void *priv = NULL) {
		assert(off >= 0);
		this->buf = buf;
		this->offset = off;
		this->size = size;
		this->thread = t;
		this->access_method = access_method;
		this->priv = priv;
	}

	int get_access_method() {
		return access_method;
	}

	thread_private *get_thread() {
		return thread;
	}

	char *get_buf() {
		return buf;
	}

	off_t get_offset() {
		return offset;
	}

	ssize_t get_size() {
		return size;
	}

	void *get_priv() {
		return priv;
	}

	void set_priv(void *priv) {
		this->priv = priv;
	}
};

class io_reply
{
	char *buf;
	off_t offset;
	ssize_t size: 32;
	int success: 1;
	int status: 16;
	int access_method: 1;
	void init(char *buf, off_t off, ssize_t size, int success,
			int status, int access_method) {
		this->buf = buf;
		this->offset = off;
		this->size = size;
		this->success = success;
		this->status = status;
		this->access_method = access_method;
	}
public:
	io_reply() {
		init(NULL, 0, 0, 0, 0, READ);
	}

	io_reply(io_request *req, int success, int status) {
		init(req->get_buf(), req->get_offset(), req->get_size(),
					success, status, req->get_access_method());
	}

	int get_status() {
		return status;
	}

	bool is_success() {
		return success;
	}

	char *get_buf() {
		return buf;
	}

	off_t get_offset() {
		return offset;
	}

	ssize_t get_size() {
		return size;
	}

	int get_access_method() {
		return access_method;
	}
};

/**
 * this is a thread-safe FIFO queue.
 * It supports bulk operations.
 */
template<class T>
class bulk_queue
{
	T *buf;
	volatile int size;
	int start;
	int num_entries;
	pthread_spinlock_t _lock;
public:
	bulk_queue(int size) {
		buf = new T[size];
		this->size = size;
		start = 0;
		num_entries = 0;
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

	~bulk_queue() {
		pthread_spin_destroy(&_lock);
		delete [] buf;
	}

	int fetch(T *entries, int num);

	int add(T *entries, int num);

	int get_num_entries() {
		return num_entries;
	}

	bool is_full() {
		return num_entries == size;
	}

	bool is_empty() {
		return num_entries == 0;
	}
};

template<class T>
class msg_sender
{
	T *buf;
	int buf_size;		// the max number of messages that can be buffered
	int num_current;	// the current number of messages in the buffer.
	bulk_queue<T> **dest_queues;
	int num_queues;
public:
	/**
	 * buf_size: the number of messages that can be buffered in the sender.
	 */
	msg_sender(int buf_size, bulk_queue<T> **queues, int num_queues);

	~msg_sender() {
		numa_free(buf, sizeof(T) * buf_size);
		numa_free(dest_queues, sizeof(bulk_queue<T> *) * num_queues);
	}

	int flush();

	int send_cached(T *msg);
};

#endif
