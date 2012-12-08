#ifndef __MESSAGING_H__
#define __MESSAGING_H__

#include <pthread.h>
#include <numa.h>
#include <assert.h>

#include "common.h"
#include "container.h"

class io_request;
typedef void (*request_callback_t)(io_request *);

class io_interface;
class io_request
{
	char *buf;
	off_t offset;
	ssize_t size: 32;
	int access_method: 1;
	io_interface *io;
	void *priv;
	io_request *next;
	request_callback_t cb;
public:
	io_request() {
		init(NULL, 0, 0, READ, NULL);
	}

	io_request(char *buf, off_t off, ssize_t size,
			int access_method, io_interface *io, void *priv = NULL) {
		init(buf, off, size, access_method, io, priv);
	}

	void init(char *buf, off_t off, ssize_t size,
			int access_method, io_interface *io, void *priv = NULL) {
		assert(off >= 0);
		this->buf = buf;
		this->offset = off;
		this->size = size;
		this->io = io;
		this->access_method = access_method & 0x1;
		this->priv = priv;
		next = NULL;
		cb = NULL;
	}

	void set_cb(request_callback_t cb) {
		this->cb = cb;
	}

	request_callback_t get_cb() const {
		return cb;
	}

	int get_access_method() {
		return access_method & 0x1;
	}

	io_interface *get_io() {
		return io;
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

	io_request *get_next_req() const {
		return next;
	}

	void set_next_req(io_request *next) {
		this->next = next;
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

template<class T>
class msg_sender
{
	T *buf;
	int buf_size;		// the max number of messages that can be buffered
	int num_current;	// the current number of messages in the buffer.
	thread_safe_FIFO_queue<T> **dest_queues;
	int num_queues;
public:
	/**
	 * buf_size: the number of messages that can be buffered in the sender.
	 */
	msg_sender(int buf_size, thread_safe_FIFO_queue<T> **queues, int num_queues);

	~msg_sender() {
		numa_free(buf, sizeof(T) * buf_size);
		numa_free(dest_queues, sizeof(thread_safe_FIFO_queue<T> *) * num_queues);
	}

	int num_msg() {
		return num_current;
	}

	int flush();

	int send_cached(T *msg);
};

#endif
