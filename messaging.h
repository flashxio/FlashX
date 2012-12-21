#ifndef __MESSAGING_H__
#define __MESSAGING_H__

#include <pthread.h>
#include <numa.h>
#include <assert.h>
#include <sys/uio.h>

#include "common.h"
#include "container.h"
#include "parameters.h"

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

class multibuf_io_request
{
	struct iovec vecs[MAX_NUM_IOVECS];
	int num_bufs;
	off_t offset;
	int access_method: 1;
	io_interface *io;
	void *priv;
public:
	multibuf_io_request() {
		memset(vecs, sizeof(vecs), 0);
		num_bufs = 0;
		offset = 0;
		access_method = 0;
		io = NULL;
		priv = NULL;
	}

	multibuf_io_request(const struct iovec vecs[], int num_vecs,
			off_t offset, int access_method, io_interface *io,
			void *priv = NULL) {
		memcpy(this->vecs, vecs, sizeof(vecs[0]) * num_vecs);
		this->num_bufs = num_vecs;
		this->offset = offset;
		this->access_method = access_method & 0x1;
		this->io = io;
		this->priv = priv;
	}

	void add_buf(char *buf, int size) {
		assert(num_bufs < MAX_NUM_IOVECS);
		vecs[num_bufs].iov_base = buf;
		vecs[num_bufs].iov_len = size;
		num_bufs++;
	}

	int get_num_bufs() const {
		return num_bufs;
	}

	char *get_buf(int idx) const {
		return (char *) vecs[idx].iov_base;
	}

	int get_buf_size(int idx) const {
		return vecs[idx].iov_len;
	}

	const struct iovec &get(int idx) const {
		return vecs[idx];
	}

	int get_access_method() const {
		return access_method & 0x1;
	}

	io_interface *get_io() const {
		return io;
	}

	off_t get_offset() const {
		return offset;
	}

	ssize_t get_size() const {
		ssize_t size = 0;
		for (int i = 0; i < num_bufs; i++)
			size += vecs[i].iov_len;
		return size;
	}

	void *get_priv() const {
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
