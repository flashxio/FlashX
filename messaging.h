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

/**
 * This class contains the info of an IO request.
 */
class io_request
{
	off_t offset;
	io_interface *io;
	io_request *orig;
	void *priv;

	int access_method: 1;
	int num_bufs: 15;
	// Is the request part of a request?
	int partial: 1;
	int vec_capacity: 15;

	/* 
	 * This is to protect the object from being removed
	 * while others are still using it.
	 */
	volatile int refcnt;

	struct iovec *vec_pointer;
	struct iovec embedded_vecs[NUM_EMBEDDED_IOVECS];
	io_request *next;
	volatile ssize_t completed_size;

	bool use_embedded() const {
		return vec_pointer == embedded_vecs;
	}

	void assign(io_request &req);

public:
	io_request() {
		init(-1, NULL, READ, NULL, NULL);
	}

	io_request(off_t off, io_interface *io, int access_method,
			io_request *orig = NULL, void *priv = NULL) {
		init(off, io, access_method, orig, priv);
	}

	io_request(char *buf, off_t off, ssize_t size, int access_method,
			io_interface *io, io_request *orig = NULL, void *priv = NULL) {
		init(buf, off, size, access_method, io, orig, priv);
	}

	io_request(io_request &req) {
		assign(req);
	}

	io_request &operator=(io_request &req) {
		assign(req);
		return *this;
	}

	~io_request() {
		if (vec_pointer != embedded_vecs)
			delete [] vec_pointer;
	}

	void init(char *buf, off_t off, ssize_t size, int access_method,
			io_interface *io, io_request *orig = NULL, void *priv = NULL) {
		init(off, io, access_method, orig, priv);
		add_buf(buf, size);
	}

	void init(off_t off, io_interface *io, int access_method,
			io_request *orig = NULL, void *priv = NULL) {
		this->offset = off;
		this->io = io;
		this->access_method = access_method & 0x1;
		this->priv = priv;
		this->partial = 0;
		this->completed_size = 0;
		this->orig = orig;
		this->refcnt = 0;
		memset(embedded_vecs, 0,
				sizeof(embedded_vecs[0]) * NUM_EMBEDDED_IOVECS);
		num_bufs = 0;
		vec_pointer = embedded_vecs;
		vec_capacity = NUM_EMBEDDED_IOVECS;
		next = NULL;
	}

	int get_access_method() const {
		return access_method & 0x1;
	}

	io_interface *get_io() const {
		return io;
	}

	io_request *get_orig() const {
		return orig;
	}

	void set_orig(io_request *orig) {
		this->orig = orig;
	}

	void set_offset(off_t offset) {
		this->offset = offset;
	}

	off_t get_offset() const {
		return offset;
	}

	void *get_priv() const {
		return priv;
	}

	void set_priv(void *priv) {
		this->priv = priv;
	}

	bool is_empty() const {
		return num_bufs == 0;
	}

	void clear() {
		memset((void *) vec_pointer, 0, sizeof(vec_pointer[0]) * vec_capacity);
		num_bufs = 0;
		set_offset(-1);
	}

	void add_buf(char *buf, int size);
	void add_buf_front(char *buf, int size);

	int get_num_bufs() const {
		return num_bufs;
	}

	/**
	 * By default, we get the first buffer. This makes sense
	 * for a single buffer request.
	 */
	char *get_buf(int idx = 0) const {
		return (char *) vec_pointer[idx].iov_base;
	}

	int get_buf_size(int idx) const {
		return vec_pointer[idx].iov_len;
	}

	const struct iovec &get(int idx) const {
		return vec_pointer[idx];
	}

	const struct iovec *get_vec() const {
		return vec_pointer;
	}

	ssize_t get_size() const {
		ssize_t size = 0;
		for (int i = 0; i < num_bufs; i++)
			size += vec_pointer[i].iov_len;
		return size;
	}

	io_request *get_next_req() const {
		return next;
	}

	void set_next_req(io_request *next) {
		this->next = next;
	}

	int inc_complete_count() {
		return __sync_add_and_fetch(&refcnt, 1);
	}

	int dec_complete_count() {
		return __sync_sub_and_fetch(&refcnt, 1);
	}

	void wait4unref() {
		while (refcnt > 0) {}
	}

	/**
	 * Maintain the completed size in this request.
	 * If the request is complete, return true;
	 */
	bool complete_size(ssize_t completed) {
		ssize_t res = __sync_add_and_fetch(&completed_size, completed);
		ssize_t size = get_size();
		assert(res <= size);
		return res == size;
	}

	void set_partial(bool partial) {
		this->partial = partial ? 1 : 0;
	}

	bool is_partial() const {
		return this->partial;
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
