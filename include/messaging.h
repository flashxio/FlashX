#ifndef __MESSAGING_H__
#define __MESSAGING_H__

#include <pthread.h>
#include <numa.h>
#include <assert.h>
#include <sys/uio.h>

#include "common.h"
#include "container.h"
#include "parameters.h"
#include "concurrency.h"

class io_interface;
io_interface *get_io(int idx);

class io_request;
class thread_safe_page;

class io_buf
{
	union {
		thread_safe_page *p;
		void *buf;
	} u;
	unsigned int size: 31;
	unsigned int is_page: 1;
public:
	io_buf() {
		u.buf = NULL;
		size = 0;
		is_page = 0;
	}

	void init(void *p, int size, bool is_page) {
		u.buf = p;
		if (is_page)
			assert(size == PAGE_SIZE);
		this->size = size;
		this->is_page = is_page;
	}

	void init(thread_safe_page *p);
	void init(void *buf, int size) {
		u.buf = buf;
		this->size = size;
		is_page = 0;
	}

	void *get_buf() const;

	int get_size() const {
		return size;
	}

	thread_safe_page *get_page() const {
		assert(is_page);
		return u.p;
	}
};

struct io_req_extension
{
	static atomic_unsigned_integer num_creates;
	io_interface *io;
	io_request *orig;
	void *priv;
	void *user_data;

	int num_bufs: 16;
	// Is the request part of a request?
	int partial: 1;
	int vec_capacity: 15;

	/* 
	 * This is to protect the object from being removed
	 * while others are still using it.
	 */
	volatile short refcnt;

	io_buf *vec_pointer;
	io_buf embedded_vecs[NUM_EMBEDDED_IOVECS];
	io_request *next;
	volatile ssize_t completed_size;

	void init() {
		this->io = NULL;
		this->orig = NULL;
		this->priv = NULL;
		this->user_data = NULL;
		this->num_bufs = 0;
		this->partial = 0;
		this->vec_capacity = NUM_EMBEDDED_IOVECS;
		this->refcnt = 0;
		vec_pointer = embedded_vecs;
		next = NULL;
		this->completed_size = 0;
	}

	io_req_extension() {
		num_creates.inc(1);
		init();
	}

	~io_req_extension() {
		if (vec_pointer != embedded_vecs)
			delete [] vec_pointer;
	}

	void add_io_buf(const io_buf &buf);
	void add_buf(char *buf, int size, bool is_page);
	void add_buf_front(char *buf, int size, bool is_page);
};

/**
 * This class contains the request info.
 * If it's in an extended form, it is created by the global cache
 * and is used within a NUMA machine.
 * It is decided when the request is created whether or not it is
 * an extended request. It can't be changed afterwards.
 */
class io_request
{
	off_t offset: 40;
	// This flag is initialized when the object is created and
	// can't be changed manually.
	// The only case that the flag of a request is changed is when
	// another request object is copied to this request object.
	// What needs to be guaranteed is that
	// the flag is true when buf_addr points to the extension object;
	// the flag is false when buf_addr points to the real buffer.
	unsigned int extended: 1;
	unsigned int high_prio: 1;
	// Is this synchronous IO?
	unsigned int sync: 1;
	/*
	 * The NUMA node id where the buffers of the request are allocated.
	 */
	unsigned int node_id: 5;
	unsigned int io_idx: 16;

	unsigned int access_method: 1;
	unsigned long buf_size: 15;
	// Linux uses 48 bit for addresses.
	unsigned long buf_addr: 48;

	io_req_extension *get_extension() const {
		assert(is_extended_req() && buf_addr);
		unsigned long addr = buf_addr;
		return (io_req_extension *) addr;
	}

	bool use_embedded() const {
		return get_extension()->vec_pointer == get_extension()->embedded_vecs;
	}

public:
	io_request() {
		extended = 0;
		buf_addr = 0;
		high_prio = 1;
		sync = false;
	}

	io_request(bool extended) {
		this->extended = extended;
		buf_addr = 0;
		high_prio = 1;
		sync = false;
		if (extended) {
			buf_addr = (long) new io_req_extension();
		}
	}

	io_request(char *buf, off_t off, ssize_t size, int access_method,
			io_interface *io, int node_id, bool sync = false) {
		extended = 0;
		this->sync = sync;
		init(buf, off, size, access_method, io, node_id);
	}

	io_request(off_t off, io_interface *io, int access_method, int node_id,
			io_request *orig, void *priv, bool sync = false) {
		extended = 1;
		buf_addr = (long) new io_req_extension();
		this->sync = sync;
		init(off, io, access_method, node_id, orig, priv, NULL);
	}

	io_request(char *buf, off_t off, ssize_t size, int access_method,
			io_interface *io, int node_id, io_request *orig,
			void *priv, bool sync = false) {
		extended = 1;
		buf_addr = (long) new io_req_extension();
		this->sync = sync;
		init(buf, off, size, access_method, io, node_id, orig, priv, NULL);
	}

	io_request(io_request &req) {
		// We need to free its own extension first.
		if (this->is_extended_req())
			delete this->get_extension();
		memcpy(this, &req, sizeof(req));
		if (req.is_extended_req()) {
			req.extended = 0;
			req.buf_addr = 0;
		}
	}

	io_request &operator=(io_request &req) {
		if (this->is_extended_req())
			delete this->get_extension();
		memcpy(this, &req, sizeof(req));
		if (req.is_extended_req()) {
			req.extended = 0;
			req.buf_addr = 0;
		}
		return *this;
	}

	~io_request() {
		if (is_extended_req() && get_extension())
			delete get_extension();
	}

	void init(const io_request &req) {
		this->sync = req.sync;
		if (!req.is_extended_req()) {
			this->init(req.get_buf(), req.get_offset(), req.get_size(),
					req.get_access_method(), req.get_io(), req.get_node_id());
		}
		else if (this->is_extended_req()) {
			this->init(req.get_offset(), req.get_io(), req.get_access_method(),
					req.get_node_id(), req.get_orig(), req.get_priv(), req.get_user_data());
			for (int i = 0; i < req.get_num_bufs(); i++) {
				this->add_io_buf(req.get_io_buf(i));
			}
		}
		else {
			// The last case is that this request doesn't have extension,
			// but the given request has extension. This request can't keep
			// all information in the given request.
			assert(!this->is_extended_req() && req.is_extended_req());
		}
	}

	void init() {
		if (is_extended_req()) {
			io_req_extension *ext = get_extension();
			assert(ext);
			ext->init();
			offset = 0;
			high_prio = 0;
			sync = 0;
			node_id = 0;
			io_idx = 0;
			access_method = 0;
			buf_size = 0;
		}
		else {
			offset = 0;
			extended = 0;
			high_prio = 0;
			sync = 0;
			node_id = 0;
			io_idx = 0;
			access_method = 0;
			buf_size = 0;
			buf_addr = 0;
		}
	}

	void init(char *buf, off_t off, ssize_t size, int access_method,
			io_interface *io, int node_id);

	void init(char *buf, off_t off, ssize_t size, int access_method,
			io_interface *io, int node_id, io_request *orig, void *priv,
			void *user_data) {
		init(off, io, access_method, node_id, orig, priv, user_data);
		add_buf(buf, size);
	}

	void init(off_t off, io_interface *io, int access_method, int node_id,
			io_request *orig, void *priv, void *user_data) {
		assert(is_extended_req());
		io_request::init(NULL, off, 0, access_method, io, node_id);
		io_req_extension *ext = get_extension();
		ext->io = io;
		ext->priv = priv;
		ext->orig = orig;
		ext->user_data = user_data;
	}

	bool is_sync() const {
		return sync;
	}

	bool is_extended_req() const {
		return extended;
	}

	off_t get_offset() const {
		return offset;
	}

	void set_offset(off_t offset) {
		this->offset = offset;
	}

	int get_access_method() const {
		return access_method & 0x1;
	}

	io_interface *get_io() const {
		if (is_extended_req()) {
			return get_extension()->io;
		}
		else
			return ::get_io(io_idx);
	}

	int get_node_id() const {
		return node_id;
	}

	void set_node_id(int node_id) {
		this->node_id = node_id;
	}

	bool is_high_prio() const {
		return (high_prio & 0x1) == 1;
	}

	void set_high_prio(bool high_prio) {
		this->high_prio = high_prio;
	}

	/**
	 * The requested data is inside a page on the disk.
	 */
	bool within_1page() const {
		return get_offset() + get_size() <= ROUND_PAGE(get_offset()) + PAGE_SIZE;
	}

	io_request *get_orig() const {
		return get_extension()->orig;
	}

	void set_orig(io_request *orig) {
		get_extension()->orig = orig;
	}

	void *get_user_data() const {
		return get_extension()->user_data;
	}

	void set_user_data(void *data) {
		get_extension()->user_data = data;
	}

	void *get_priv() const {
		return get_extension()->priv;
	}

	void set_priv(void *priv) {
		get_extension()->priv = priv;
	}

	bool is_empty() const {
		return get_extension()->num_bufs == 0;
	}

	bool is_valid() const {
		return get_offset() != -1;
	}

	ssize_t get_size() const {
		if (!is_extended_req()) {
			return buf_size;
		}
		else {
			ssize_t size = 0;
			for (int i = 0; i < get_extension()->num_bufs; i++)
				size += get_extension()->vec_pointer[i].get_size();
			return size;
		}
	}

	/**
	 * By default, we get the first buffer. This makes sense
	 * for a single buffer request.
	 */
	char *get_buf(int idx = 0) const {
		if (!is_extended_req()) {
			unsigned long addr = buf_addr;
			return (char *) addr;
		}
		else {
			return (char *) get_extension()->vec_pointer[idx].get_buf();
		}
	}

	thread_safe_page *get_page(int idx) const {
		return get_extension()->vec_pointer[idx].get_page();
	}

	void add_buf(char *buf, int size) {
		get_extension()->add_buf(buf, size, 0);
	}

	void add_page(thread_safe_page *p) {
		get_extension()->add_buf((char *) p, PAGE_SIZE, 1);
	}

	void add_io_buf(const io_buf &buf) {
		get_extension()->add_io_buf(buf);
	}

	void add_buf_front(char *buf, int size) {
		get_extension()->add_buf_front(buf, size, 0);
	}

	void add_page_front(thread_safe_page *p) {
		get_extension()->add_buf_front((char *) p, PAGE_SIZE, 1);
	}

	int get_num_bufs() const {
		if (is_extended_req())
			return get_extension()->num_bufs;
		else
			return 1;
	}

	int get_buf_size(int idx) const {
		return get_extension()->vec_pointer[idx].get_size();
	}

	const io_buf &get_io_buf(int idx) const {
		return get_extension()->vec_pointer[idx];
	}

	const struct iovec get(int idx) const {
		struct iovec ret;
		struct io_req_extension *ext = get_extension();
		ret.iov_base = ext->vec_pointer[idx].get_buf();
		ret.iov_len = ext->vec_pointer[idx].get_size();
		return ret;
	}

	const int get_vec(struct iovec *vec, int num) const {
		num = min(get_num_bufs(), num);
		for (int i = 0; i < num; i++) {
			struct io_req_extension *ext = get_extension();
			vec[i].iov_base = ext->vec_pointer[i].get_buf();
			vec[i].iov_len = ext->vec_pointer[i].get_size();
		}
		return num;
	}

	io_request *get_next_req() const {
		return get_extension()->next;
	}

	void set_next_req(io_request *next) {
		get_extension()->next = next;
	}

	int inc_complete_count() {
		return __sync_add_and_fetch(&get_extension()->refcnt, 1);
	}

	int dec_complete_count() {
		return __sync_sub_and_fetch(&get_extension()->refcnt, 1);
	}

	void wait4unref() {
		while (get_extension()->refcnt > 0) {}
	}

	/**
	 * Maintain the completed size in this request.
	 * If the request is complete, return true;
	 */
	bool complete_size(ssize_t completed) {
		ssize_t res = __sync_add_and_fetch(&get_extension()->completed_size, completed);
		ssize_t size = get_size();
		assert(res <= size);
		return res == size;
	}

	void set_partial(bool partial) {
		get_extension()->partial = partial ? 1 : 0;
	}

	bool is_partial() const {
		return get_extension()->partial;
	}
};

class io_reply
{
	off_t offset: 40;
	unsigned int success: 1;
	unsigned int status: 16;

	unsigned int access_method: 1;
	unsigned long buf_size: 15;
	unsigned long buf_addr: 48;

	void init(char *buf, io_interface *io, off_t off, ssize_t size, int success,
			int status, int access_method) {
		long addr = (long) buf;
		this->buf_addr = addr;
		this->offset = off;
		this->buf_size = size;
		this->success = success;
		this->status = status;
		this->access_method = access_method;
	}
public:
	io_reply() {
		init(NULL, NULL, 0, 0, 0, 0, READ);
	}

	io_reply(io_request *req, int success, int status) {
		init(req->get_buf(), req->get_io(), req->get_offset(), req->get_size(),
					success, status, req->get_access_method());
	}

	int get_status() const {
		return status;
	}

	bool is_success() const {
		return success;
	}

	char *get_buf() const {
		long addr = buf_addr;
		return (char *) addr;
	}

	off_t get_offset() const {
		return offset;
	}

	ssize_t get_size() const {
		return buf_size;
	}

	int get_access_method() const {
		return access_method;
	}
};

template<class T>
class msg_sender
{
	T *buf;
	int buf_size;		// the max number of messages that can be buffered
	int num_current;	// the current number of messages in the buffer.
	fifo_queue<T> **dest_queues;
	int num_queues;
public:
	/**
	 * buf_size: the number of messages that can be buffered in the sender.
	 */
	msg_sender(int buf_size, fifo_queue<T> **queues,
			int num_queues);

	~msg_sender() {
		numa_free(buf, sizeof(T) * buf_size);
		numa_free(dest_queues, sizeof(fifo_queue<T> *) * num_queues);
	}

	int num_msg() {
		return num_current;
	}

	int flush();

	void flush_all() {
		while (num_msg() > 0)
			flush();
	}

	int send_cached(T *msg);
};

template<class T>
class thread_safe_msg_sender
{
	thread_safe_FIFO_queue<T> buf;
	std::vector<fifo_queue<T> *> dest_queues;
public:
	/**
	 * buf_size: the number of messages that can be buffered in the sender.
	 */
	thread_safe_msg_sender(int buf_size, fifo_queue<T> **queues,
			int num_queues): buf(buf_size), dest_queues(num_queues) {
		for (int i = 0; i < num_queues; i++)
			dest_queues[i] = queues[i];
	}

	int num_msg() {
		return buf.get_num_entries();
	}

	int flush();

	void flush_all() {
		while (num_msg() > 0)
			flush();
	}

	int send_cached(T *msg);
	int send_cached(T *msg, int num);
	int send(T *msg, int num);
};

template<class T>
class simple_msg_sender
{
	fifo_queue<T> buf;
	blocking_FIFO_queue<T> *queue;
public:
	/**
	 * buf_size: the number of messages that can be buffered in the sender.
	 */
	simple_msg_sender(blocking_FIFO_queue<T> *queue,
			int init_queue_size): buf(init_queue_size, true) {
		this->queue = queue;
	}

	int flush(bool blocking) {
		if (buf.is_empty()) {
			return 0;
		}
		if (blocking)
			return queue->add_partial(&buf);
		else
			return queue->non_blocking_add(&buf);
	}

	void flush_all() {
		while (!buf.is_empty())
			queue->add(&buf);
	}

	int get_num_remaining() {
		return buf.get_num_entries();
	}

	int send_cached(T *msg) {
		if (buf.is_full())
			buf.expand_queue(buf.get_size() * 2);
		buf.push_back(*msg);
		return 1;
	}

	int send_cached(T *msgs, int num) {
		int orig_num = num;
		int ret = buf.add(msgs, num);
		while (ret < num) {
			buf.expand_queue(buf.get_size() * 2);
			msgs += ret;
			num -= ret;
			ret = buf.add(msgs, num);
		}
		return orig_num;
	}

	blocking_FIFO_queue<T> *get_queue() const {
		return queue;
	}
};

class request_sender: public simple_msg_sender<io_request>
{
public:
	request_sender(blocking_FIFO_queue<io_request> *queue,
			int init_queue_size): simple_msg_sender(queue, init_queue_size) {
	}
};

#endif
