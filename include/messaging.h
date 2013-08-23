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

const int MAX_INLINE_SIZE=128;

class user_compute
{
public:
	virtual int serialize(char *buf, int size) const = 0;
	virtual int get_serialized_size() const = 0;
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
	enum {
		BASIC_REQ,
		EXT_REQ,
		USER_COMPUTE,
	};

	off_t offset: 40;
	static const int MAX_BUF_SIZE = (1 << 24) - 1;
	unsigned long buf_size: 24;

	// These two flags decide how the payload is interpreted, so they are
	// initialized when the object is created and can't be changed manually.
	// The only case that payload_type is changed is when
	// another request object is copied to this request object.
	// What needs to be guaranteed is that the pointer in the payload always
	// points to the object of the right type.
	// However, data_inline is never changed. It is usually false. It is
	// used when to interpret the io requests in a message.
	unsigned int payload_type: 2;
	unsigned int data_inline: 1;

	unsigned int access_method: 1;
	// Is this synchronous IO?
	unsigned int sync: 1;
	unsigned int high_prio: 1;
	static const int MAX_NODE_ID = (1 << 10) - 1;
	unsigned int node_id: 10;
	// Linux uses 48 bit for addresses.
	unsigned long io_addr: 48;

	union {
		void *buf_addr;
		user_compute *compute;
		io_req_extension *ext;
		char buf[0];
	} payload;

	io_req_extension *get_extension() const {
		assert(is_extended_req() && payload.ext);
		if (data_inline)
			return (io_req_extension *) payload.buf;
		else
			return payload.ext;
	}

	bool use_embedded() const {
		return get_extension()->vec_pointer == get_extension()->embedded_vecs;
	}

public:
	// By default, a request is initialized as a flush request.
	io_request() {
		payload_type = BASIC_REQ;
		payload.buf_addr = NULL;
		high_prio = 1;
		sync = 1;
		data_inline = 0;
	}

	io_request(bool extended) {
		if (extended) {
			payload_type = EXT_REQ;
			payload.ext = new io_req_extension();
		}
		else {
			payload_type = BASIC_REQ;
			payload.buf_addr = NULL;
		}
		high_prio = 1;
		sync = false;
		data_inline = 0;
	}

	io_request(char *buf, off_t off, ssize_t size, int access_method,
			io_interface *io, int node_id, bool sync = false) {
		payload_type = BASIC_REQ;
		data_inline = 0;
		init(buf, off, size, access_method, io, node_id, sync);
	}

	io_request(off_t off, int access_method, io_interface *io, int node_id,
			io_request *orig, void *priv, bool sync = false) {
		payload_type = EXT_REQ;
		data_inline = 0;
		payload.ext = new io_req_extension();
		init(off, access_method, io, node_id, orig, priv, NULL, sync);
	}

	io_request(user_compute *compute, off_t off, ssize_t size,
			int access_method, io_interface *io, int node_id,
			bool sync = false) {
		payload_type = USER_COMPUTE;
		data_inline = 0;
		init(NULL, off, size, access_method, io, node_id, sync);
		payload.compute = compute;
	}

	io_request(char *buf, off_t off, ssize_t size, int access_method,
			io_interface *io, int node_id, io_request *orig,
			void *priv, bool sync = false) {
		payload_type = EXT_REQ;
		data_inline = 0;
		payload.ext = new io_req_extension();
		init(buf, off, size, access_method, io, node_id, orig, priv, NULL, sync);
	}

	io_request(io_request &req) {
		assert(!req.data_inline);
		memcpy(this, &req, sizeof(req));
		if (req.is_extended_req()) {
			req.payload_type = BASIC_REQ;
			req.payload.buf_addr = NULL;
		}
	}

	io_request &operator=(io_request &req) {
		assert(!data_inline && !req.data_inline);
		// We need to free its own extension first.
		if (this->is_extended_req())
			delete this->get_extension();
		memcpy(this, &req, sizeof(req));
		if (req.is_extended_req()) {
			req.payload_type = BASIC_REQ;
			req.payload.buf_addr = NULL;
		}
		return *this;
	}

	~io_request() {
		if (is_extended_req() && get_extension())
			delete get_extension();
	}

	void init(const io_request &req) {
		assert(!data_inline);
		this->sync = req.sync;
		if (req.payload_type == USER_COMPUTE
				|| this->payload_type == USER_COMPUTE) {
			assert(req.payload_type == USER_COMPUTE
					&& this->payload_type == USER_COMPUTE);
			// TODO
		}
		else if (!req.is_extended_req()) {
			this->init(req.get_buf(), req.get_offset(), req.get_size(),
					req.get_access_method(), req.get_io(), req.get_node_id());
		}
		else if (this->is_extended_req()) {
			this->init(req.get_offset(), req.get_access_method(), req.get_io(),
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
		data_inline = 0;
		if (is_extended_req()) {
			io_req_extension *ext = get_extension();
			assert(ext);
			ext->init();
			offset = 0;
			high_prio = 0;
			sync = 0;
			node_id = 0;
			io_addr = 0;
			access_method = 0;
			buf_size = 0;
		}
		else {
			offset = 0;
			payload_type = BASIC_REQ;
			high_prio = 0;
			sync = 0;
			node_id = 0;
			io_addr = 0;
			access_method = 0;
			buf_size = 0;
			payload.buf_addr = NULL;
		}
	}

	void init(char *buf, off_t off, ssize_t size, int access_method,
			io_interface *io, int node_id, int sync = false);

	void init(char *buf, off_t off, ssize_t size, int access_method,
			io_interface *io, int node_id, io_request *orig, void *priv,
			void *user_data, int sync = false) {
		init(off, access_method, io, node_id, orig, priv, user_data, sync);
		add_buf(buf, size);
	}

	void init(off_t off, int access_method, io_interface *io, int node_id,
			io_request *orig, void *priv, void *user_data, int sync = false) {
		assert(is_extended_req());
		io_request::init(NULL, off, 0, access_method, io, node_id, sync);
		io_req_extension *ext = get_extension();
		ext->priv = priv;
		ext->orig = orig;
		ext->user_data = user_data;
	}

	int get_file_id() const;

	/**
	 * Test whether the request is a flush request.
	 * The difference of a sync'd request and a flush request is that
	 * a flush request isn't a valid request for accessing data.
	 */
	bool is_flush() const {
		return sync && high_prio && (payload.buf_addr == NULL);
	}

	bool is_sync() const {
		return sync;
	}

	bool is_extended_req() const {
		return payload_type == EXT_REQ;
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

	void set_io(io_interface *io) {
		this->io_addr = (long) io;
	}

	io_interface *get_io() const {
		return (io_interface *) (long) io_addr;
	}

	int get_node_id() const {
		return node_id;
	}

	void set_node_id(int node_id) {
		assert(node_id <= MAX_NODE_ID);
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
			return (char *) payload.buf_addr;
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
class thread_safe_msg_sender
{
	thread_safe_FIFO_queue<T> buf;
	fifo_queue<T> *dest_queue;

	/**
	 * buf_size: the number of messages that can be buffered in the sender.
	 */
	thread_safe_msg_sender(int node_id, int buf_size, fifo_queue<T> *queue): buf(
			node_id, buf_size) {
		dest_queue = queue;
	}

public:
	static thread_safe_msg_sender<T> *create(int node_id, int buf_size,
			fifo_queue<T> *queue) {
		assert(node_id >= 0);
		void *addr = numa_alloc_onnode(sizeof(thread_safe_msg_sender<T>),
				node_id);
		return new(addr) thread_safe_msg_sender<T>(node_id, buf_size, queue);
	}

	static void destroy(thread_safe_msg_sender<T> *s) {
		s->~thread_safe_msg_sender();
		numa_free(s, sizeof(*s));
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

protected:
	/**
	 * buf_size: the number of messages that can be buffered in the sender.
	 */
	simple_msg_sender(int node_id, blocking_FIFO_queue<T> *queue,
			int init_queue_size): buf(node_id, init_queue_size, true) {
		this->queue = queue;
	}

public:
	static simple_msg_sender<T> *create(int node_id, blocking_FIFO_queue<T> *queue,
			int init_queue_size) {
		assert(node_id >= 0);
		void *addr = numa_alloc_onnode(sizeof(simple_msg_sender<T>), node_id);
		return new(addr) simple_msg_sender<T>(queue, init_queue_size);
	}

	static void destroy(simple_msg_sender<T> *s) {
		s->~simple_msg_sender();
		numa_free(s, sizeof(*s));
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
	request_sender(int node_id, blocking_FIFO_queue<io_request> *queue,
			int init_queue_size): simple_msg_sender(node_id, queue,
				init_queue_size) {
	}

public:
	static request_sender *create(int node_id,
			blocking_FIFO_queue<io_request> *queue, int init_queue_size) {
		assert(node_id >= 0);
		void *addr = numa_alloc_onnode(sizeof(request_sender), node_id);
		return new(addr) request_sender(node_id, queue, init_queue_size);
	}

	static void destroy(request_sender *s) {
		s->~request_sender();
		numa_free(s, sizeof(*s));
	}
};

#endif
