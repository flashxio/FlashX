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
		// TODO we should remove the ownership to the request extension
		// from the IO request. i.e., the extension should be allocated and
		// deallocated by someone else.
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

	bool is_data_inline() const {
		return data_inline == 1;
	}

	/**
	 * We need to serialize an io request to a buffer so it can be sent to
	 * another thread.
	 * @accept_inline: indicates whether the IO request can inline data
	 * to the buffer.
	 */
	int serialize(char *buf, int size, bool accept_inline) {
		int serialized_size;
		if (is_data_inline()) {
			assert(accept_inline);
			serialized_size = get_serialized_size();
			assert(serialized_size <= size);
			memcpy(buf, this, serialized_size);
		}
		else if (payload_type == EXT_REQ) {
			// We never serialize the io request extension to the message.
			serialized_size = sizeof(io_request);
			assert(serialized_size <= size);
			memcpy(buf, this, sizeof(*this));
			payload.ext = NULL;
			payload_type = BASIC_REQ;
		}
		else if (payload_type == BASIC_REQ) {
			// We only serialize the data buffer to the message for write
			// requests. The size of the data buffer has to be small.
			if (get_size() > MAX_INLINE_SIZE || access_method == READ
					|| !accept_inline) {
				serialized_size = sizeof(io_request);
				assert(serialized_size <= size);
				memcpy(buf, this, sizeof(*this));
			}
			else {
				serialized_size = sizeof(io_request) - sizeof(this->payload)
					+ get_size();
				assert(serialized_size <= size);
				memcpy(buf, this, sizeof(*this));
				io_request *p = (io_request *) buf;
				p->data_inline = 1;
				memcpy(p->payload.buf, (char *) this->payload.buf_addr,
						this->get_size());
			}
		}
		else {
			assert(accept_inline);
			// The user compute object is always serialized to the message.
			user_compute *compute = this->payload.compute;
			serialized_size = sizeof(io_request) - sizeof(this->payload)
				+ compute->get_serialized_size();
			assert(serialized_size <= size && sizeof(io_request)
					<= (unsigned) size);
			memcpy(buf, this, sizeof(*this));
			io_request *p = (io_request *) buf;
			p->data_inline = 1;
			compute->serialize((char *) p->payload.buf,
					size - (sizeof(io_request) - sizeof(this->payload)));
		}
		return serialized_size;
	}

	/**
	 * This method returns the size of an IO request after it is serialized.
	 */
	int get_serialized_size() const {
		if (payload_type == EXT_REQ)
			return sizeof(io_request);
		else if (payload_type == BASIC_REQ && (get_size() > MAX_INLINE_SIZE
					|| access_method == READ))
			return sizeof(io_request);
		else if (payload_type == BASIC_REQ)
			return sizeof(io_request) - sizeof(this->payload) + get_size();
		else if (!is_data_inline()) {
			user_compute *compute = this->payload.compute;
			return sizeof(io_request) - sizeof(this->payload)
				+ compute->get_serialized_size();
		}
		else {
			user_compute *compute = (user_compute *) this->payload.buf;
			return sizeof(io_request) - sizeof(this->payload)
				+ compute->get_serialized_size();
		}
	}

	/**
	 * This method deserialize an request from the buffer.
	 * If the request data is inline in the buffer, instead of allocating
	 * memory for the extra objects of the IO request, the extra objects
	 * will stay in the buffer and the created request will point to the buffer.
	 */
	static void deserialize(io_request &req, char *buf, int size) {
		assert((unsigned) size >= sizeof(io_request));
		io_request *p = (io_request *) buf;
		memcpy(&req, buf, sizeof(req));
		if (req.is_data_inline()) {
			assert(req.payload_type != EXT_REQ);
			switch(req.payload_type) {
				case BASIC_REQ:
					req.payload.buf_addr = p->payload.buf;
					break;
				case USER_COMPUTE:
					req.payload.compute = (user_compute *) p->payload.buf;
					break;
				default:
					assert(0);
			}
			req.data_inline = 0;
		}
	}

	static io_request *deserialize(char *buf, int size) {
		assert((unsigned) size >= sizeof(io_request));
		io_request *ret = (io_request *) buf;
		assert(ret->get_serialized_size() <= size);
		return ret;
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

	bool is_data_inline() const {
		return false;
	}

	int serialize(char *buf, int size, bool accept_inline) {
		assert((unsigned) size >= sizeof(*this));
		memcpy(buf, this, sizeof(*this));
		return get_serialized_size();
	}

	int get_serialized_size() const {
		return sizeof(*this);
	}

	static void deserialize(io_reply &reply, char *buf, int size) {
		assert((unsigned) size >= sizeof(io_reply));
		reply = *(io_reply *) buf;
	}

	static io_reply *deserialize(char *buf, int size) {
		assert((unsigned) size >= sizeof(io_reply));
		io_reply *ret = (io_reply *) buf;
		assert(ret->get_serialized_size() <= size);
		return ret;
	}
};

class slab_allocator;

/**
 * It is an object container used for message passing.
 * Instead of sending objects to another thread directly, we add objects
 * to a message, and send the message to the thread.
 *
 * The message supports objects of variant sizes, but objects need to
 * be able to serialize and deserialize itself from the message.
 *
 * The objects that can be added to a message need to support
 *	serialize(),
 *	get_serialized_size(),
 *	deserialize().
 * If the message allows objects inline in the message buffer, after objects
 * are serialized in the message, they shouldn't own any memory, and
 * the message isn't responsible for destroying them when the message
 * is destroyed.
 * If the message doesn't allow objects inline in the message buffer,
 * the objects stored in the message may still own their own memory. Thus,
 * when the message is destroyed, all objects will be destroyed.
 *
 * The ways of fetching objects are different in these two modes:
 * If the message accepts inline objects, get_next_inline() should be used.
 * Otherwise, get_next() should be used.
 *
 * There are multiple benefits of using the message object:
 * It reduces lock contention.
 * When it works with message senders, we can reduce the number of memory
 * copies. There are at most two memory copies: add an object to
 * the message; (optionally) the receiver thread copies the message to
 * the local memory if it runs on a NUMA node different from the message sender.
 * Previously, we need to copy an object to the sender's local buffer,
 * copy the local buffer to the queue, copy objects in the queue to the
 * remote buffer.
 */
template<class T>
class message
{
	// The allocator of the message buffer.
	slab_allocator *alloc;

	char *buf;
	short curr_get_off;
	short curr_add_off;
	short num_objs;
	// Indicate whether the data of an object can be inline in the message.
	short accept_inline: 1;

	void init() {
		alloc = NULL;
		buf = NULL;
		curr_get_off = 0;
		curr_add_off = 0;
		num_objs = 0;
		accept_inline = 0;
	}

	void destroy();

	/**
	 * This just returns the address of the next object.
	 */
	T *get_next_addr() {
		int remaining = size() - curr_get_off;
		assert(num_objs > 0);
		num_objs--;
		T *ret = T::deserialize(&buf[curr_get_off], remaining);
		curr_get_off += ret->get_serialized_size();
		return ret;
	}
public:
	message() {
		init();
	}

	message(slab_allocator *alloc, bool accept_inline);

	~message() {
		destroy();
	}

	message(message<T> &msg) {
		memcpy(this, &msg, sizeof(msg));
		msg.init();
	}

	message<T> &operator=(message<T> &msg) {
		destroy();
		memcpy(this, &msg, sizeof(msg));
		msg.init();
		return *this;
	}

	void clear() {
		destroy();
		init();
	}

	/**
	 * It's actually the number of remaining objects in the message.
	 */
	int get_num_objs() const {
		return num_objs;
	}

	bool is_empty() const {
		return get_num_objs() == 0;
	}

	int size() const;

	bool has_next() const {
		return num_objs > 0;
	}

	int get_next_inline(T objs[], int num_objs) {
		/* If the message accepts inline objects, there are no ownership
		 * problems. The memory owned by the objects has been embedded
		 * in the message buffer.
		 */
		assert(accept_inline);
		int i = 0;
		while (has_next()) {
			int remaining = size() - curr_get_off;
			assert(num_objs > 0);
			num_objs--;
			T::deserialize(objs[i], &buf[curr_get_off], remaining);
			curr_get_off += objs[i].get_serialized_size();
			i++;
		}
		return i;
	}

	bool get_next(T &obj) {
		T *ret = get_next_addr();
		assert(!accept_inline && !ret->is_data_inline());
		// The copy constructor of the object will remove the ownership
		// of any memory pointed by the object, so we don't need to call
		// the deconstructor of the object.
		obj = *ret;
		return true;
	}

	int get_next_objs(T objs[], int num) {
		for (int i = 0; i < num; i++) {
			if (has_next())
				get_next(objs[i]);
			else
				return i;
		}
		return num;
	}

	int add(T *objs, int num = 1) {
		int num_added = 0;
		for (int i = 0; i < num; i++) {
			int remaining = size() - curr_add_off;
			if (remaining < objs[i].get_serialized_size())
				return num_added;
			curr_add_off += objs[i].serialize(&buf[curr_add_off], remaining,
					accept_inline);
			num_objs++;
			num_added++;
		}
		return num_added;
	}

	bool copy_to(message<T> &msg) {
		assert(msg.alloc);
		assert(msg.size() >= this->size());
		memcpy(msg.buf, this->buf, curr_add_off);
		// It probably makes more sense to reset the get offset.
		// I expect the user will iterate all objects the message later.
		msg.curr_get_off = 0;
		msg.curr_add_off = this->curr_add_off;
		msg.num_objs = this->num_objs;
		msg.accept_inline = this->accept_inline;
		// After we copy all objects to another message, the current
		// message doesn't contain objects.
		this->num_objs = 0;
		this->curr_get_off = this->curr_add_off;
		return true;
	}
};

/**
 * It contains multiple messages. It basically helps construct messages.
 */
template<class T>
class msg_buffer: public fifo_queue<message<T> >
{
	static const int INIT_MSG_BUF_SIZE = 16;

	message<T> curr;
	slab_allocator *alloc;
	bool accept_inline;

	void add_msg(message<T> &msg) {
		if (fifo_queue<message<T> >::is_full()) {
			fifo_queue<message<T> >::expand_queue(
					fifo_queue<message<T> >::get_size() * 2);
		}
		fifo_queue<message<T> >::add(&msg, 1);
	}

public:
	msg_buffer(int node_id, slab_allocator *alloc,
			bool accpet_inline): fifo_queue<message<T> >(
			node_id, INIT_MSG_BUF_SIZE, true), curr(alloc, accept_inline) {
		this->alloc = alloc;
		this->accept_inline = accept_inline;
	}

	int add_objs(T *objs, int num = 1) {
		int num_added = 0;
		while (num > 0) {
			int ret = curr.add(objs, num);
			// The current message is full. We need to add the current
			// message to the queue and create a new one.
			if (ret == 0) {
				add_msg(curr);
				message<T> tmp(alloc, accept_inline);
				curr = tmp;
			}
			else {
				num_added += ret;
				objs += ret;
				num -= ret;
			}
		}
		return num_added;
	}
};

template<class T>
class msg_queue: public blocking_FIFO_queue<message<T> >
{
	// TODO I may need to make sure all messages are compatible with the flag.
	bool accept_inline;

public:
	msg_queue(int node_id, const std::string name, int init_size, int max_size,
			bool accept_inline): blocking_FIFO_queue<message<T> >(node_id,
				name, init_size, max_size) {
		this->accept_inline = accept_inline;
	}

	static msg_queue<T> *create(int node_id, const std::string name,
			int init_size, int max_size, bool accept_inline) {
		void *addr;
		if (node_id < 0)
			addr = numa_alloc_local(sizeof(msg_queue<T>));
		else
			addr = numa_alloc_onnode(sizeof(msg_queue<T>), node_id);
		return new(addr) msg_queue<T>(node_id, name, init_size, max_size,
				accept_inline);
	}

	static void destroy(msg_queue<T> *q) {
		q->~msg_queue();
		numa_free(q, sizeof(*q));
	}

	bool is_accept_inline() const {
		return accept_inline;
	}
};

template<class T>
class thread_safe_msg_sender
{
	pthread_spinlock_t _lock;
	message<T> buf;

	slab_allocator *alloc;
	msg_queue<T> *dest_queue;

	/**
	 * buf_size: the number of messages that can be buffered in the sender.
	 */
	thread_safe_msg_sender(slab_allocator *alloc,
			msg_queue<T> *queue): buf(alloc, queue->is_accept_inline()) {
		this->alloc = alloc;
		dest_queue = queue;
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

public:
	static thread_safe_msg_sender<T> *create(int node_id, slab_allocator *alloc,
			msg_queue<T> *queue) {
		assert(node_id >= 0);
		void *addr = numa_alloc_onnode(sizeof(thread_safe_msg_sender<T>),
				node_id);
		return new(addr) thread_safe_msg_sender<T>(alloc, queue);
	}

	static void destroy(thread_safe_msg_sender<T> *s) {
		s->~thread_safe_msg_sender();
		numa_free(s, sizeof(*s));
	}

	/**
	 * flush the entries in the buffer to the queues.
	 * return the number of entries that have been flushed.
	 */
	int flush() {
		pthread_spin_lock(&_lock);
		if (!buf.is_empty()) {
			message<T> tmp = buf;
			message<T> tmp1(alloc, dest_queue->is_accept_inline());
			buf = tmp1;
			pthread_spin_unlock(&_lock);
			return dest_queue->add(&tmp, 1);
		}
		else {
			pthread_spin_unlock(&_lock);
			return 0;
		}
	}

	void flush_all() {
		// flush_all() now is the same as flush().
		flush();
	}

	int send_cached(T *msg, int num = 1);
	int send(T *msg, int num);
};

template<class T>
class simple_msg_sender
{
	slab_allocator *alloc;
	msg_buffer<T> buf;
	msg_queue<T> *queue;

protected:
	/**
	 * buf_size: the number of messages that can be buffered in the sender.
	 */
	simple_msg_sender(int node_id, slab_allocator *alloc,
			msg_queue<T> *queue): buf(node_id, alloc, queue->is_accept_inline()) {
		this->alloc = alloc;
		this->queue = queue;
	}

public:
	static simple_msg_sender<T> *create(int node_id, slab_allocator *alloc,
			msg_queue<T> *queue) {
		assert(node_id >= 0);
		void *addr = numa_alloc_onnode(sizeof(simple_msg_sender<T>), node_id);
		return new(addr) simple_msg_sender<T>(node_id, alloc, queue);
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
			queue->add(&buf);
		else
			queue->non_blocking_add(&buf);
		return 1;
	}

	void flush_all() {
		while (!buf.is_empty())
			queue->add(&buf);
	}

	/**
	 * This returns the number of remaining messages instead of the number
	 * of remaining objects.
	 */
	int get_num_remaining() {
		return buf.get_num_entries();
	}

	int send_cached(T *msgs, int num = 1) {
		return buf.add_objs(msgs, num);
	}

	msg_queue<T> *get_queue() const {
		return queue;
	}
};

class request_sender: public simple_msg_sender<io_request>
{
	request_sender(int node_id, slab_allocator *alloc,
			msg_queue<io_request> *queue): simple_msg_sender(
				node_id, alloc, queue) {
	}

public:
	static request_sender *create(int node_id, slab_allocator *alloc,
			msg_queue<io_request> *queue) {
		assert(node_id >= 0);
		void *addr = numa_alloc_onnode(sizeof(request_sender), node_id);
		return new(addr) request_sender(node_id, alloc, queue);
	}

	static void destroy(request_sender *s) {
		s->~request_sender();
		numa_free(s, sizeof(*s));
	}
};

template<class T>
class simple_sender
{
	fifo_queue<T> buf;
	blocking_FIFO_queue<T> *queue;

protected:
	/**
	 *      * buf_size: the number of messages that can be buffered in the sender.
	 *           */
	simple_sender(int node_id, blocking_FIFO_queue<T> *queue,
			int init_queue_size): buf(node_id, init_queue_size, true) {
		this->queue = queue;
	}

public:
	static simple_sender<T> *create(int node_id, blocking_FIFO_queue<T> *queue,
			int init_queue_size) {
		assert(node_id >= 0);
		void *addr = numa_alloc_onnode(sizeof(simple_sender<T>), node_id);
		return new(addr) simple_sender<T>(queue, init_queue_size);
	}

	static void destroy(simple_sender<T> *s) {
		s->~simple_sender();
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

#endif
