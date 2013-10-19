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
#include "io_request.h"

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
			node_id, INIT_MSG_BUF_SIZE, true) {
		this->alloc = alloc;
		this->accept_inline = accept_inline;
	}

	int add_objs(T *objs, int num = 1) {
		int num_added = 0;
		if (fifo_queue<message<T> >::is_empty()) {
			message<T> tmp(alloc, accept_inline);
			add_msg(tmp);
		}
		while (num > 0) {
			int ret = fifo_queue<message<T> >::back().add(objs, num);
			// The last message is full. We need to add a new message
			// to the queue.
			if (ret == 0) {
				message<T> tmp(alloc, accept_inline);
				add_msg(tmp);
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
class msg_queue: public thread_safe_FIFO_queue<message<T> >
{
	// TODO I may need to make sure all messages are compatible with the flag.
	bool accept_inline;
	const std::string name;

public:
	msg_queue(int node_id, const std::string _name, int init_size, int max_size,
			bool accept_inline): thread_safe_FIFO_queue<message<T> >(node_id,
				init_size, max_size), name(_name) {
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

	/**
	 * This method needs to be used with caution.
	 * It may change the behavior of other threads if they also access
	 * the queue, so it's better to use it when no other threads are
	 * using it.
	 * It is also a heavy operation.
	 */
	int get_num_objs() {
		int num = thread_safe_FIFO_queue<message<T> >::get_num_entries();
#ifdef MEMCHECK
		message<T> *msgs = new message<T>[num];
#else
		message<T> msgs[num];
#endif
		int ret = thread_safe_FIFO_queue<message<T> >::fetch(msgs, num);
		int num_objs = 0;
		for (int i = 0; i < ret; i++) {
			num_objs += msgs[i].get_num_objs();
		}
		thread_safe_FIFO_queue<message<T> >::add(msgs, ret);
#ifdef MEMCHECK
		delete [] msgs;
#endif
		return num_objs;
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

	int get_num_remaining() const {
		return buf.get_num_objs();
	}
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

	int flush() {
		if (buf.is_empty()) {
			return 0;
		}
		queue->add(&buf);
		// We have to make sure all messages have been sent to the queue.
		assert(buf.is_empty());
		return 1;
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

#endif
