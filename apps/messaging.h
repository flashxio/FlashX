#ifndef __GRAPH_MESSAGING_H__
#define __GRAPH_MESSAGING_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SA-GraphLib.
 *
 * SA-GraphLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SA-GraphLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SA-GraphLib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "slab_allocator.h"

#include "vertex.h"
#include "partitioner.h"

class message
{
	// The allocator of the message buffer.
	slab_allocator *alloc;

	char *buf;
	int curr_get_off;
	int curr_add_off;
	int num_objs;

	void init() {
		alloc = NULL;
		buf = NULL;
		curr_get_off = 0;
		curr_add_off = 0;
		num_objs = 0;
	}

	void destroy() {
		if (buf) {
			alloc->free(buf);
		}
	}

public:
	message() {
		init();
	}

	message(slab_allocator *alloc) {
		init();
		this->alloc = alloc;
		this->buf = alloc->alloc();
		assert(this->buf);
	}

	~message() {
		destroy();
	}

	message(message &msg) {
		memcpy(this, &msg, sizeof(msg));
		msg.init();
	}

	message &operator=(message &msg) {
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

	int size() const {
		return alloc->get_obj_size();
	}

	bool has_next() const {
		return num_objs > 0;
	}

	template<class T>
	int get_next(T *objs[], int num) {
		int i;
		for (i = 0; has_next() && i < num; i++) {
			int remaining = size() - curr_get_off;
			objs[i] = T::deserialize(&buf[curr_get_off], remaining);
			curr_get_off += objs[i]->get_serialized_size();
			num_objs--;
		}
		return i;
	}

	template<class T>
	T *add(T &obj) {
		int remaining = size() - curr_add_off;
		if (remaining < obj.get_serialized_size())
			return NULL;
		T *obj_p = (T *) &buf[curr_add_off];
		curr_add_off += obj.serialize(&buf[curr_add_off], remaining);
		num_objs++;
		return obj_p;
	}

	int inc_msg_size(int msg_size) {
		int remaining = size() - curr_add_off;
		if (remaining < msg_size)
			return 0;
		curr_add_off += msg_size;
		return msg_size;
	}

	template<class T>
	int add(T *objs, int num = 1) {
		int num_added = 0;
		for (int i = 0; i < num; i++) {
			int remaining = size() - curr_add_off;
			if (remaining < objs[i].get_serialized_size())
				return num_added;
			curr_add_off += objs[i].serialize(&buf[curr_add_off], remaining);
			num_objs++;
			num_added++;
		}
		return num_added;
	}

	template<class T>
	void remove_last_obj(T &obj) {
		assert(num_objs > 0);
		assert(curr_add_off >= obj.get_serialized_size());
		T *obj_p = (T *) &buf[curr_add_off - obj.get_serialized_size()];
		assert(obj_p->get_serialized_size() == obj.get_serialized_size());
		num_objs--;
		curr_add_off -= obj.get_serialized_size();
	}

	bool copy_to(message &msg) {
		assert(msg.alloc);
		assert(msg.size() >= this->size());
		memcpy(msg.buf, this->buf, curr_add_off);
		// It probably makes more sense to reset the get offset.
		// I expect the user will iterate all objects the message later.
		msg.curr_get_off = 0;
		msg.curr_add_off = this->curr_add_off;
		msg.num_objs = this->num_objs;
		// After we copy all objects to another message, the current
		// message doesn't contain objects.
		this->num_objs = 0;
		this->curr_get_off = this->curr_add_off;
		return true;
	}
};

class msg_queue: public thread_safe_FIFO_queue<message>
{
public:
	msg_queue(int node_id, const std::string _name, int init_size,
			int max_size): thread_safe_FIFO_queue<message>(_name,
				node_id, init_size, max_size) {
	}

	static msg_queue *create(int node_id, const std::string name,
			int init_size, int max_size) {
		return new msg_queue(node_id, name, init_size, max_size);
	}

	static void destroy(msg_queue *q) {
		delete q;
	}

	/**
	 * This method needs to be used with caution.
	 * It may change the behavior of other threads if they also access
	 * the queue, so it's better to use it when no other threads are
	 * using it.
	 * It is also a heavy operation.
	 */
	int get_num_objs() {
		int num = thread_safe_FIFO_queue<message>::get_num_entries();
		stack_array<message> msgs(num);
		int ret = thread_safe_FIFO_queue<message>::fetch(msgs.data(), num);
		int num_objs = 0;
		for (int i = 0; i < ret; i++) {
			num_objs += msgs[i].get_num_objs();
		}
		int tmp = thread_safe_FIFO_queue<message>::add(msgs.data(), ret);
		assert(ret == tmp);
		return num_objs;
	}
};

class simple_msg_sender
{
	slab_allocator *alloc;
	message buf;
	msg_queue *queue;
	int num_objs;

protected:
	/**
	 * buf_size: the number of messages that can be buffered in the sender.
	 */
	simple_msg_sender(int node_id, slab_allocator *alloc,
			msg_queue *queue): buf(alloc) {
		this->alloc = alloc;
		this->queue = queue;
		num_objs = 0;
	}

public:
	static simple_msg_sender *create(int node_id, slab_allocator *alloc,
			msg_queue *queue) {
		assert(node_id >= 0);
		return new simple_msg_sender(node_id, alloc, queue);
	}

	static void destroy(simple_msg_sender *s) {
		delete s;
	}

	int flush() {
		num_objs = 0;
		if (buf.is_empty()) {
			return 0;
		}
		queue->add(&buf, 1);
		// We have to make sure all messages have been sent to the queue.
		assert(buf.is_empty());
		message tmp(alloc);
		buf = tmp;
		return 1;
	}

	/**
	 * This returns the number of remaining messages instead of the number
	 * of remaining objects.
	 */
	int get_num_remaining() {
		return num_objs;
	}

	template<class T>
	int send_cached(T &msg) {
		num_objs++;
		T *ret = buf.add(msg);
		if (ret == NULL) {
			flush();
			ret = buf.add(msg);
			assert(ret != NULL);
		}
		return 1;
	}

	msg_queue *get_queue() const {
		return queue;
	}
};

/**
 * The following part is more vertex-specific message passing.
 */

class vertex_message
{
protected:
	unsigned multicast: 1;
	// This is a flag to indicate that a user also wants to activate
	// the destination vertex.
	unsigned activate: 1;
	// This is a flag to indicate that it's a system's activation message.
	unsigned activation_msg: 1;
	unsigned size: 29;
	union {
		vertex_id_t dest;
		int num_dests;
	} u;
public:
	static vertex_message *deserialize(char *buf, int size) {
		vertex_message *msg = (vertex_message *) buf;
		assert(msg->size <= size);
		return msg;
	}

	vertex_message(int size, bool activate) {
		this->activate = activate;
		this->multicast = 0;
		this->activation_msg = 0;
		this->u.dest = -1;
		this->size = size;
		assert(size % 4 == 0);
	}

	void set_dest(local_vid_t id) {
		this->u.dest = id.id;
	}

	bool is_activation_msg() const {
		return activation_msg;
	}

	bool is_activate() const {
		return activate;
	}

	bool is_multicast() const {
		return multicast;
	}

	local_vid_t get_dest() const {
		return local_vid_t(u.dest);
	}

	int get_serialized_size() const {
		return size;
	}

	bool is_empty() const {
		return size == sizeof(vertex_message);
	}

	int serialize(char *buf, int size) const {
		assert(this->size <= size);
		memcpy(buf, this, this->size);
		return this->size;
	}
};

class multicast_message;

class multicast_dest_list
{
	multicast_message *msg;
	vertex_id_t *dest_list;
public:
	multicast_dest_list() {
		msg = NULL;
		dest_list = NULL;
	}

	void clear() {
		msg = NULL;
		dest_list = NULL;
	}

	multicast_dest_list(multicast_message *msg);
	void add_dest(local_vid_t id);
	int get_num_dests() const;
	local_vid_t get_dest(int idx) const;
};

class multicast_message: public vertex_message
{
	const vertex_id_t *get_dest_begin() const {
		return (vertex_id_t *) (((char *) this) + get_orig_msg_size());
	}

	vertex_id_t *get_dest_begin() {
		return (vertex_id_t *) (((char *) this) + get_orig_msg_size());
	}

	int get_orig_msg_size() const {
		return size - u.num_dests * sizeof(vertex_id_t);
	}
public:
	static multicast_message *convert2multicast(vertex_message *msg) {
		multicast_message *mmsg = (multicast_message *) msg;
		mmsg->u.num_dests = 0;
		mmsg->multicast = 1;
		return mmsg;
	}

	static multicast_message *cast2multicast(vertex_message *msg) {
		multicast_message *mmsg = (multicast_message *) msg;
		assert(mmsg->multicast);
		return mmsg;
	}

	multicast_message(int size, bool activate): vertex_message(size,
			activate) {
		this->multicast = 1;
		this->u.num_dests = 0;
	}

	bool is_empty() const {
		return (size_t) get_orig_msg_size() == sizeof(vertex_message);
	}

	bool is_multicast() const {
		return multicast;
	}

	int get_num_dests() const {
		return u.num_dests;
	}

	multicast_dest_list get_dest_list() {
		return multicast_dest_list(this);
	}

	int get_serialized_size() const {
		return size;
	}

	/**
	 * This is the body size of a multicast message.
	 * The destination list is excluded, but it includes the header
	 * of vertex_message.
	 */
	int get_body_size() const {
		return get_serialized_size() - get_num_dests() * sizeof(vertex_id_t);
	}

	friend class multicast_dest_list;
};

class activation_message: public multicast_message
{
public:
	activation_message(): multicast_message(sizeof(vertex_message), true) {
		this->activation_msg = 1;
	}
};

inline multicast_dest_list::multicast_dest_list(multicast_message *msg)
{
	this->msg = msg;
	dest_list = msg->get_dest_begin();
}

inline void multicast_dest_list::add_dest(local_vid_t id)
{
	dest_list[msg->u.num_dests++] = id.id;
	msg->size += sizeof(id.id);
}

inline int multicast_dest_list::get_num_dests() const
{
	return msg->u.num_dests;
}

inline local_vid_t multicast_dest_list::get_dest(int idx) const
{
	return local_vid_t(dest_list[idx]);
}

class multicast_msg_sender
{
	const static int MMSG_BUF_SIZE = PAGE_SIZE;
	slab_allocator *alloc;
	// The local buffer of vertex messages.
	message buf;
	// The destination queue of the sender.
	msg_queue *queue;
	// The initialized multicast vertex message.
	multicast_message *mmsg;
	// The number of destinations in the multicast message.
	int num_dests;
	// This stores the body of the multicast message.
	char mmsg_temp_buf[MMSG_BUF_SIZE];
	multicast_dest_list dest_list;

	multicast_msg_sender(slab_allocator *alloc,
			msg_queue *queue): buf(alloc) {
		this->alloc = alloc;
		this->queue = queue;
		this->mmsg = NULL;
		this->num_dests = 0;
	}
public:
	static multicast_msg_sender *create(slab_allocator *alloc,
			msg_queue *queue) {
		return new multicast_msg_sender(alloc, queue);
	}

	static void destroy(multicast_msg_sender *s) {
		delete s;
	}

	int flush();

	template<class T>
	void init(const T &msg) {
		assert(mmsg == NULL);
		num_dests = 0;
		vertex_message *p = (vertex_message *) buf.add(msg);
		if (p == NULL) {
			flush();
			p = (vertex_message *) buf.add(msg);
			assert(p);
		}
		this->mmsg = multicast_message::convert2multicast(p);
		assert(this->mmsg->get_serialized_size() <= MMSG_BUF_SIZE);
		this->mmsg->serialize(mmsg_temp_buf, MMSG_BUF_SIZE);
		dest_list = this->mmsg->get_dest_list();
	}

	bool add_dest(local_vid_t id);

	void end_multicast() {
		if (num_dests == 0) {
			multicast_message *mmsg_template
				= (multicast_message *) mmsg_temp_buf;
			buf.remove_last_obj(*mmsg_template);
		}
		mmsg = NULL;
		num_dests = 0;
		dest_list.clear();
	}
};

#endif
