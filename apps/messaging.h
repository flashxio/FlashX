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
	int add(T &obj) {
		int remaining = size() - curr_add_off;
		if (remaining < obj.get_serialized_size())
			return 0;
		curr_add_off += obj.serialize(&buf[curr_add_off], remaining);
		num_objs++;
		return 1;
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
		void *addr;
		if (node_id < 0)
			addr = numa_alloc_local(sizeof(msg_queue));
		else
			addr = numa_alloc_onnode(sizeof(msg_queue), node_id);
		return new(addr) msg_queue(node_id, name, init_size, max_size);
	}

	static void destroy(msg_queue *q) {
		q->~msg_queue();
		numa_free(q, sizeof(*q));
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
		void *addr = numa_alloc_onnode(sizeof(simple_msg_sender), node_id);
		return new(addr) simple_msg_sender(node_id, alloc, queue);
	}

	static void destroy(simple_msg_sender *s) {
		s->~simple_msg_sender();
		numa_free(s, sizeof(*s));
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
		int ret = buf.add(msg);
		if (ret == 0) {
			flush();
			ret = buf.add(msg);
			assert(ret == 1);
		}
		return 1;
	}

	msg_queue *get_queue() const {
		return queue;
	}
};

#endif
