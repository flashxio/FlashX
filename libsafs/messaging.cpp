/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>

#include <algorithm>

#include "messaging.h"
#include "io_interface.h"
#include "slab_allocator.h"
#include "cache.h"

template<class T>
message<T>::message(slab_allocator *alloc, bool accept_inline)
{
	init();
	this->alloc = alloc;
	this->buf = alloc->alloc();
	assert(this->buf);
	this->accept_inline = accept_inline;
}

template<class T>
void message<T>::destroy()
{
	if (buf) {
		// we need to destroy the remaining objects in the message buffer
		// only when objects in the message aren't inline.
		if (!accept_inline) {
			while (has_next()) {
				T *obj = get_next_addr();
				obj->~T();
			}
		}
		alloc->free(buf);
	}
}

template<class T>
int message<T>::size() const
{
	return alloc->get_obj_size();
}

template<class T>
int thread_safe_msg_sender<T>::send_cached(T *msg, int num)
{
	int num_added = 0;
	// We expect the method is always successful.
	// so we try again and again until we succeed.
	while (true) {
		pthread_spin_lock(&_lock);
		int ret = buf.add(msg, num);
		pthread_spin_unlock(&_lock);
		msg += ret;
		num -= ret;
		num_added += ret;
		if (num == 0)
			return num_added;
		// If the buffer is full, we should flush the buffer.
		flush();
	}
}

// Send msgs to the destinatiion queue directly without caching.
template<class T>
int thread_safe_msg_sender<T>::send(T *msg, int num)
{
	// We should flush the msgs in the cache first.
	// but it won't flush all msgs in the cache.
	flush();

	int num_sent = 0;
	while (num > 0) {
		message<T> tmp(alloc, dest_queue->is_accept_inline());
		int ret = tmp.add(msg, num);
		msg += ret;
		num -= ret;
		num_sent += ret;

		// We need to make sure the message is added to the queue.
		// TODO we should add multiple messages together.
		while (dest_queue->add(&tmp, 1) <= 0);
	}

	return num_sent;
}

/**
 * these are to force to instantiate the templates
 * for io_request and io_reply.
 */
template class thread_safe_FIFO_queue<message<io_request> >;
template class thread_safe_FIFO_queue<message<io_reply> >;
template class blocking_FIFO_queue<message<io_request> >;
template class blocking_FIFO_queue<message<io_reply> >;
template class thread_safe_msg_sender<io_reply>;
template class message<io_request>;
template class message<io_reply>;

atomic_unsigned_integer io_req_extension::num_creates;
