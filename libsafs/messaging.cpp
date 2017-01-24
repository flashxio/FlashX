/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdio.h>

#include <algorithm>

#include "messaging.h"
#include "io_interface.h"
#include "slab_allocator.h"
#include "cache.h"

namespace safs
{

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
		_lock.lock();
		int ret = buf.add(msg, num);
		_lock.unlock();
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

template class thread_safe_msg_sender<safs::io_reply>;

}

/**
 * these are to force to instantiate the templates
 * for io_request and io_reply.
 */
template class thread_safe_FIFO_queue<safs::message<safs::io_request> >;
template class thread_safe_FIFO_queue<safs::message<safs::io_reply> >;
template class blocking_FIFO_queue<safs::message<safs::io_request> >;
template class blocking_FIFO_queue<safs::message<safs::io_reply> >;
template class safs::message<safs::io_request>;
template class safs::message<safs::io_reply>;
