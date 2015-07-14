/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
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
#include <malloc.h>

#include "local_mem_buffer.h"
#include "matrix_config.h"
#include "local_matrix_store.h"

namespace fm
{

namespace detail
{

local_mem_buffer::~local_mem_buffer()
{
	size_t num_bufs = 0;
	for (auto it = bufs.begin(); it != bufs.end(); it++)
		num_bufs += it->second.size();
	assert(num_bufs + num_frees == num_allocs);
}

namespace
{

class local_deleter
{
	// TODO maybe I need to use a smart pointer for this.
	std::deque<char *> &q;
public:
	local_deleter(std::deque<char *> &_q): q(_q) {
	}

	void operator()(char *addr) {
		q.push_back(addr);
	}
};

}

std::shared_ptr<char> local_mem_buffer::_alloc(size_t num_bytes)
{
	std::shared_ptr<char> ret;
	auto it = bufs.find(num_bytes);
	std::deque<char *> *q;
	if (it == bufs.end()) {
		bufs.insert(std::pair<size_t, std::deque<char *> >(num_bytes,
					std::deque<char *>()));
		it = bufs.find(num_bytes);
		assert(it != bufs.end());
		q = &it->second;
	}
	else
		q = &it->second;
	if (q->empty()) {
		num_allocs++;
		// When the piece of memory is deallocated, it'll be pushed back to
		// the deque.
		ret = std::shared_ptr<char>((char *) memalign(PAGE_SIZE, num_bytes),
				local_deleter(*q));
	}
	else {
		char *tmp = q->front();
		q->pop_front();
		ret = std::shared_ptr<char>(tmp, local_deleter(*q));
	}
	return ret;
}

void local_mem_buffer::clear_local_bufs()
{
	// This method may be called in another thread.
	// TODO maybe I should use a lock to protect the per-thread queue.
	for (auto it = bufs.begin(); it != bufs.end(); it++) {
		std::deque<char *> &q = it->second;
		num_frees += q.size();
		while (!q.empty()) {
			char *buf = q.front();
			q.pop_front();
			free(buf);
		}
	}
	portions.clear();
	// TODO we should also clear all entries in `buf'. But it causes
	// memory deallocation error. Why?
}

std::shared_ptr<char> local_mem_buffer::alloc(size_t num_bytes)
{
	if (!initialized)
		return NULL;

	void *addr = pthread_getspecific(mem_key);
	if (addr)
		return ((local_mem_buffer *) addr)->_alloc(num_bytes);
	else {
		local_mem_buffer *buf = new local_mem_buffer();
		pthread_setspecific(mem_key, buf);
		mem_lock.lock();
		mem_set.push_back(buf);
		mem_lock.unlock();
		return buf->_alloc(num_bytes);
	}
}

bool local_mem_buffer::init()
{
	int ret = pthread_key_create(&mem_key, NULL);
	if (ret == 0) {
		assert(mem_set.empty());
		initialized = true;
		return true;
	}
	else {
		BOOST_LOG_TRIVIAL(error)
			<< "can't create the pthread key for local memory buffer";
		return false;
	}
}

void local_mem_buffer::destroy()
{
	mem_lock.lock();
	// We should destroy all local memory buffers because we don't know
	// if we need memory buffers of the same size next time.
	// Destroy all local mem buffers.
	for (auto it = mem_set.begin(); it != mem_set.end(); it++)
		delete *it;
	mem_set.clear();
	mem_lock.unlock();
	initialized = false;
	pthread_key_delete(mem_key);
}

void local_mem_buffer::clear_bufs()
{
	mem_lock.lock();
	for (auto it = mem_set.begin(); it != mem_set.end(); it++)
		(*it)->clear_local_bufs();
	mem_lock.unlock();
}

void local_mem_buffer::_cache_portion(long key,
		local_matrix_store::const_ptr portion)
{
	auto it = portions.find(key);
	if (it == portions.end())
		portions.insert(std::pair<long, local_matrix_store::const_ptr>(
					key, portion));
	// We only buffer the most recent portion for a matrix.
	// This is enough for the current mapply operations.
	else
		it->second = portion;
}

local_matrix_store::const_ptr local_mem_buffer::_get_mat_portion(long key)
{
	auto it = portions.find(key);
	if (it == portions.end())
		return NULL;
	else
		return it->second;
}

void local_mem_buffer::cache_portion(long key,
		local_matrix_store::const_ptr portion)
{
	if (!initialized)
		return;

	void *addr = pthread_getspecific(mem_key);
	if (addr)
		((local_mem_buffer *) addr)->_cache_portion(key, portion);
	else {
		local_mem_buffer *buf = new local_mem_buffer();
		pthread_setspecific(mem_key, buf);
		mem_lock.lock();
		mem_set.push_back(buf);
		mem_lock.unlock();
		buf->_cache_portion(key, portion);
	}
}

local_matrix_store::const_ptr local_mem_buffer::get_mat_portion(long key)
{
	if (!initialized)
		return NULL;

	void *addr = pthread_getspecific(mem_key);
	if (addr)
		return ((local_mem_buffer *) addr)->_get_mat_portion(key);
	else
		return NULL;
}

spin_lock local_mem_buffer::mem_lock;
std::vector<local_mem_buffer *> local_mem_buffer::mem_set;
std::atomic<bool> local_mem_buffer::initialized;
pthread_key_t local_mem_buffer::mem_key;

}

}
