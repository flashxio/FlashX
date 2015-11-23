#ifndef __SLAB_ALLOCATOR_H__
#define __SLAB_ALLOCATOR_H__

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

#include <stdlib.h>
#include <assert.h>

#include <memory>

#include "concurrency.h"
#include "aligned_allocator.h"
#include "container.h"

static const int SLAB_LOCAL_BUF_SIZE = 100;

class slab_allocator
{
public:
	class linked_obj {
		linked_obj *next;
	public:
		linked_obj() {
			next = NULL;
		}

		void add(linked_obj *obj) {
			obj->next = this->next;
			this->next = obj;
		}

		linked_obj *get_next() {
			return next;
		}

		void set_next(linked_obj *obj) {
			next = obj;
		}
	};

	class linked_obj_list {
		linked_obj head;
		linked_obj *end;		// point to the last element or NULL.
		int size;
	public:
		linked_obj_list() {
			size = 0;
			end = NULL;
		}

		void add(linked_obj *obj) {
			if (get_size() == 0) {
				end = obj;
				assert(head.get_next() == NULL);
			}
			head.add(obj);
			size++;
		}

		void add_list(linked_obj_list *list) {
			if (list->get_size() == 0)
				return;
			if (get_size() == 0) {
				size = list->get_size();
				head = list->head;
				end = list->end;
				list->head.set_next(NULL);
				list->end = NULL;
				list->size = 0;
			}
			else {
				list->end->set_next(head.get_next());
				head = list->head;
				size += list->size;
				list->head.set_next(NULL);
				list->end = NULL;
				list->size = 0;
			}
		}

		/*
		 * It removes the specified number of objects from the list
		 * and returns them in a form of linked list..
		 */
		linked_obj *pop(int num) {
			linked_obj *obj = head.get_next();
			int num_pop = 0;
			linked_obj *prev_obj = NULL;
			for (int i = 0; i < num && obj != NULL; i++) {
				prev_obj = obj;
				obj = obj->get_next();
				num_pop++;
			}
			linked_obj *ret = head.get_next();
			if (prev_obj)
				prev_obj->set_next(NULL);

			head.set_next(obj);
			size -= num_pop;
			assert(size >= 0);
			if (size == 0)
				end = NULL;
			return ret;
		}

		int get_size() {
			return size;
		}
	};

private:
	const int obj_size;
	// the size to increase each time there aren't enough objects
	const long increase_size;
	// the maximal size of all objects
	const long max_size;
	const int node_id;
	const int local_buf_size;
	const bool thread_safe;

	linked_obj_list list;
	// the current size of memory used by the allocator.
	atomic_number<long> curr_size;
	bool init;
	bool pinned;

	std::vector<char *> alloc_bufs;

	pthread_spinlock_t lock;
	// The buffers pre-allocated to serve allocation requests
	// from the local threads.
	pthread_key_t local_buf_key;

	thread_safe_FIFO_queue<fifo_queue<char *> *> per_thread_queues;

	std::string name;
	static atomic_integer alloc_counter;

	fifo_queue<char *> *get_local_buf() {
		fifo_queue<char *> *local_buf_refs
			= (fifo_queue<char *> *) pthread_getspecific(local_buf_key);
		if (local_buf_refs == NULL) {
			local_buf_refs = fifo_queue<char *>::create(node_id, local_buf_size);
			per_thread_queues.add(&local_buf_refs, 1);
			pthread_setspecific(local_buf_key, local_buf_refs);
		}
		return local_buf_refs;
	}
#ifdef MEMCHECK
	aligned_allocator allocator;
#endif
public:
	slab_allocator(const std::string &name, int _obj_size, long _increase_size,
			// We allow pages to be pinned when allocated.
			long _max_size, int _node_id, bool init = false, bool pinned = false,
			int _local_buf_size = SLAB_LOCAL_BUF_SIZE, bool thread_safe = true);

	virtual ~slab_allocator();

	int get_obj_size() const {
		return obj_size;
	}

	int alloc(char **objs, int num);

	void free(char **objs, int num);

	char *alloc();

	void free(char *obj);

	// Use it carefully. It's not thread-safe.
	bool contains(const char *addr) const {
		for (unsigned i = 0; i < alloc_bufs.size(); i++) {
			if (addr >= alloc_bufs[i] && addr < alloc_bufs[i]
					+ increase_size)
				return true;
		}
		return false;
	}

	long get_max_size() const {
		return max_size;
	}

	long get_curr_size() const {
		return curr_size.get();
	}

	const std::string &get_name() const {
		return name;
	}
};

template<class T>
class obj_initiator
{
public:
	typedef typename std::unique_ptr<obj_initiator<T> > ptr;
	virtual void init(T *obj) = 0;
};

template<class T>
class obj_destructor
{
public:
	typedef typename std::unique_ptr<obj_destructor<T> > ptr;
	virtual void destroy(T *obj) = 0;
};

template<class T>
class default_obj_initiator: public obj_initiator<T>
{
public:
	void init(T *obj) {
		T tmp;
		*obj = tmp;
	}
};

template<class T>
class default_obj_destructor: public obj_destructor<T>
{
public:
	void destroy(T *obj) {
	}
};

template<class T>
class obj_allocator: public slab_allocator
{
	typename obj_initiator<T>::ptr initiator;
	typename obj_destructor<T>::ptr destructor;
public:
	obj_allocator(const std::string &name, int node_id, bool thread_safe,
			long increase_size, long max_size = INT_MAX,
			typename obj_initiator<T>::ptr initiator = typename obj_initiator<T>::ptr(
				new default_obj_initiator<T>()),
			typename obj_destructor<T>::ptr destructor = typename obj_destructor<T>::ptr(
				new default_obj_destructor<T>())
				// leave some space for linked_obj, so the values in an object
				// won't be modified.
				): slab_allocator(name, sizeof(T) + sizeof(slab_allocator::linked_obj),
				increase_size, max_size, node_id, true, false, SLAB_LOCAL_BUF_SIZE,
				thread_safe) {
		assert(increase_size <= max_size);
		this->initiator = std::move(initiator);
		this->destructor = std::move(destructor);
	}

	virtual int alloc_objs(T **objs, int num) {
		char *addrs[num];
		int ret = slab_allocator::alloc(addrs, num);
		for (int i = 0; i < ret; i++) {
			objs[i] = (T *) (addrs + sizeof(slab_allocator::linked_obj));
			initiator->init(objs[i]);
		}
		return ret;
	}

	virtual T *alloc_obj() {
		char *addr = slab_allocator::alloc();
		if (addr == NULL) {
			return NULL;
		}

		T *obj = (T *) (addr + sizeof(slab_allocator::linked_obj));
		initiator->init(obj);
		return obj;
	}

	virtual void free(T **objs, int num) {
		char *addrs[num];
		for (int i = 0; i < num; i++) {
			destructor->destroy(objs[i]);
			addrs[i] = ((char *) objs[i]) - sizeof(slab_allocator::linked_obj);
		}
		slab_allocator::free(addrs, num);
	}

	virtual void free(T *obj) {
		char *addr = ((char *) obj) - sizeof(slab_allocator::linked_obj);
		destructor->destroy(obj);
		slab_allocator::free(addr);
	}
};

#endif
