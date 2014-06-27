#ifndef __SIMPLE_KV_STORE_H__
#define __SIMPLE_KV_STORE_H__

/**
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

#include "io_interface.h"
#include "container.h"
#include "cache.h"
#include "slab_allocator.h"

template<class ValueType, class ValueTaskType>
class KV_compute: public user_compute
{
	embedded_array<ValueTaskType> tasks;
	int num_tasks;
	bool has_run;
public:
	KV_compute(compute_allocator *alloc): user_compute(alloc) {
		num_tasks = 0;
		has_run = false;
	}

	void init(ValueTaskType tasks[], int num) {
		this->tasks.resize(num);
		for (int i = 0; i < num; i++) {
			this->tasks[i] = tasks[i];
		}
		this->num_tasks = num;
		has_run = false;
	}

	virtual int serialize(char *buf, int size) const {
		return 0;
	}

	virtual int get_serialized_size() const {
		return 0;
	}

	virtual void run(page_byte_array &arr) {
		off_t start_off = arr.get_offset() / sizeof(ValueType);
		off_t end_off = (arr.get_offset() + arr.get_size()) / sizeof(ValueType);
		for (int i = 0; i < num_tasks; i++) {
			off_t idx = tasks[i].get_idx();
			assert(idx >= start_off && idx < end_off);
			ValueType v = arr.get<ValueType>(idx - start_off);
			tasks[i].run(v);
		}
		has_run = true;
	}

	virtual bool has_completed() const {
		return has_run;
	}

	virtual int has_requests() const {
		return false;
	}

	virtual request_range get_next_request() {
		assert(0);
	}
};

template<class ValueType, class ValueTaskType>
class KV_compute_allocator: public compute_allocator
{
	class compute_initializer: public obj_initiator<KV_compute<ValueType, ValueTaskType> >
	{
		KV_compute_allocator<ValueType, ValueTaskType> *alloc;
	public:
		compute_initializer(
				KV_compute_allocator<ValueType, ValueTaskType> *alloc) {
			this->alloc = alloc;
		}

		virtual void init(KV_compute<ValueType, ValueTaskType> *obj) {
			new (obj) KV_compute<ValueType, ValueTaskType>(alloc);
		}
	};

	class compute_destructor: public obj_destructor<KV_compute<ValueType, ValueTaskType> >
	{
	public:
		void destroy(KV_compute<ValueType, ValueTaskType> *obj) {
			obj->~KV_compute<ValueType, ValueTaskType>();
		}
	};

	obj_allocator<KV_compute<ValueType, ValueTaskType> > allocator;
public:
	KV_compute_allocator(int node_id): allocator("KV_compute_allocator",
			node_id, 1024 * 1024, params.get_max_obj_alloc_size(),
			typename obj_initiator<KV_compute<ValueType, ValueTaskType> >::ptr(
				new compute_initializer(this)),
			typename obj_destructor<KV_compute<ValueType, ValueTaskType> >::ptr(
				new compute_destructor())) {
	}

	virtual user_compute *alloc() {
		return allocator.alloc_obj();
	}

	virtual void free(user_compute *obj) {
		allocator.free((KV_compute<ValueType, ValueTaskType> *) obj);
	}
};

/**
 * This is a simple key-value store over a single file.
 * It supports only one type of user-defined tasks on values and can be used
 * in one thread.
 * User tasks are executed asynchronously.
 */
template<class ValueType, class TaskType>
class simple_KV_store
{
	io_interface::ptr io;
	KV_compute_allocator<ValueType, TaskType> alloc;

	embedded_array<io_request> req_buf;
	int num_reqs;

	void add_io_request(io_request &req) {
		if (req_buf.get_capacity() <= num_reqs)
			req_buf.resize(req_buf.get_capacity() * 2);
		req_buf[num_reqs] = req;
		num_reqs++;
	}

	void flush_io_requests() {
		printf("issue %d requests\n", num_reqs);
		io->access(req_buf.data(), num_reqs);
		num_reqs = 0;
	}

	simple_KV_store(io_interface::ptr io): alloc(io->get_node_id()) {
		this->io = io;
		num_reqs = 0;
		assert(PAGE_SIZE % sizeof(ValueType) == 0);
	}
public:
	typedef std::shared_ptr<simple_KV_store<ValueType, TaskType> > ptr;

	static ptr create(io_interface::ptr io) {
		return ptr(new simple_KV_store<ValueType, TaskType>(io));
	}

	/**
	 * Serve user requests asynchronously.
	 * It assumes that all user tasks have been sorted based on the locations
	 * of the entries accessed by the user tasks.
	 */
	void async_request(TaskType tasks[], int num) {
		if (num == 0)
			return;

		// Each time we issue a single request to serve as many user tasks
		// as possible. Each time we issue a request to read at least one
		// page. We'll merge requests if the pages touched by user tasks
		// in the input array are adjacent to each other.

		// The offset of the first page accessed by the I/O request.
		off_t first_page_off = ROUND_PAGE(tasks[0].get_idx() * sizeof(ValueType));
		// The offset of the last page accessed by the I/O request.
		// The page is excluded by the I/O request.
		off_t last_page_off = first_page_off + PAGE_SIZE;
		// The first task covered by the I/O request.
		int first_task_idx = 0;
		int i = 1;
		for (; i < num; i++) {
			off_t task_off = tasks[i].get_idx() * sizeof(ValueType);
			off_t page_off = ROUND_PAGE(task_off);
			assert(page_off >= first_page_off);
			// If the task access the page covered by the I/O request.
			if (page_off < last_page_off)
				continue;
			// If the task access the page right behind the range covered
			// by the I/O request.
			if (page_off == last_page_off) {
				last_page_off += PAGE_SIZE;
				continue;
			}

			// The user task accesses a page far away from the range covered
			// by the current I/O request.
			// Issue a request to serve tasks between `first_task_idx` and `i`.
			KV_compute<ValueType, TaskType> *compute
				= (KV_compute<ValueType, TaskType> *) alloc.alloc();
			compute->init(&tasks[first_task_idx], i - first_task_idx);
			data_loc_t loc(io->get_file_id(), first_page_off);
			io_request req(compute, loc, last_page_off - first_page_off, READ);
			add_io_request(req);

			// Re-initialize the range covered by the new I/O request.
			first_page_off = page_off;
			last_page_off = first_page_off + PAGE_SIZE;
			first_task_idx = i;
		}

		assert(first_task_idx < num);
		KV_compute<ValueType, TaskType> *compute
			= (KV_compute<ValueType, TaskType> *) alloc.alloc();
		compute->init(&tasks[first_task_idx], num - first_task_idx);
		data_loc_t loc(io->get_file_id(), first_page_off);
		io_request req(compute, loc, last_page_off - first_page_off, READ);
		add_io_request(req);
		flush_io_requests();
	}
};

#endif
