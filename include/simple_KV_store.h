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

template<class ValueType, class TaskType>
class simple_KV_store;

template<class ValueType, class ValueTaskType>
class KV_compute: public user_compute
{
	embedded_array<ValueTaskType> tasks;
	int num_tasks;
	int num_added_tasks;
	bool has_run;
	simple_KV_store<ValueType, ValueTaskType> *store;
public:
	KV_compute(simple_KV_store<ValueType, ValueTaskType> *store,
			compute_allocator *alloc): user_compute(alloc) {
		this->store = store;
		num_tasks = 0;
		num_added_tasks = 0;
		has_run = false;
	}

	bool has_tasks() const {
		return num_tasks > 0;
	}

	void add_task(const ValueTaskType &task) {
		num_added_tasks++;
		if (tasks.get_capacity() <= num_tasks)
			tasks.resize(num_tasks * 2);
		if (num_tasks == 0 || !tasks[num_tasks - 1].merge(task))
			tasks[num_tasks++] = task;
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
			int num_entries = tasks[i].get_num_entries();
			off_t it_start = (idx - start_off) * sizeof(ValueType);
			off_t it_end = it_start + sizeof(ValueType) * num_entries;
			page_byte_array::seq_const_iterator<ValueType> it
				= arr.get_seq_iterator<ValueType>(it_start, it_end);
			tasks[i].run(it);
		}
		has_run = true;
		store->complete_tasks(num_added_tasks);
	}

	virtual bool has_completed() {
		return has_run;
	}

	virtual int has_requests() {
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
		simple_KV_store<ValueType, ValueTaskType> *store;
	public:
		compute_initializer(
				simple_KV_store<ValueType, ValueTaskType> *store,
				KV_compute_allocator<ValueType, ValueTaskType> *alloc) {
			this->store = store;
			this->alloc = alloc;
		}

		virtual void init(KV_compute<ValueType, ValueTaskType> *obj) {
			new (obj) KV_compute<ValueType, ValueTaskType>(store, alloc);
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
	KV_compute_allocator(simple_KV_store<ValueType, ValueTaskType> *store,
			int node_id): allocator("KV_compute_allocator",
			node_id, false, 1024 * 1024, params.get_max_obj_alloc_size(),
			typename obj_initiator<KV_compute<ValueType, ValueTaskType> >::ptr(
				new compute_initializer(store, this)),
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

	struct task_comp_larger {
		bool operator()(const TaskType task1, const TaskType task2) {
			return task1.get_idx() >= task2.get_idx();
		}
	};
	struct task_comp_smaller {
		bool operator()(const TaskType &task1, const TaskType &task2) {
			return task1.get_idx() <= task2.get_idx();
		}
	};
	struct task_less {
		bool operator()(const TaskType &task1, const TaskType &task2) {
			return task1.get_idx() < task2.get_idx();
		}
	};

	std::deque<TaskType> task_buf;
	ssize_t num_pending_tasks;

	embedded_array<io_request> req_buf;
	int num_reqs;

	void add_io_request(io_request &req) {
		if (req_buf.get_capacity() <= num_reqs)
			req_buf.resize(req_buf.get_capacity() * 2);
		req_buf[num_reqs] = req;
		num_reqs++;
	}

	void flush_io_requests() {
		io->access(req_buf.data(), num_reqs);
		num_reqs = 0;
	}

	void sort_tasks() {
		if (task_buf.size() < 2)
			return;
		if (task_buf[0].get_idx() > task_buf[1].get_idx()) {
			if (std::is_sorted(task_buf.begin(), task_buf.end(), task_comp_larger()))
				return;
		}
		else
			if (std::is_sorted(task_buf.begin(), task_buf.end(), task_comp_smaller()))
				return;
		std::sort(task_buf.begin(), task_buf.end(), task_less());
	}

	simple_KV_store(io_interface::ptr io): alloc(this, io->get_node_id()) {
		this->io = io;
		num_reqs = 0;
		assert(PAGE_SIZE % sizeof(ValueType) == 0);
		num_pending_tasks = 0;
	}
public:
	typedef std::shared_ptr<simple_KV_store<ValueType, TaskType> > ptr;

	static ptr create(io_interface::ptr io) {
		return ptr(new simple_KV_store<ValueType, TaskType>(io));
	}

	void flush_requests() {
		if (task_buf.empty())
			return;
		// We can't issue any I/O requests to SAFS.
		if (io->num_pending_ios() >= io->get_max_num_pending_ios())
			return;

		// Each time we issue a single request to serve as many user tasks
		// as possible. Each time we issue a request to read at least one
		// page. We'll merge requests if the pages touched by user tasks
		// in the input array are adjacent to each other.

		sort_tasks();
		int avail_io_slots = io->get_remaining_io_slots();
		// The offset of the first page accessed by the I/O request.
		const TaskType task = task_buf.front();
		task_buf.pop_front();
		off_t first_page_off = ROUND_PAGE(task.get_idx() * sizeof(ValueType));
		// The offset of the last page accessed by the I/O request.
		// The page is excluded by the I/O request.
		off_t last_page_off = ROUNDUP_PAGE((task.get_idx()
					+ task.get_num_entries()) * sizeof(ValueType));
		KV_compute<ValueType, TaskType> *compute
			= (KV_compute<ValueType, TaskType> *) alloc.alloc();
		compute->add_task(task);
		while (!task_buf.empty()
				// We have one more request outside of the while loop.
				&& num_reqs < avail_io_slots - 1) {
			const TaskType task = task_buf.front();
			task_buf.pop_front();
			off_t page_off = ROUND_PAGE(task.get_idx() * sizeof(ValueType));
			off_t end_page_off = ROUNDUP_PAGE((task.get_idx()
					+ task.get_num_entries()) * sizeof(ValueType));
			// If the range overlaps.
			if ((page_off >= first_page_off && page_off <= last_page_off)
					|| (end_page_off >= first_page_off && end_page_off <= last_page_off)) {
				compute->add_task(task);
				first_page_off = std::min(page_off, first_page_off);
				last_page_off = std::max(end_page_off, last_page_off);
				continue;
			}

			// The user task accesses a page far away from the range covered
			// by the current I/O request.
			data_loc_t loc(io->get_file_id(), first_page_off);
			io_request req(compute, loc, last_page_off - first_page_off, READ);
			add_io_request(req);

			// Re-initialize the range covered by the new I/O request.
			compute = (KV_compute<ValueType, TaskType> *) alloc.alloc();
			compute->add_task(task);
			first_page_off = ROUND_PAGE(task.get_idx() * sizeof(ValueType));
			last_page_off = end_page_off;
		}

		assert(compute->has_tasks());
		data_loc_t loc(io->get_file_id(), first_page_off);
		io_request req(compute, loc, last_page_off - first_page_off, READ);
		add_io_request(req);
		flush_io_requests();
	}

	/**
	 * Serve user requests asynchronously.
	 */
	void async_request(TaskType &task) {
		if (task_buf.empty() || !task_buf.back().merge(task)) {
			task_buf.push_back(task);
			num_pending_tasks++;
		}
	}

	void complete_tasks(int num_tasks) {
		num_pending_tasks -= num_tasks;
		assert(num_pending_tasks >= 0);
	}

	size_t get_num_pending_tasks() const {
		return num_pending_tasks;
	}
};

#endif
