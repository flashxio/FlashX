#ifndef __COMP_IO_SCHEDULER_H__
#define __COMP_IO_SCHEDULER_H__

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

/**
 * This class defines an interface of a I/O scheduler for requests generated
 * by user tasks. It is only used in the page cache.
 * Right now, an I/O scheduler is associated with an I/O instance and can
 * only schedule I/O requests issued by the I/O instance.
 */
class comp_io_scheduler
{
	// If user computation has generated too many requests, we may not
	// complete a user computation (we can't fetch all requests generated
	// by it). We have to keep the incomplete computation here, and we
	// will try to fetch more requests from it later.
	fifo_queue<user_compute *> incomplete_computes;
	// The I/O instance where the I/O scheduler works on.
	io_interface *io;

public:
	/**
	 * This class iterates the user tasks managed by the I/O scheduler.
	 */
	class compute_iterator {
		fifo_queue<user_compute *>::const_iterator it;
	public:
		compute_iterator(const fifo_queue<user_compute *> &computes,
				bool end): it(&computes) {
			if (end)
				it = computes.get_end();
		}

		/**
		 * This method gets the current user task.
		 * \return the current user task.
		 */
		user_compute *operator*() const {
			return *it;
		}

		/**
		 * This method moves the iterator forward by one.
		 * \return the reference to this iterator.
		 */
		compute_iterator &operator++() {
			++it;
			return *this;
		}

		/**
		 * This method tests whether the two iterators are the same.
		 * \param it the other iterator.
		 * \return true if they are the same.
		 */
		bool operator==(const compute_iterator &it) const {
			return this->it == it.it;
		}

		/**
		 * This method tests whether the two iterator aren't the same.
		 * \param it the other iterator.
		 * \return true if they aren't the same.
		 */
		bool operator!=(const compute_iterator &it) const {
			return this->it != it.it;
		}
	};

	/**
	 * TODO
	 * The iterator should only iterate on the user tasks with requests.
	 */

	/**
	 * This method gets the iterator pointing to the first user task
	 * managed by the I/O scheduler.
	 * \return the iterator to the beginning.
	 */
	compute_iterator get_begin() const {
		return compute_iterator(incomplete_computes, false);
	}

	/**
	 * This method gets the iterator pointing to the last user task
	 * managed by the I/O scheduler.
	 * \return the iterator to the end.
	 */
	compute_iterator get_end() const {
		return compute_iterator(incomplete_computes, true);
	}

	/**
	 * The constructor of an I/O scheduler.
	 * \param node_id the NUMA node where the I/O scheduler runs.
	 */
	comp_io_scheduler(int node_id);

	virtual ~comp_io_scheduler() {
		assert(incomplete_computes.is_empty());
	}

	/**
	 * This method gets multiple I/O requests.
	 * \param reqs the buffer where the fetched requests are stored.
	 * \return the number of fetched requests.
	 */
	virtual size_t get_requests(fifo_queue<io_request> &reqs) = 0;

	/**
	 * This method sets the I/O instance that the I/O scheduler is
	 * associated with.
	 * \param io the I/O instance.
	 */
	void set_io(io_interface *io) {
		this->io = io;
	}

	/**
	 * This method gets the I/O instance that the I/O scheduler is
	 * associated with.
	 * \return the I/O instance.
	 */
	io_interface *get_io() const {
		return io;
	}

	/**
	 * This method adds a user task to the I/O scheduler so that the I/O
	 * scheduler can schedule it.
	 * \param compute the user task.
	 */
	void add_compute(user_compute *compute) {
		// We have to make sure the computation has requested new data
		// successfully, otherwise, it may not be executed again.
		if (!compute->test_flag(user_compute::IN_QUEUE)) {
			compute->inc_ref();
			incomplete_computes.push_back(compute);
			compute->set_flag(user_compute::IN_QUEUE, true);
		}
	}

	/**
	 * This method is a static method. It destroys a user task.
	 * \param compute the user task.
	 */
	static void delete_compute(user_compute *compute) {
		assert(compute->test_flag(user_compute::IN_QUEUE));
		compute->set_flag(user_compute::IN_QUEUE, false);
		compute->dec_ref();
		// the user compute isn't free'd here, it's the user's responsiblity
		// of freeing it.
		if (compute->get_ref() == 0) {
			compute_allocator *alloc = compute->get_allocator();
			alloc->free(compute);
		}
	}

	/**
	 * This method gets the number of user tasks that haven't been completed.
	 * \return the number of user tasks that haven't been completed.
	 */
	size_t get_num_incomplete_computes() {
		return incomplete_computes.get_num_entries();
	}

	/**
	 * This method tests whether there are user tasks that still haven't
	 * been completed.
	 * \return indicates whether there are user tasks that still haven't
	 * been completed.
	 */
	bool is_empty() {
		return incomplete_computes.is_empty();
	}

	/**
	 * This method garbage collect user compute tasks that have been
	 * completed.
	 */
	void gc_computes();
};

#endif
