#ifndef __EM_VECTOR_H__
#define __EM_VECTOR_H__

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

#include <unordered_map>
#include <memory>
#include <atomic>

#include "io_interface.h"

#include "vec_store.h"
#include "local_vec_store.h"
#include "mem_worker_thread.h"
#include "EM_object.h"

namespace fm
{

namespace detail
{

class matrix_store;
class EM_vec_store;

/*
 * This sorts the first external-memory vector in `vecs' and shuffles
 * the remaining vectors accordingly.
 */
std::vector<std::shared_ptr<EM_vec_store> > sort(
		const std::vector<std::shared_ptr<const EM_vec_store> > &vecs);

class EM_vec_store: public vec_store, public EM_object
{
	struct empty_free {
		void operator()(EM_vec_store *) {
		}
	};

	class file_holder {
		std::string file_name;
		file_holder(const std::string &name) {
			this->file_name = name;
		}
	public:
		typedef std::shared_ptr<file_holder> ptr;
		static ptr create_temp(size_t num_bytes);

		~file_holder();
		std::string get_name() const {
			return file_name;
		}
	};

	class io_set {
		safs::file_io_factory::shared_ptr factory;
		// This keeps an I/O instance for each thread.
		std::unordered_map<thread *, safs::io_interface::ptr> thread_ios;
		pthread_key_t io_key;
		pthread_spinlock_t io_lock;
	public:
		typedef std::shared_ptr<io_set> ptr;
		io_set(safs::file_io_factory::shared_ptr factory);
		~io_set();

		safs::io_interface::ptr create_io();
		// This returns the I/O instance for the curr thread.
		safs::io_interface &get_curr_io() const;
		// Test if the current thread has an I/O instance for the vector.
		bool has_io() const;
	};

	file_holder::ptr holder;
	io_set::ptr ios;

	EM_vec_store(size_t length, const scalar_type &type);
	EM_vec_store(const EM_vec_store &store);
public:
	typedef std::shared_ptr<EM_vec_store> ptr;
	typedef std::shared_ptr<const EM_vec_store> const_ptr;

	static ptr cast(vec_store::ptr vec);
	static const_ptr cast(vec_store::const_ptr vec);

	static ptr create(size_t length, const scalar_type &type) {
		return ptr(new EM_vec_store(length, type));
	}

	~EM_vec_store();

	size_t get_byte_off(size_t entry_off) const {
		return entry_off * get_entry_size();
	}

	virtual bool resize(size_t length);

	virtual bool append(std::vector<vec_store::const_ptr>::const_iterator vec_it,
			std::vector<vec_store::const_ptr>::const_iterator vec_end);
	virtual bool append(const vec_store &vec);
	virtual vec_store::ptr deep_copy() const;
	virtual vec_store::ptr shallow_copy();
	virtual vec_store::const_ptr shallow_copy() const;

	virtual size_t get_portion_size() const;
	virtual local_vec_store::const_ptr get_portion(off_t loc, size_t size) const;
	virtual local_vec_store::ptr get_portion(off_t loc, size_t size);
	/*
	 * This is different from the one used in the memory data container.
	 * This interface accepts a portion compute object, which is invoked
	 * when the data of the portion is ready in memory. The returned
	 * local vector store doesn't have invalid data right after it's
	 * returned.
	 */
	virtual local_vec_store::ptr get_portion_async(off_t start,
			size_t size, portion_compute::ptr compute) const;
	virtual std::vector<local_vec_store::ptr> get_portion_async(
			const std::vector<std::pair<off_t, size_t> > &locs,
			portion_compute::ptr compute) const;
	/*
	 * Write the data in the local buffer to some portion in the vector.
	 * The location is indicated in the local buffer. However, a user
	 * can redirect the location by providing `off' as a parameter.
	 */
	virtual void write_portion(local_vec_store::const_ptr portion, off_t off = -1);

	virtual void reset_data();
	virtual void set_data(const set_vec_operate &op);

	virtual vec_store::ptr sort_with_index();
	virtual void sort();
	virtual bool is_sorted() const;

	virtual std::shared_ptr<const matrix_store> conv2mat(size_t nrow,
			size_t ncol, bool byrow) const;

	virtual safs::io_interface::ptr create_io();

	friend std::vector<EM_vec_store::ptr> sort(
			const std::vector<EM_vec_store::const_ptr> &vecs);
};

namespace EM_sort_detail
{

/*
 * This priority queue helps to sort data in the ascending order.
 */
class anchor_prio_queue
{
	struct anchor_struct
	{
		local_buf_vec_store::const_ptr local_anchors;
		int id;
		off_t curr_off;
	};

	class anchor_ptr_less {
		const bulk_operate *gt;
	public:
		anchor_ptr_less(const scalar_type &type) {
			gt = type.get_basic_ops().get_op(basic_ops::op_idx::GT);;
		}

		bool operator()(const anchor_struct *anchor1,
				const anchor_struct *anchor2) const {
			bool ret;
			gt->runAA(1, anchor1->local_anchors->get(anchor1->curr_off),
					anchor2->local_anchors->get(anchor2->curr_off), &ret);
			return ret;
		}
	};

	const size_t sort_buf_size;
	const size_t anchor_gap_size;
	std::vector<anchor_struct> anchor_bufs;
	typedef std::priority_queue<anchor_struct *, std::vector<anchor_struct *>,
			anchor_ptr_less> anchor_queue_t;
	anchor_queue_t queue;

	off_t get_anchor_off(const anchor_struct &anchor) const;
public:
	typedef std::shared_ptr<anchor_prio_queue> ptr;

	anchor_prio_queue(const std::vector<local_buf_vec_store::ptr> &anchor_vals,
			size_t sort_buf_size, size_t anchor_gap_size);
	scalar_variable::ptr get_min_frontier() const;
	size_t get_anchor_gap_size() const {
		return anchor_gap_size;
	}

	/*
	 * Here we pop a set of chunks of data whose values are the potentially
	 * the smallest.
	 */
	std::vector<off_t> pop(size_t size);
	/*
	 * Fetch the first anchors from all queues.
	 */
	std::vector<off_t> fetch_all_first();
};

/*
 * This class keeps some summaries of sorted vectors in memory, so later on
 * we can use these summaries to pinpoint the right data for merging.
 * To generate a summary of a vector, we put some conceptional anchors in
 * the vector and keep the values in the anchor location.
 */
class sort_portion_summary
{
	const size_t sort_buf_size;
	const size_t anchor_gap_size;
	std::vector<local_buf_vec_store::ptr> anchor_vals;
public:
	sort_portion_summary(size_t num_sort_bufs, size_t sort_buf_size,
			size_t anchor_gap_size);
	void add_portion(local_buf_vec_store::const_ptr sorted_buf);
	anchor_prio_queue::ptr get_prio_queue() const;

	local_vec_store::const_ptr get_anchor_vals(off_t idx) const {
		return anchor_vals[idx];
	}
	size_t get_num_bufs() const {
		return anchor_vals.size();
	}
};

/*
 * This class sorts a portion of data in the vector and writes the sorted
 * results to disks.
 */
class EM_vec_sort_compute: public portion_compute
{
	// The portions are from different vectors.
	// Each portion has to be read from the disks.
	std::vector<local_buf_vec_store::ptr> portions;
	// Where the sorted portion written to.
	std::vector<EM_vec_store::ptr> to_vecs;
	sort_portion_summary &summary;
	// The number of portions that have been read from the disks.
	size_t num_completed;
public:
	EM_vec_sort_compute(const std::vector<EM_vec_store::ptr> &vecs,
			sort_portion_summary &_summary): summary(_summary) {
		this->to_vecs = vecs;
		num_completed = 0;
	}
	virtual void run(char *buf, size_t size);
	void set_bufs(std::vector<local_buf_vec_store::ptr> portions) {
		this->portions = portions;
	}
};

class EM_vec_merge_dispatcher;

/*
 * This class merges data read from disks and writes it back to disks.
 */
class EM_vec_merge_compute: public portion_compute
{
	// This defines the container with all the portions used for merging
	// a vector. The first buffer in the set may be the leftover from
	// the previous merge.
	typedef std::vector<local_buf_vec_store::const_ptr> merge_set_t;
	std::vector<merge_set_t> stores;
	EM_vec_merge_dispatcher &dispatcher;
	size_t num_completed;
	// It's not the same as the number of local buffers in `stores' because
	// `stores' may contain the buffer with the leftdata from the previous
	// merge.
	size_t num_expected;
public:
	EM_vec_merge_compute(
			const std::vector<local_buf_vec_store::ptr> &prev_leftovers,
			EM_vec_merge_dispatcher &_dispatcher);
	virtual void run(char *buf, size_t size);
	void set_bufs(const std::vector<merge_set_t> &bufs);
};

/*
 * The two functions compute the sort buffer size and anchor gap size.
 */
std::pair<size_t, size_t> cal_sort_buf_size(const scalar_type &type);
std::pair<size_t, size_t> cal_sort_buf_size(
		const std::vector<const scalar_type *> &types);

}

}

}

#endif
