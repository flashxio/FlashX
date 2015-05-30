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

#include <memory>
#include <atomic>

#include "io_interface.h"

#include "vec_store.h"
#include "local_vec_store.h"
#include "mem_worker_thread.h"

namespace fm
{

namespace detail
{

class matrix_store;

class EM_vec_store: public vec_store
{
	safs::file_io_factory::shared_ptr factory;

	EM_vec_store(size_t length, const scalar_type &type);
public:
	typedef std::shared_ptr<EM_vec_store> ptr;

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

	virtual void reset_data();
	virtual void set_data(const set_vec_operate &op);

	virtual vec_store::ptr sort_with_index();
	virtual void sort();
	virtual bool is_sorted() const;

	virtual std::shared_ptr<const matrix_store> conv2mat(size_t nrow,
			size_t ncol, bool byrow) const;
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
		const bulk_operate *gt;

		bool operator<(const anchor_struct &anchor) const {
			bool ret;
			gt->runAA(1, local_anchors->get(curr_off),
					anchor.local_anchors->get(anchor.curr_off), &ret);
			return ret;
		}
	};

	std::priority_queue<anchor_struct> queue;
	size_t anchor_gap_size;
	size_t sort_buf_size;

	off_t get_anchor_off(const anchor_struct &anchor) const;
public:
	typedef std::shared_ptr<anchor_prio_queue> ptr;

	anchor_prio_queue(const std::vector<local_buf_vec_store::ptr> &anchor_vals);
	scalar_variable::ptr get_min_frontier() const;

	/*
	 * Here we pop a set of chunks of data whose values are the potentially
	 * the smallest.
	 */
	std::vector<off_t> pop(size_t size);
};

/*
 * This class keeps some summaries of sorted vectors in memory, so later on
 * we can use these summaries to pinpoint the right data for merging.
 * To generate a summary of a vector, we put some conceptional anchors in
 * the vector and keep the values in the anchor location.
 */
class sort_portion_summary
{
	size_t anchor_gap_size;
	size_t sort_buf_size;
	std::vector<local_buf_vec_store::ptr> anchor_vals;
public:
	sort_portion_summary(const scalar_type &type, size_t num_sort_bufs);
	void add_portion(local_buf_vec_store::ptr sorted_buf);
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
	off_t off_in_bytes;
	local_buf_vec_store::ptr store;
	// We have to make sure that the I/O is alive before all subvec
	// computation is destroyed.
	safs::io_interface &io;
	const EM_vec_store &vec;
	sort_portion_summary &summary;
public:
	EM_vec_sort_compute(local_buf_vec_store::ptr store,
			safs::io_interface &_io, const EM_vec_store &_vec,
			sort_portion_summary &_summary);
	virtual void run(char *buf, size_t size);
};

/*
 * This I/O task is to issue an I/O request to read the portion from
 * disks for sorting.
 */
class EM_vec_sort_issue_task: public portion_io_task
{
	local_buf_vec_store::ptr store;
	const EM_vec_store &vec;
	sort_portion_summary &summary;
public:
	EM_vec_sort_issue_task(off_t global_start, size_t length,
			const EM_vec_store &_vec,
			sort_portion_summary &_summary);
	virtual void run(safs::io_interface::ptr io1, safs::io_interface::ptr io2);
};

class EM_vec_merge_dispatcher;

/*
 * This I/O task is to issue an I/O request to read portions from
 * disks for merging.
 */
class EM_vec_merge_issue_task: public portion_io_task
{
	std::vector<local_buf_vec_store::ptr> stores;
	local_buf_vec_store::ptr prev_leftover;
	const EM_vec_store &vec;
	EM_vec_merge_dispatcher &dispatcher;
public:
	EM_vec_merge_issue_task(const std::vector<off_t> &_anchor_locs,
			local_buf_vec_store::ptr prev_leftover, const EM_vec_store &_vec,
			EM_vec_merge_dispatcher &dispatcher);

	virtual void run(safs::io_interface::ptr read_io,
			safs::io_interface::ptr write_io);
};

/*
 * This class merges data read from disks and writes it back to disks.
 */
class EM_vec_merge_compute: public portion_compute
{
	// The last buffer may be the leftover from the previous merge.
	std::vector<local_buf_vec_store::const_ptr> stores;
	// We have to make sure that the I/O is alive before all subvec
	// computation is destroyed.
	safs::io_interface &write_io;
	const EM_vec_store &vec;
	EM_vec_merge_dispatcher &dispatcher;
	size_t num_completed;
	size_t expected_ios;
public:
	EM_vec_merge_compute(const std::vector<local_buf_vec_store::ptr> &_stores,
			local_buf_vec_store::const_ptr prev_leftover,
			safs::io_interface &_io, const EM_vec_store &_vec,
			EM_vec_merge_dispatcher &_dispatcher);
	virtual void run(char *buf, size_t size);
};

}

}

}

#endif
