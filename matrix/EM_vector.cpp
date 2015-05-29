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

#include <libgen.h>
#include <malloc.h>

#include <unordered_map>

#include "safs_file.h"
#include "io_request.h"
#include "io_interface.h"
#include "in_mem_io.h"

#include "EM_vector.h"
#include "matrix_config.h"
#include "mem_worker_thread.h"
#include "local_vec_store.h"
#include "matrix_store.h"

namespace fm
{

namespace detail
{

static safs::file_io_factory::shared_ptr create_temp_file(size_t num_bytes)
{
	char *tmp = tempnam(".", "vec");
	std::string tmp_name = basename(tmp);
	safs::safs_file f(safs::get_sys_RAID_conf(), tmp_name);
	assert(!f.exist());
	bool ret = f.create_file(num_bytes);
	assert(ret);
	safs::file_io_factory::shared_ptr factory
		= safs::create_io_factory(tmp_name, safs::REMOTE_ACCESS);
	free(tmp);
	return factory;
}

EM_vec_store::EM_vec_store(size_t length, const scalar_type &type): vec_store(
		length, type, false)
{
	factory = create_temp_file(length * type.get_size());
}

EM_vec_store::~EM_vec_store()
{
	safs::safs_file f(safs::get_sys_RAID_conf(), factory->get_name());
	assert(f.exist());
	f.delete_file();
}

bool EM_vec_store::resize(size_t length)
{
	// TODO
	assert(0);
	return false;
}

bool EM_vec_store::append(
		std::vector<vec_store::const_ptr>::const_iterator vec_it,
		std::vector<vec_store::const_ptr>::const_iterator vec_end)
{
	assert(0);
}

bool EM_vec_store::append(const vec_store &vec)
{
	assert(0);
}

vec_store::ptr EM_vec_store::deep_copy() const
{
	assert(0);
}

vec_store::ptr EM_vec_store::shallow_copy()
{
	assert(0);
}

vec_store::const_ptr EM_vec_store::shallow_copy() const
{
	assert(0);
}

size_t EM_vec_store::get_portion_size() const
{
	return matrix_conf.get_anchor_gap_size() / get_entry_size();
}

void EM_vec_store::reset_data()
{
	assert(0);
}

vec_store::ptr EM_vec_store::sort_with_index()
{
	assert(0);
}

class EM_vec_dispatcher: public task_dispatcher
{
	const EM_vec_store &store;
	off_t portion_idx;
	pthread_spinlock_t lock;
	size_t portion_size;
public:
	EM_vec_dispatcher(const EM_vec_store &_store,
			size_t portion_size = 0): store(_store) {
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		portion_idx = 0;
		if (portion_size == 0)
			this->portion_size = store.get_portion_size();
		else
			this->portion_size = portion_size;
	}

	virtual portion_io_task::ptr get_task() {
		pthread_spin_lock(&lock);
		off_t global_start = portion_idx * portion_size;
		if ((size_t) global_start >= store.get_length()) {
			pthread_spin_unlock(&lock);
			return portion_io_task::ptr();
		}
		size_t length = std::min(portion_size, store.get_length() - global_start);
		portion_idx++;
		pthread_spin_unlock(&lock);
		return create_vec_task(global_start, length, store);
	}

	virtual portion_io_task::ptr create_vec_task(off_t global_start, size_t length,
			const EM_vec_store &store) = 0;
};

///////////////////////////// Sort the vector /////////////////////////////////

namespace EM_sort_detail
{

anchor_prio_queue::anchor_prio_queue(const std::vector<local_buf_vec_store::ptr> &anchor_vals)
{
	for (size_t i = 0; i < anchor_vals.size(); i++) {
		anchor_struct anchor;
		anchor.local_anchors = anchor_vals[i];
		anchor.id = i;
		anchor.curr_off = 0;
		const scalar_type &type = anchor_vals.front()->get_type();
		anchor.gt = type.get_basic_ops().get_op(basic_ops::op_idx::GT);;
		queue.push(anchor);
	}
	const scalar_type &type = anchor_vals.front()->get_type();
	anchor_gap_size = matrix_conf.get_anchor_gap_size() / type.get_size();
	sort_buf_size = matrix_conf.get_sort_buf_size() / type.get_size();
}

off_t anchor_prio_queue::get_anchor_off(const anchor_struct &anchor) const
{
	return anchor.id * sort_buf_size + anchor.curr_off * anchor_gap_size;
}

scalar_variable::ptr anchor_prio_queue::get_min_frontier() const
{
	if (queue.empty())
		return scalar_variable::ptr();
	else {
		local_buf_vec_store::const_ptr local_anchors = queue.top().local_anchors;
		const scalar_type &type = local_anchors->get_type();
		scalar_variable::ptr var = type.create_scalar();
		assert((size_t) queue.top().curr_off < local_anchors->get_length());
		var->set_raw(local_anchors->get(queue.top().curr_off),
				type.get_size());
		return var;
	}
}

/*
 * Here we pop a set of chunks of data whose values are the potentially
 * the smallest.
 */
std::vector<off_t> anchor_prio_queue::pop(size_t size)
{
	std::vector<off_t> chunks;
	long remaining_size = size;
	while (remaining_size > 0 && !queue.empty()) {
		anchor_struct anchor = queue.top();
		off_t off = get_anchor_off(anchor);
		chunks.push_back(off);
		remaining_size -= anchor_gap_size;
		queue.pop();

		// If there are still anchors left in the partition, we should
		// update it and put it back to the priority.
		anchor.curr_off++;
		if (anchor.local_anchors->get_length() > (size_t) anchor.curr_off)
			queue.push(anchor);
	}
	return chunks;
}

sort_portion_summary::sort_portion_summary(const scalar_type &type,
		size_t num_sort_bufs)
{
	size_t entry_size = type.get_size();
	anchor_gap_size = matrix_conf.get_anchor_gap_size() / entry_size;
	sort_buf_size = matrix_conf.get_sort_buf_size() / entry_size;
	anchor_vals.resize(num_sort_bufs);
}

void sort_portion_summary::add_portion(local_buf_vec_store::ptr sorted_buf)
{
	std::vector<off_t> idxs;
	for (size_t i = 0; i < sorted_buf->get_length(); i += anchor_gap_size)
		idxs.push_back(i);
	//		assert(idxs.back() < sorted_buf->get_length());
	//		if ((size_t) idxs.back() != sorted_buf->get_length() - 1)
	//			idxs.push_back(sorted_buf->get_length() - 1);
	//		assert(idxs.back() < sorted_buf->get_length());
	off_t idx = sorted_buf->get_global_start() / sort_buf_size;
	assert(anchor_vals[idx] == NULL);
	anchor_vals[idx] = sorted_buf->get(idxs);
	if ((size_t) idx == anchor_vals.size() - 1)
		assert(sorted_buf->get_length() <= sort_buf_size);
	else
		assert(sorted_buf->get_length() == sort_buf_size);
}

anchor_prio_queue::ptr sort_portion_summary::get_prio_queue() const
{
	return anchor_prio_queue::ptr(new anchor_prio_queue(anchor_vals));
}

EM_vec_sort_compute::EM_vec_sort_compute(local_buf_vec_store::ptr store,
			safs::io_interface &_io, const EM_vec_store &_vec,
			sort_portion_summary &_summary): io(_io), vec(_vec), summary(
				_summary)
{
	this->store = store;
	off_in_bytes = vec.get_byte_off(store->get_global_start());
}

void EM_vec_sort_compute::run(char *buf, size_t size)
{
	assert(store->get_raw_arr() == buf);
	assert(store->get_length() * store->get_entry_size() == size);
	// Sort each portion in parallel.
	// Here we rely on OpenMP to sort the data in the buffer in parallel.
	store->get_type().get_sorter().sort(buf, store->get_length(), false);
	summary.add_portion(store);

	// Write the sorting result to disks.
	safs::data_loc_t loc(io.get_file_id(), off_in_bytes);
	safs::io_request req(store->get_raw_arr(), loc,
			store->get_length() * store->get_entry_size(), WRITE);
	portion_compute::ptr compute(new portion_write_complete(store));
	register_portion_compute(req, compute);
	io.access(&req, 1);
	io.flush_requests();
}

EM_vec_sort_issue_task::EM_vec_sort_issue_task(off_t global_start, size_t length,
		const EM_vec_store &_vec,
		sort_portion_summary &_summary): vec(_vec), summary(_summary)
{
	store = local_buf_vec_store::ptr(new local_buf_vec_store(
				global_start, length, vec.get_type(), -1));
}

void EM_vec_sort_issue_task::run(safs::io_interface::ptr io1,
		safs::io_interface::ptr io2)
{
	// In this case, we read data from a file and write it back to
	// the same file.
	assert(io1 == io2);
	off_t off = vec.get_byte_off(store->get_global_start());
	safs::data_loc_t loc(io1->get_file_id(), off);
	safs::io_request req(store->get_raw_arr(), loc,
			store->get_length() * store->get_entry_size(), READ);
	portion_compute::ptr compute(new EM_vec_sort_compute(store, *io1, vec,
				summary));
	register_portion_compute(req, compute);
	io1->access(&req, 1);
	io1->flush_requests();
}

class EM_vec_sort_dispatcher: public EM_vec_dispatcher
{
	std::shared_ptr<sort_portion_summary> summary;
public:
	typedef std::shared_ptr<EM_vec_sort_dispatcher> ptr;

	EM_vec_sort_dispatcher(const EM_vec_store &store): EM_vec_dispatcher(store,
			// We want to have a larger buffer for sorting.
			matrix_conf.get_sort_buf_size() / store.get_entry_size()) {
		size_t sort_buf_size
			= matrix_conf.get_sort_buf_size() / store.get_entry_size();
		size_t num_sort_bufs
			= ceil(((double) store.get_length()) / sort_buf_size);
		summary = std::shared_ptr<sort_portion_summary>(
				new sort_portion_summary(store.get_type(), num_sort_bufs));
	}

	const sort_portion_summary &get_sort_summary() const {
		return *summary;
	}

	virtual portion_io_task::ptr create_vec_task(off_t global_start,
			size_t length, const EM_vec_store &store) {
		return portion_io_task::ptr(new EM_vec_sort_issue_task(global_start,
					length, store, *summary));
	}
};

/////////////////// Merge portions ///////////////////////

class merge_writer
{
	const size_t local_buf_size;	// In the number of elements
	off_t merge_end;		// In the number of bytes
	local_buf_vec_store::ptr buf;
	size_t data_size_in_buf;		// In the number of elements
public:
	merge_writer(const scalar_type &type): local_buf_size(
			matrix_conf.get_write_io_buf_size() / type.get_size()) {
		merge_end = 0;
		buf = local_buf_vec_store::ptr(new local_buf_vec_store(
					0, local_buf_size, type, -1));
		data_size_in_buf = 0;
	}

	void flush_buffer_data(safs::io_interface &write_io,
			io_worker_task &worker_task) {
		buf->resize(data_size_in_buf);
		const scalar_type &type = buf->get_type();
		safs::data_loc_t loc(write_io.get_file_id(), merge_end);
		safs::io_request req(buf->get_raw_arr(), loc,
				data_size_in_buf * buf->get_entry_size(), WRITE);
		portion_compute::ptr compute(new portion_write_complete(buf));
		worker_task.register_portion_compute(req, compute);
		write_io.access(&req, 1);
		write_io.flush_requests();
		merge_end += req.get_size();

		// The data is written asynchronously, we need to allocate a new buffer.
		buf = local_buf_vec_store::ptr(new local_buf_vec_store(
					0, local_buf_size, type, -1));
		data_size_in_buf = 0;
	}

	void append(safs::io_interface &write_io, local_vec_store::ptr data,
			io_worker_task &worker_task) {
		// In the number of elements.
		off_t off_in_new_data = 0;
		size_t new_data_size = data->get_length();

		while (new_data_size > 0) {
			size_t copy_data_size = std::min(new_data_size,
					// The size of the available space in the buffer.
					buf->get_length() - data_size_in_buf);
			// We always copy the data to the local buffer and write
			// the local buffer to disks. The reason of doing so is to
			// make sure the size and offset of data written to disks
			// are aligned to the I/O block size.
			// TODO maybe we should remove this extra memory copy.
			memcpy(buf->get(data_size_in_buf), data->get(off_in_new_data),
					copy_data_size * buf->get_entry_size());
			data_size_in_buf += copy_data_size;
			// Update the amount of data in the incoming data buffer.
			off_in_new_data += copy_data_size;
			new_data_size -= copy_data_size;

			// If the buffer is full, we need to write it out first.
			if (data_size_in_buf == buf->get_length())
				flush_buffer_data(write_io, worker_task);
		}
	}
};

class EM_vec_merge_dispatcher: public task_dispatcher
{
	const EM_vec_store &store;
	anchor_prio_queue::ptr anchors;
	local_buf_vec_store::ptr prev_leftover;
	size_t sort_buf_size;
	merge_writer writer;
public:
	EM_vec_merge_dispatcher(const EM_vec_store &_store,
			anchor_prio_queue::ptr anchors);
	void set_prev_leftover(local_buf_vec_store::ptr prev_leftover) {
		this->prev_leftover = prev_leftover;
	}

	merge_writer &get_merge_writer() {
		return writer;
	}

	const anchor_prio_queue &get_anchors() const {
		return *anchors;
	}

	virtual portion_io_task::ptr get_task();
};

EM_vec_merge_compute::EM_vec_merge_compute(
		const std::vector<local_buf_vec_store::ptr> &_stores,
		local_buf_vec_store::const_ptr prev_leftover,
		safs::io_interface &_io, const EM_vec_store &_vec,
		EM_vec_merge_dispatcher &_dispatcher): stores(_stores.begin(),
			_stores.end()), write_io(_io), vec(_vec), dispatcher(_dispatcher)
{
	expected_ios = this->stores.size();
	if (prev_leftover)
		this->stores.push_back(prev_leftover);
	num_completed = 0;
}

EM_vec_merge_dispatcher::EM_vec_merge_dispatcher(const EM_vec_store &_store,
		anchor_prio_queue::ptr anchors): store(_store), writer(store.get_type())
{
	this->anchors = anchors;
	sort_buf_size = matrix_conf.get_sort_buf_size() / store.get_entry_size();
}

portion_io_task::ptr EM_vec_merge_dispatcher::get_task()
{
	size_t leftover = 0;
	if (prev_leftover)
		leftover = prev_leftover->get_length();
	assert(sort_buf_size > leftover);
	const std::vector<off_t> anchor_locs = anchors->pop(
			sort_buf_size - leftover);
	if (anchor_locs.empty() && prev_leftover == NULL) {
		assert(anchors->get_min_frontier() == NULL);
		assert(leftover == 0);
		return portion_io_task::ptr();
	}
	else {
		return portion_io_task::ptr(new EM_vec_merge_issue_task(anchor_locs,
					prev_leftover, store, *this));
	}
}

EM_vec_merge_issue_task::EM_vec_merge_issue_task(
		const std::vector<off_t> &_anchor_locs,
		local_buf_vec_store::ptr prev_leftover,
		const EM_vec_store &_vec, EM_vec_merge_dispatcher &_dispatcher): vec(
			_vec), dispatcher(_dispatcher)
{
	this->prev_leftover = prev_leftover;
	std::vector<size_t> anchor_locs(_anchor_locs.begin(), _anchor_locs.end());
	std::sort(anchor_locs.begin(), anchor_locs.end());
	size_t anchor_gap_size
		= matrix_conf.get_anchor_gap_size() / vec.get_type().get_size();
	for (size_t i = 0; i < anchor_locs.size(); i++) {
		size_t num_eles = std::min(anchor_gap_size,
				vec.get_length() - anchor_locs[i]);
		size_t off = anchor_locs[i];
		// If the anchors are contiguous, we merge them.
		while (i + 1 < anchor_locs.size()
				&& anchor_locs[i + 1] == anchor_locs[i] + anchor_gap_size) {
			i++;
			num_eles += std::min(anchor_gap_size,
					vec.get_length() - anchor_locs[i]);
		}
		stores.push_back(local_buf_vec_store::ptr(new local_buf_vec_store(
						off, num_eles, vec.get_type(), -1)));
	}
}

void EM_vec_merge_issue_task::run(safs::io_interface::ptr read_io,
		safs::io_interface::ptr write_io)
{
	// In this case, we need to read some data from the disks first and then
	// merge with the data left from the previous merge.
	if (!stores.empty()) {
		portion_compute::ptr compute(new EM_vec_merge_compute(stores,
					prev_leftover, *write_io, vec, dispatcher));
		safs::io_request reqs[stores.size()];
		for (size_t i = 0; i < stores.size(); i++) {
			off_t off = vec.get_byte_off(stores[i]->get_global_start());
			safs::data_loc_t loc(read_io->get_file_id(), off);
			reqs[i] = safs::io_request(stores[i]->get_raw_arr(), loc,
					stores[i]->get_length() * stores[i]->get_entry_size(), READ);
			register_portion_compute(reqs[i], compute);
		}
		read_io->access(reqs, stores.size());
		read_io->flush_requests();
	}
	// In this case, we don't need to read data from disks any more.
	// We only need to write the data left from the previous merge.
	else {
		dispatcher.get_merge_writer().append(*write_io, prev_leftover,
				get_worker_task());
		// This is the last write. we should flush everything to disks.
		dispatcher.get_merge_writer().flush_buffer_data(*write_io,
				get_worker_task());
		// No more leftover.
		dispatcher.set_prev_leftover(NULL);
	}
}

void EM_vec_merge_compute::run(char *buf, size_t size)
{
	num_completed++;
	// If all data in the buffers is ready, we should merge all the buffers.
	if (num_completed == expected_ios) {
		// Find the min values among the last elements in the buffers.
		const scalar_type &type = vec.get_type();
		scalar_variable::ptr min_val
			= dispatcher.get_anchors().get_min_frontier();

		// Breaks the local buffers into two parts. The first part is to
		// merge with others; we have to keep the second part for further
		// merging.
		std::vector<std::pair<const char *, const char *> > merge_data(
				stores.size());
		std::vector<std::pair<const char *, const char *> > leftovers(
				stores.size());
		size_t leftover_size = 0;
		size_t merge_size = 0;
		for (size_t i = 0; i < stores.size(); i++) {
			const char *start = stores[i]->get_raw_arr();
			const char *end = stores[i]->get_raw_arr()
				+ stores[i]->get_length() * stores[i]->get_entry_size();
			off_t leftover_start;
			if (min_val != NULL)
				leftover_start = type.get_stl_algs().lower_bound(
						start, end, min_val->get_raw());
			else
				leftover_start = stores[i]->get_length();
			merge_size += leftover_start;
			leftover_size += (stores[i]->get_length() - leftover_start);
			merge_data[i] = std::pair<const char *, const char *>(
					stores[i]->get(0), stores[i]->get(leftover_start));
			leftovers[i] = std::pair<const char *, const char *>(
					stores[i]->get(leftover_start),
					stores[i]->get(stores[i]->get_length()));
		}

		// Here we rely on OpenMP to merge the data in the buffer in parallel.
		local_buf_vec_store::ptr merge_res(new local_buf_vec_store(0,
					merge_size, type, -1));
		type.get_sorter().merge(merge_data, merge_res->get_raw_arr(),
				merge_size);

		// Write the merge result to disks.
		dispatcher.get_merge_writer().append(write_io, merge_res,
				get_worker_task());
		// Keep the leftover and merge them into a single buffer.
		local_buf_vec_store::ptr leftover_buf = local_buf_vec_store::ptr(
				new local_buf_vec_store(0, leftover_size, type, -1));
		type.get_sorter().merge(leftovers, leftover_buf->get_raw_arr(),
				leftover_size);
		dispatcher.set_prev_leftover(leftover_buf);
	}
}

}

void EM_vec_store::sort()
{
	assert(matrix_conf.get_sort_buf_size() % get_entry_size() == 0);
	size_t sort_buf_size = matrix_conf.get_sort_buf_size() / get_entry_size();
	size_t portion_size = get_portion_size();
	assert(sort_buf_size >= portion_size);
	assert(sort_buf_size % portion_size == 0);
	assert(matrix_conf.get_anchor_gap_size() % get_entry_size() == 0);
	size_t anchor_gap_size = matrix_conf.get_anchor_gap_size() / get_entry_size();
	assert(anchor_gap_size >= portion_size);
	assert(anchor_gap_size % portion_size == 0);

	/*
	 * Divide the vector into multiple large parts and sort each part in parallel.
	 */
	EM_sort_detail::EM_vec_sort_dispatcher::ptr sort_dispatcher(
			new EM_sort_detail::EM_vec_sort_dispatcher(*this));
	io_worker_task sort_worker(factory, factory, sort_dispatcher, 2);
	sort_worker.run();

	/* Merge all parts.
	 * Here we assume that one level of merging is enough and we rely on
	 * OpenMP to parallelize merging.
	 */
	EM_sort_detail::EM_vec_merge_dispatcher::ptr merge_dispatcher(
			new EM_sort_detail::EM_vec_merge_dispatcher(*this,
				sort_dispatcher->get_sort_summary().get_prio_queue()));
	safs::file_io_factory::shared_ptr new_factory = create_temp_file(
			get_length() * get_entry_size());
	// TODO let's not use asynchornous I/O for now.
	io_worker_task merge_worker(factory, new_factory, merge_dispatcher, 0);
	merge_worker.run();

	// In the end, we points to the new file.
	factory = new_factory;
}

////////////////////////// Set data of the vector ////////////////////////////

namespace
{

/*
 * This class records the summary of the array. It stores three pieces of
 * information from each portion: the first and last elements in the portion
 * and a flag that indicates whether the portion is sorted.
 * This implementation assumes that two levels to test whether an array
 * is sorted. It is enough for a very large vector.
 * Although this data structure is shared by all threads, each portion owns
 * its own elements in the shared vector, so locking isn't needed.
 */
class issorted_summary
{
	// This stores the elements in the each end of a portion.
	local_buf_vec_store::ptr ends;
	// This vector is modified by multiple threads in parallel, so we can't
	// use a boolean vector, which isn't thread-safe. We should use
	// thread-safe bitmap.
	std::vector<int> issorted;
	size_t portion_size;
public:
	issorted_summary(const EM_vec_store &vec) {
		size_t num_portions = vec.get_num_portions();
		ends = local_vec_store::ptr(new local_buf_vec_store(0,
					num_portions * 2, vec.get_type(), -1));
		issorted.resize(num_portions);
		portion_size = vec.get_portion_size();
	}

	void set_portion_result(local_buf_vec_store::ptr store) {
		bool sorted = store->get_type().get_sorter().is_sorted(
				store->get_raw_arr(), store->get_length(), false);
		off_t portion_idx = store->get_global_start() / portion_size;
		assert(portion_idx >= 0 && (size_t) portion_idx < issorted.size());
		issorted[portion_idx] = sorted;
		// Save the elements in each end.
		assert((size_t) portion_idx * 2 + 1 < ends->get_length());
		ends->set_raw(portion_idx * 2, store->get(0));
		ends->set_raw(portion_idx * 2 + 1, store->get(store->get_length() - 1));
	}

	bool is_sorted() const {
		for (size_t i = 0; i < issorted.size(); i++)
			if (!issorted[i])
				return false;
		bool ret = ends->get_type().get_sorter().is_sorted(
				ends->get_raw_arr(), ends->get_length(), false);
		return ret;
	}
};

class EM_vec_issorted_compute: public portion_compute
{
	local_buf_vec_store::ptr store;
	issorted_summary &summary;
public:
	EM_vec_issorted_compute(local_buf_vec_store::ptr store,
			issorted_summary &_summary): summary(_summary) {
		this->store = store;
	}

	virtual void run(char *buf, size_t size) {
		assert(store->get_raw_arr() == buf);
		assert(store->get_length() * store->get_entry_size() == size);
		summary.set_portion_result(store);
	}
};

/*
 * This I/O task is to issue an I/O request to read the portion from
 * disks for testing if the array is sorted.
 */
class EM_vec_issorted_issue_task: public portion_io_task
{
	off_t off_in_bytes;
	local_buf_vec_store::ptr store;
	issorted_summary &summary;
public:
	EM_vec_issorted_issue_task(off_t global_start, size_t length,
			const EM_vec_store &vec,
			issorted_summary &_summary): summary(_summary) {
		store = local_buf_vec_store::ptr(new local_buf_vec_store(
					global_start, length, vec.get_type(), -1));
		off_in_bytes = vec.get_byte_off(store->get_global_start());
	}
	virtual void run(safs::io_interface::ptr read_io, safs::io_interface::ptr) {
		safs::data_loc_t loc(read_io->get_file_id(), off_in_bytes);
		safs::io_request req(store->get_raw_arr(), loc,
				store->get_length() * store->get_entry_size(), READ);
		portion_compute::ptr compute(new EM_vec_issorted_compute(store, summary));
		register_portion_compute(req, compute);
		read_io->access(&req, 1);
		read_io->flush_requests();
	}
};

}

class EM_vec_issorted_dispatcher: public EM_vec_dispatcher
{
	issorted_summary summary;
public:
	typedef std::shared_ptr<EM_vec_issorted_dispatcher> ptr;

	EM_vec_issorted_dispatcher(const EM_vec_store &store): EM_vec_dispatcher(
			store), summary(store) {
	}

	const issorted_summary &get_summary() const {
		return summary;
	}

	virtual portion_io_task::ptr create_vec_task(off_t global_start,
			size_t length, const EM_vec_store &store) {
		return portion_io_task::ptr(new EM_vec_issorted_issue_task(global_start,
					length, store, summary));
	}
};

bool EM_vec_store::is_sorted() const
{
	mem_thread_pool::ptr threads = mem_thread_pool::get_global_mem_threads();
	EM_vec_issorted_dispatcher::ptr dispatcher(
			new EM_vec_issorted_dispatcher(*this));
	for (size_t i = 0; i < threads->get_num_threads(); i++) {
		io_worker_task *task = new io_worker_task(factory, NULL, dispatcher);
		threads->process_task(i % threads->get_num_nodes(), task);
	}
	threads->wait4complete();
	return dispatcher->get_summary().is_sorted();
}

////////////////////////// Set data of the vector ////////////////////////////

namespace
{

class EM_vec_setdata_task: public portion_io_task
{
	off_t off_in_bytes;
	local_buf_vec_store::ptr buf;
public:
	EM_vec_setdata_task(off_t global_start, size_t length,
			const EM_vec_store &vec, const set_operate &op) {
		buf = local_buf_vec_store::ptr(new local_buf_vec_store(
					global_start, length, vec.get_type(), -1));
		off_in_bytes = vec.get_byte_off(buf->get_global_start());
		buf->set_data(op);
	}

	virtual void run(safs::io_interface::ptr, safs::io_interface::ptr write_io) {
		safs::data_loc_t loc(write_io->get_file_id(), off_in_bytes);
		safs::io_request req(buf->get_raw_arr(), loc,
				buf->get_length() * buf->get_entry_size(), WRITE);
		portion_compute::ptr compute(new portion_write_complete(buf));
		register_portion_compute(req, compute);
		write_io->access(&req, 1);
		write_io->flush_requests();
	}
};

class EM_vec_setdata_dispatcher: public EM_vec_dispatcher
{
	const set_operate &op;
public:
	EM_vec_setdata_dispatcher(const EM_vec_store &store,
			const set_operate &_op): EM_vec_dispatcher(store), op(_op) {
	}

	virtual portion_io_task::ptr create_vec_task(off_t global_start,
			size_t length, const EM_vec_store &store) {
		return portion_io_task::ptr(new EM_vec_setdata_task(global_start,
					length, store, op));
	}
};

}

void EM_vec_store::set_data(const set_operate &op)
{
	mem_thread_pool::ptr threads = mem_thread_pool::get_global_mem_threads();
	EM_vec_setdata_dispatcher::ptr dispatcher(
			new EM_vec_setdata_dispatcher(*this, op));
	for (size_t i = 0; i < threads->get_num_threads(); i++) {
		io_worker_task *task = new io_worker_task(NULL, factory, dispatcher);
		threads->process_task(i % threads->get_num_nodes(), task);
	}
	threads->wait4complete();
}

matrix_store::const_ptr EM_vec_store::conv2mat(size_t nrow, size_t ncol,
			bool byrow) const
{
	BOOST_LOG_TRIVIAL(error)
		<< "can't convert a NUMA vector to a matrix";
	return matrix_store::ptr();
}

}

}
