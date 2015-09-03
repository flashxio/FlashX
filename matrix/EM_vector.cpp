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
#include <boost/math/common_factor.hpp>
#include <boost/format.hpp>

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

namespace
{
/*
 * When we write data to disks, we need to have something to hold the buffer.
 * This holds the local buffer until the write completes.
 */
class portion_write_complete: public portion_compute
{
	local_buf_vec_store::const_ptr store;
public:
	portion_write_complete(local_buf_vec_store::const_ptr store) {
		this->store = store;
	}

	virtual void run(char *buf, size_t size) {
		assert(store->get_raw_arr() == buf);
		assert(store->get_length() * store->get_entry_size() == size);
	}
};

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

	virtual bool issue_task() {
		pthread_spin_lock(&lock);
		off_t global_start = portion_idx * portion_size;
		if ((size_t) global_start >= store.get_length()) {
			pthread_spin_unlock(&lock);
			return false;
		}
		size_t length = std::min(portion_size, store.get_length() - global_start);
		portion_idx++;
		pthread_spin_unlock(&lock);
		create_vec_task(global_start, length, store);
		return true;
	}

	virtual void create_vec_task(off_t global_start, size_t length,
			const EM_vec_store &from_vec) = 0;
};

class seq_writer
{
	size_t local_buf_size;	// In the number of elements
	off_t merge_end;		// In the number of elements
	local_buf_vec_store::ptr buf;
	size_t data_size_in_buf;		// In the number of elements
	EM_vec_store::ptr to_vec;
public:
	seq_writer(EM_vec_store::ptr vec, off_t write_start) {
		this->local_buf_size
			= matrix_conf.get_write_io_buf_size() / vec->get_type().get_size();
		this->to_vec = vec;
		merge_end = write_start;
		buf = local_buf_vec_store::ptr(new local_buf_vec_store(
					-1, local_buf_size, to_vec->get_type(), -1));
		data_size_in_buf = 0;
	}

	~seq_writer() {
		if (data_size_in_buf > 0)
			flush_buffer_data(true);
	}

	void flush_buffer_data(bool last);
	void append(local_vec_store::const_ptr data);
};

void seq_writer::flush_buffer_data(bool last)
{
	if (!last)
		assert((data_size_in_buf * buf->get_length()) % PAGE_SIZE == 0);
	else
		data_size_in_buf = ROUNDUP((data_size_in_buf * buf->get_entry_size()),
				PAGE_SIZE) / buf->get_entry_size();
	buf->resize(data_size_in_buf);
	const scalar_type &type = buf->get_type();
	to_vec->write_portion_async(buf, merge_end);
	merge_end += data_size_in_buf;

	if (!last)
		// The data is written asynchronously, we need to allocate
		// a new buffer.
		buf = local_buf_vec_store::ptr(new local_buf_vec_store(
					-1, local_buf_size, type, -1));
	else
		buf = NULL;
	data_size_in_buf = 0;
}

void seq_writer::append(local_vec_store::const_ptr data)
{
	assert(buf);
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
			flush_buffer_data(false);
	}
}

EM_vec_store::ptr EM_vec_store::cast(vec_store::ptr vec)
{
	if (vec->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error) << "Can't cast an in-mem vector to EM_vec_store";
		return EM_vec_store::ptr();
	}
	return std::static_pointer_cast<EM_vec_store>(vec);
}

EM_vec_store::const_ptr EM_vec_store::cast(vec_store::const_ptr vec)
{
	if (vec->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error) << "Can't cast an in-mem vector to EM_vec_store";
		return EM_vec_store::const_ptr();
	}
	return std::static_pointer_cast<const EM_vec_store>(vec);
}

EM_vec_store::EM_vec_store(safs::file_io_factory::shared_ptr factory): vec_store(
		// Without giving any information, we assume this is a byte array.
		factory->get_file_size(), get_scalar_type<char>(), false)
{
	holder = file_holder::create(factory->get_name());
	ios = io_set::ptr(new io_set(factory));
}

EM_vec_store::EM_vec_store(const EM_vec_store &store): vec_store(
		store.get_length(), store.get_type(), false)
{
	holder = store.holder;
	ios = store.ios;
}

EM_vec_store::EM_vec_store(size_t length, const scalar_type &type): vec_store(
		length, type, false)
{
	holder = file_holder::create_temp("vec", length * type.get_size());
	safs::file_io_factory::shared_ptr factory = safs::create_io_factory(
			holder->get_name(), safs::REMOTE_ACCESS);
	ios = io_set::ptr(new io_set(factory));
}

EM_vec_store::~EM_vec_store()
{
}

bool EM_vec_store::resize(size_t length)
{
	// TODO
	assert(0);
	return false;
}

namespace
{

class EM_vec_append_dispatcher: public task_dispatcher
{
	std::vector<vec_store::const_ptr>::const_iterator vec_it;
	std::vector<vec_store::const_ptr>::const_iterator vec_end;
	seq_writer &writer;
	size_t portion_size;
	off_t local_off;
public:
	EM_vec_append_dispatcher(seq_writer &_writer,
			std::vector<vec_store::const_ptr>::const_iterator vec_start,
			std::vector<vec_store::const_ptr>::const_iterator vec_end,
			size_t portion_size): writer(_writer) {
		this->vec_it = vec_start;
		this->vec_end = vec_end;
		this->portion_size = portion_size;
		this->local_off = 0;
	}
	virtual bool issue_task();
};

bool EM_vec_append_dispatcher::issue_task()
{
	if (vec_it == vec_end) {
		writer.flush_buffer_data(true);
		return false;
	}
	vec_store::const_ptr vec = *vec_it;
	size_t size;
	off_t local_curr_off = local_off;
	if (portion_size >= vec->get_length() - local_off) {
		size = vec->get_length() - local_off;
		vec_it++;
		local_off = 0;
	}
	else {
		size = portion_size;
		local_off += portion_size;
	}
	// TODO we might want to read portion asynchronously.
	local_vec_store::const_ptr lstore = vec->get_portion(local_curr_off, size);
	writer.append(lstore);
	return true;
}

}

bool EM_vec_store::append(
		std::vector<vec_store::const_ptr>::const_iterator vec_start,
		std::vector<vec_store::const_ptr>::const_iterator vec_end)
{
	size_t tot_size = 0;
	for (auto it = vec_start; it != vec_end; it++) {
		tot_size += (*it)->get_length();
		if (get_type() != (*it)->get_type()) {
			BOOST_LOG_TRIVIAL(error)
				<< "can't append a vector with different type";
			return false;
		}
	}

	/*
	 * If the last page that stores the elements in the vector isn't full,
	 * we need to read the last page first.
	 */

	// In the number of bytes.
	off_t off = ROUND(get_length() * get_entry_size(), PAGE_SIZE);
	size_t size = get_length() * get_entry_size() - off;
	assert(off % get_entry_size() == 0);
	seq_writer writer(EM_vec_store::ptr(this, empty_free()),
			off / get_entry_size());
	if (size > 0) {
		local_vec_store::ptr portion = get_portion(off / get_entry_size(),
				size / get_entry_size());
		writer.append(portion);
	}

	size_t portion_size = matrix_conf.get_stream_io_size() / get_type().get_size();
	task_dispatcher::ptr dispatcher(new EM_vec_append_dispatcher(writer,
				vec_start, vec_end, portion_size));
	io_worker_task worker(dispatcher, 1);
	worker.register_EM_obj(this);
	worker.run();

	return vec_store::resize(get_length() + tot_size);
}

bool EM_vec_store::append(const vec_store &vec)
{
	struct deleter {
		void operator()(const vec_store *) {
		}
	};
	std::vector<vec_store::const_ptr> vecs(1);
	vecs[0] = vec_store::const_ptr(&vec, deleter());
	return append(vecs.begin(), vecs.end());
}

namespace
{

class EM_vec_copy_write: public portion_compute
{
	local_buf_vec_store::const_ptr store;
	EM_vec_store &to_vec;
public:
	EM_vec_copy_write(EM_vec_store &_to_vec): to_vec(_to_vec) {
	}

	void set_buf(local_buf_vec_store::const_ptr store) {
		this->store = store;
	}

	virtual void run(char *buf, size_t size) {
		assert(store);
		assert(store->get_raw_arr() == buf);
		assert(store->get_length() * store->get_entry_size() == size);
		to_vec.write_portion_async(store);
	}
};

class EM_vec_copy_dispatcher: public EM_vec_dispatcher
{
	EM_vec_store &to_vec;
public:
	EM_vec_copy_dispatcher(const EM_vec_store &from_store,
			EM_vec_store &_to_vec, size_t portion_size): EM_vec_dispatcher(
				from_store, portion_size), to_vec(_to_vec) {
	}

	virtual void create_vec_task(off_t global_start, size_t orig_length,
			const EM_vec_store &from_vec) {
		// If the length of the vector isn't aligned with the page size,
		// the underlying file is extended a little.
		size_t entry_size = from_vec.get_type().get_size();
		assert(round_ele(global_start, PAGE_SIZE, entry_size) == global_start);
		size_t length = roundup_ele(orig_length, PAGE_SIZE, entry_size);
		EM_vec_copy_write *compute = new EM_vec_copy_write(to_vec);
		local_buf_vec_store::ptr buf = from_vec.get_portion_async(
				global_start, length, portion_compute::ptr(compute));
		compute->set_buf(buf);
	}
};

}

vec_store::ptr EM_vec_store::deep_copy() const
{
	// TODO we might need to give users the optional to config it.
	const size_t copy_portion_size = 128 * 1024 * 1024;
	EM_vec_store::ptr new_vec = EM_vec_store::create(get_length(), get_type());
	EM_vec_copy_dispatcher::ptr copy_dispatcher(
			new EM_vec_copy_dispatcher(*this, *new_vec, copy_portion_size));
	// The buffer size is large, we don't need some many async I/Os.
	io_worker_task sort_worker(copy_dispatcher, 1);
	sort_worker.register_EM_obj(const_cast<EM_vec_store *>(this));
	sort_worker.register_EM_obj(new_vec.get());
	sort_worker.run();
	return new_vec;
}

vec_store::ptr EM_vec_store::shallow_copy()
{
	return vec_store::ptr(new EM_vec_store(*this));
}

vec_store::const_ptr EM_vec_store::shallow_copy() const
{
	return vec_store::ptr(new EM_vec_store(*this));
}

size_t EM_vec_store::get_portion_size() const
{
	// TODO
	return 1024 * 1024;
}

local_vec_store::ptr EM_vec_store::get_portion_async(off_t orig_start,
		size_t orig_size, portion_compute::ptr compute) const
{
	size_t entry_size = get_type().get_size();
	off_t start = round_ele(orig_start, PAGE_SIZE, entry_size);
	size_t size = orig_size + (orig_start - start);
	size = roundup_ele(size, PAGE_SIZE, entry_size);

	// TODO fix the bug if `start' and `size' aren't aligned with PAGE_SIZE.
	safs::io_interface &io = ios->get_curr_io();
	local_buf_vec_store::ptr buf(new local_buf_vec_store(start, size,
				get_type(), -1));
	off_t off = start * entry_size;
	safs::data_loc_t loc(io.get_file_id(), off);
	safs::io_request req(buf->get_raw_arr(), loc,
			buf->get_length() * buf->get_entry_size(), READ);
	static_cast<portion_callback &>(io.get_callback()).add(req, compute);
	io.access(&req, 1);
	io.flush_requests();

	if (orig_start != start || orig_size != size)
		buf->expose_sub_vec(orig_start - start, orig_size);
	return buf;
}

std::vector<local_vec_store::ptr> EM_vec_store::get_portion_async(
		const std::vector<std::pair<off_t, size_t> > &locs,
		portion_compute::ptr compute) const
{
	size_t entry_size = get_type().get_size();
	// TODO fix the bug if `locs' aren't aligned with PAGE_SIZE.
	safs::io_interface &io = ios->get_curr_io();
	std::vector<safs::io_request> reqs(locs.size());
	std::vector<local_vec_store::ptr> ret_bufs(locs.size());
	for (size_t i = 0; i < locs.size(); i++) {
		off_t start = round_ele(locs[i].first, PAGE_SIZE, entry_size);
		size_t size = locs[i].second + (locs[i].first - start);
		size = roundup_ele(size, PAGE_SIZE, entry_size);

		local_buf_vec_store::ptr buf(new local_buf_vec_store(start, size,
					get_type(), -1));
		off_t off = start * entry_size;
		safs::data_loc_t loc(io.get_file_id(), off);
		reqs[i] = safs::io_request(buf->get_raw_arr(), loc,
				buf->get_length() * buf->get_entry_size(), READ);
		static_cast<portion_callback &>(io.get_callback()).add(reqs[i], compute);

		if (locs[i].first != start || locs[i].second != size)
			buf->expose_sub_vec(locs[i].first - start, locs[i].second);
		ret_bufs[i] = buf;
	}
	io.access(reqs.data(), reqs.size());
	io.flush_requests();
	return ret_bufs;
}

bool EM_vec_store::set_portion(std::shared_ptr<const local_vec_store> store,
		off_t loc)
{
	if (store->get_type() != get_type()) {
		BOOST_LOG_TRIVIAL(error) << "The input store has a different type";
		return false;
	}
	if (loc + store->get_length() > get_length()) {
		BOOST_LOG_TRIVIAL(error) << "out of boundary";
		return false;
	}

	write_portion_async(store, loc);
	safs::io_interface &io = ios->get_curr_io();
	io.wait4complete(1);
	return true;
}

local_vec_store::const_ptr EM_vec_store::get_portion(off_t loc, size_t size) const
{
	return const_cast<EM_vec_store *>(this)->get_portion(loc, size);
}

local_vec_store::ptr EM_vec_store::get_portion(off_t orig_loc, size_t orig_size)
{
	if (orig_loc + orig_size > get_length()) {
		BOOST_LOG_TRIVIAL(error) << "Out of boundary";
		return local_vec_store::ptr();
	}

	size_t entry_size = get_type().get_size();
	off_t loc = round_ele(orig_loc, PAGE_SIZE, entry_size);
	size_t size = orig_size + (orig_loc - loc);
	size = roundup_ele(size, PAGE_SIZE, entry_size);
	safs::io_interface &io = ios->get_curr_io();
	bool ready = false;
	portion_compute::ptr compute(new sync_read_compute(ready));
	local_vec_store::ptr ret = get_portion_async(loc, size, compute);
	while (!ready)
		io.wait4complete(1);
	if (orig_loc != loc || orig_size != size)
		ret->expose_sub_vec(orig_loc - loc, orig_size);
	return ret;
}

void EM_vec_store::write_portion_async(local_vec_store::const_ptr store,
		off_t off)
{
	off_t start = off;
	if (start < 0)
		start = store->get_global_start();
	assert(start >= 0);

	safs::io_interface &io = ios->get_curr_io();
	off_t off_in_bytes = start * get_type().get_size();
	safs::data_loc_t loc(io.get_file_id(), off_in_bytes);
	safs::io_request req(const_cast<char *>(store->get_raw_arr()), loc,
			store->get_length() * store->get_entry_size(), WRITE);
	portion_compute::ptr compute(new portion_write_complete(store));
	static_cast<portion_callback &>(io.get_callback()).add(req, compute);
	io.access(&req, 1);
	// TODO I might want to flush requests later.
	io.flush_requests();
}

void EM_vec_store::reset_data()
{
	assert(0);
}

vec_store::ptr EM_vec_store::sort_with_index()
{
	assert(0);
}

std::vector<safs::io_interface::ptr> EM_vec_store::create_ios() const
{
	std::vector<safs::io_interface::ptr> ret(1);
	ret[0] = ios->create_io();
	return ret;
}

///////////////////////////// Sort the vector /////////////////////////////////

namespace EM_sort_detail
{

anchor_prio_queue::anchor_prio_queue(
		const std::vector<local_buf_vec_store::ptr> &anchor_vals,
		size_t _sort_buf_size, size_t _anchor_gap_size): sort_buf_size(
			_sort_buf_size), anchor_gap_size(_anchor_gap_size), queue(
			anchor_ptr_less(anchor_vals.front()->get_type()))
{
	anchor_bufs.resize(anchor_vals.size());
	for (size_t i = 0; i < anchor_vals.size(); i++) {
		anchor_struct anchor;
		anchor.local_anchors = anchor_vals[i];
		anchor.id = i;
		anchor.curr_off = 0;
		anchor_bufs[i] = anchor;
		queue.push(&anchor_bufs[i]);
	}
}

off_t anchor_prio_queue::get_anchor_off(const anchor_struct &anchor) const
{
	return anchor.id * sort_buf_size + anchor.curr_off * anchor_gap_size;
}

/*
 * By looking into the values in the anchor locations, we can know immediately
 * the minimal value among the data that hasn't been read.
 */
scalar_variable::ptr anchor_prio_queue::get_min_frontier() const
{
	if (queue.empty())
		return scalar_variable::ptr();
	else {
		local_buf_vec_store::const_ptr local_anchors = queue.top()->local_anchors;
		const scalar_type &type = local_anchors->get_type();
		scalar_variable::ptr var = type.create_scalar();
		assert((size_t) queue.top()->curr_off < local_anchors->get_length());
		var->set_raw(local_anchors->get(queue.top()->curr_off),
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
		anchor_struct *anchor = queue.top();
		assert((size_t) anchor->curr_off < anchor->local_anchors->get_length());
		off_t off = get_anchor_off(*anchor);
		chunks.push_back(off);
		remaining_size -= anchor_gap_size;
		queue.pop();

		// If there are still anchors left in the partition, we should
		// update it and put it back to the priority.
		anchor->curr_off++;
		if (anchor->local_anchors->get_length() > (size_t) anchor->curr_off)
			queue.push(anchor);
	}
	return chunks;
}

std::vector<off_t> anchor_prio_queue::fetch_all_first()
{
	std::vector<off_t> chunks;
	for (size_t i = 0; i < anchor_bufs.size(); i++) {
		anchor_struct &anchor = anchor_bufs[i];
		if (anchor.local_anchors->get_length() > (size_t) anchor.curr_off) {
			off_t off = get_anchor_off(anchor);
			chunks.push_back(off);
			anchor.curr_off++;
		}
	}

	// We have change the anchor structs, now we have reconstruct
	// the priority queue again.
	queue = anchor_queue_t(anchor_ptr_less(
				anchor_bufs.front().local_anchors->get_type()));
	assert(queue.empty());
	for (size_t i = 0; i < anchor_bufs.size(); i++) {
		anchor_struct *anchor = &anchor_bufs[i];
		if (anchor->local_anchors->get_length() > (size_t) anchor->curr_off)
			queue.push(anchor);
	}

	return chunks;
}

sort_portion_summary::sort_portion_summary(size_t num_sort_bufs,
		size_t _sort_buf_size, size_t _anchor_gap_size): sort_buf_size(
			_sort_buf_size), anchor_gap_size(_anchor_gap_size)
{
	anchor_vals.resize(num_sort_bufs);
}

void sort_portion_summary::add_portion(local_buf_vec_store::const_ptr sorted_buf)
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
	return anchor_prio_queue::ptr(new anchor_prio_queue(anchor_vals,
				sort_buf_size, anchor_gap_size));
}

void EM_vec_sort_compute::run(char *buf, size_t size)
{
	num_completed++;
	if (num_completed == portions.size()) {
		// Sort each portion in parallel.
		// Here we rely on OpenMP to sort the data in the buffer in parallel.
		local_buf_vec_store::ptr sort_buf = portions.front();
		std::vector<off_t> orig_offs(sort_buf->get_length());
		sort_buf->get_type().get_sorter().sort_with_index(
				sort_buf->get_raw_arr(), orig_offs.data(),
				sort_buf->get_length(), false);
		summary.add_portion(sort_buf);

		// It might be a sub vector, we should reset its exposed part, so
		// the size of data written to disks is aligned to the page size.
		sort_buf->reset_expose();
		// Write the sorting result to disks.
		to_vecs.front()->write_portion_async(sort_buf);
		for (size_t i = 1; i < portions.size(); i++) {
			portions[i]->reset_expose();
			// If the element size is different in each array, the padding
			// size might also be different. We should make the elements
			// in the padding area are mapped to the same locations.
			size_t off = orig_offs.size();
			orig_offs.resize(portions[i]->get_length());
			if (off < orig_offs.size()) {
				for (; off < orig_offs.size(); off++)
					orig_offs[off] = off;
			}

			local_vec_store::ptr shuffle_buf = portions[i]->get(orig_offs);
			to_vecs[i]->write_portion_async(shuffle_buf,
					portions[i]->get_global_start());
		}
	}
}

class EM_vec_sort_dispatcher: public EM_vec_dispatcher
{
	std::shared_ptr<sort_portion_summary> summary;
	std::vector<EM_vec_store::const_ptr> from_vecs;
	std::vector<EM_vec_store::ptr> to_vecs;
public:
	typedef std::shared_ptr<EM_vec_sort_dispatcher> ptr;

	EM_vec_sort_dispatcher(const std::vector<EM_vec_store::const_ptr> &from_vecs,
			const std::vector<EM_vec_store::ptr> &to_vecs,
			size_t sort_buf_size, size_t anchor_gap_size);

	const sort_portion_summary &get_sort_summary() const {
		return *summary;
	}

	virtual void create_vec_task(off_t global_start,
			size_t length, const EM_vec_store &from_vec);
};

EM_vec_sort_dispatcher::EM_vec_sort_dispatcher(
		const std::vector<EM_vec_store::const_ptr> &from_vecs,
		const std::vector<EM_vec_store::ptr> &to_vecs,
		size_t sort_buf_size, size_t anchor_gap_size): EM_vec_dispatcher(
			*from_vecs.front(), sort_buf_size)
{
	EM_vec_store::const_ptr sort_vec = from_vecs.front();
	size_t num_sort_bufs
		= ceil(((double) sort_vec->get_length()) / sort_buf_size);
	summary = std::shared_ptr<sort_portion_summary>(
			new sort_portion_summary(num_sort_bufs, sort_buf_size,
				anchor_gap_size));
	this->from_vecs = from_vecs;
	this->to_vecs = to_vecs;
}

void EM_vec_sort_dispatcher::create_vec_task(off_t global_start,
			size_t orig_length, const EM_vec_store &from_vec)
{
	EM_vec_sort_compute *sort_compute = new EM_vec_sort_compute(to_vecs,
			*summary);
	portion_compute::ptr compute(sort_compute);
	std::vector<local_vec_store::ptr> from_portions(from_vecs.size());
	for (size_t i = 0; i < from_portions.size(); i++) {
		size_t entry_size = from_vecs[i]->get_type().get_size();
		assert(round_ele(global_start, PAGE_SIZE, entry_size) == global_start);
		size_t length = roundup_ele(orig_length, PAGE_SIZE, entry_size);

		from_portions[i] = from_vecs[i]->get_portion_async(
				global_start, length, compute);
		from_portions[i]->expose_sub_vec(0, orig_length);
	}
	sort_compute->set_bufs(from_portions);
}

/////////////////// Merge portions ///////////////////////

class EM_vec_merge_dispatcher: public task_dispatcher
{
	std::vector<EM_vec_store::const_ptr> from_vecs;
	std::vector<local_buf_vec_store::ptr> prev_leftovers;
	anchor_prio_queue::ptr anchors;
	size_t sort_buf_size;
	std::vector<seq_writer> writers;
public:
	EM_vec_merge_dispatcher(const std::vector<EM_vec_store::const_ptr> &from_vecs,
			const std::vector<EM_vec_store::ptr> &to_vecs,
			anchor_prio_queue::ptr anchors, size_t sort_buf_size);
	void set_prev_leftovers(
			const std::vector<local_buf_vec_store::ptr> &prev_leftovers) {
		this->prev_leftovers = prev_leftovers;
	}

	local_buf_vec_store::const_ptr get_prev_leftover(off_t idx) const {
		return prev_leftovers[idx];
	}

	seq_writer &get_merge_writer(int idx) {
		return writers[idx];
	}

	const anchor_prio_queue &get_anchors() const {
		return *anchors;
	}

	virtual bool issue_task();
};

EM_vec_merge_compute::EM_vec_merge_compute(
		const std::vector<local_buf_vec_store::ptr> &prev_leftovers,
		EM_vec_merge_dispatcher &_dispatcher): dispatcher(_dispatcher)
{
	stores.resize(prev_leftovers.size());
	for (size_t i = 0; i < prev_leftovers.size(); i++) {
		// If there is a leftover for a vector from the previous merge,
		// all vectors should have the same number of leftover elements.
		if (prev_leftovers[0]) {
			assert(prev_leftovers[i]);
			assert(prev_leftovers[0]->get_length()
					== prev_leftovers[i]->get_length());
		}
		if (prev_leftovers[i])
			this->stores[i].push_back(prev_leftovers[i]);
	}
	num_completed = 0;
}

void EM_vec_merge_compute::set_bufs(const std::vector<merge_set_t> &bufs)
{
	assert(bufs.size() == stores.size());
	num_expected = 0;
	for (size_t i = 0; i < bufs.size(); i++) {
		// If all vectors should have the same number of buffers to merge.
		assert(bufs[0].size() == bufs[i].size());
		num_expected += bufs[i].size();
		this->stores[i].insert(this->stores[i].end(), bufs[i].begin(),
				bufs[i].end());
	}
}

EM_vec_merge_dispatcher::EM_vec_merge_dispatcher(
		const std::vector<EM_vec_store::const_ptr> &from_vecs,
		const std::vector<EM_vec_store::ptr> &to_vecs,
		anchor_prio_queue::ptr anchors, size_t _sort_buf_size): sort_buf_size(
			_sort_buf_size)
{
	this->from_vecs = from_vecs;
	assert(from_vecs.size() == to_vecs.size());
	for (size_t i = 0; i < from_vecs.size(); i++)
		assert(from_vecs[i]->get_type() == to_vecs[i]->get_type());
	this->anchors = anchors;
	for (size_t i = 0; i < to_vecs.size(); i++)
		writers.emplace_back(to_vecs[i], 0);
	prev_leftovers.resize(from_vecs.size());
}

bool EM_vec_merge_dispatcher::issue_task()
{
	typedef std::vector<local_buf_vec_store::const_ptr> merge_set_t;
	size_t leftover = 0;
	assert(!prev_leftovers.empty());
	std::vector<off_t> anchor_locs;
	size_t anchor_gap_size = anchors->get_anchor_gap_size();
	if (prev_leftovers[0]) {
		leftover = prev_leftovers[0]->get_length();
		if (sort_buf_size <= leftover)
			BOOST_LOG_TRIVIAL(info)
				<< boost::format("leftover (%1%) is larger than sort buf size (%2%)")
				% leftover % sort_buf_size;
		size_t pop_size;
		if (sort_buf_size > leftover)
			pop_size = sort_buf_size - leftover;
		else
			pop_size = anchor_gap_size;
		anchor_locs = anchors->pop(pop_size);
	}
	else {
		anchor_locs = anchors->fetch_all_first();
		size_t fetch_size = anchor_locs.size() * anchor_gap_size;
		if (fetch_size < sort_buf_size) {
			std::vector<off_t> more_locs = anchors->pop(
					sort_buf_size - fetch_size);
			anchor_locs.insert(anchor_locs.end(), more_locs.begin(),
					more_locs.end());
		}
	}
	// If there isn't any data to merge and there isn't leftover from
	// the previous merge.
	if (anchor_locs.empty() && prev_leftovers[0] == NULL) {
		assert(anchors->get_min_frontier() == NULL);
		assert(leftover == 0);
		// If there is nothing left to merge and there isn't leftover, we still
		// need to flush the buffered data.
		for (size_t i = 0; i < writers.size(); i++)
			writers[i].flush_buffer_data(true);
		return false;
	}
	else {
		// Merge the anchors.
		std::vector<std::pair<off_t, size_t> > data_locs;
		std::sort(anchor_locs.begin(), anchor_locs.end());
		for (size_t i = 0; i < anchor_locs.size(); i++) {
			size_t num_eles = std::min(anchor_gap_size,
					from_vecs[0]->get_length() - anchor_locs[i]);
			size_t off = anchor_locs[i];
			// If the anchors are contiguous, we merge them.
			while (i + 1 < anchor_locs.size()
					&& (size_t) anchor_locs[i + 1] == anchor_locs[i] + anchor_gap_size) {
				i++;
				num_eles += std::min(anchor_gap_size,
						from_vecs[0]->get_length() - anchor_locs[i]);
			}
			data_locs.push_back(std::pair<off_t, size_t>(off, num_eles));
		}

		// In this case, we need to read some data from the disks first and then
		// merge with the data left from the previous merge.
		if (!data_locs.empty()) {
			EM_vec_merge_compute *_compute = new EM_vec_merge_compute(
					prev_leftovers, *this);
			portion_compute::ptr compute(_compute);
			std::vector<merge_set_t> merge_sets(from_vecs.size());
			for (size_t j = 0; j < from_vecs.size(); j++) {
				std::vector<local_vec_store::ptr> portions
					= from_vecs[j]->get_portion_async(data_locs, compute);
				merge_sets[j].insert(merge_sets[j].end(), portions.begin(),
						portions.end());
			}
			_compute->set_bufs(merge_sets);
		}
		// In this case, we don't need to read data from disks any more.
		// We only need to write the data left from the previous merge.
		else {
			for (size_t i = 0; i < writers.size(); i++) {
				if (prev_leftovers[i])
					writers[i].append(prev_leftovers[i]);
				// This is the last write. we should flush everything to disks.
				writers[i].flush_buffer_data(true);
				// No more leftover.
				prev_leftovers[i] = NULL;
			}
		}
		return true;
	}
}

void EM_vec_merge_compute::run(char *buf, size_t size)
{
	num_completed++;
	// If all data in the buffers is ready, we should merge all the buffers.
	if (num_completed == num_expected) {
		assert(stores.size() > 0);
		merge_set_t &merge_bufs = stores[0];
		// Find the min values among the last elements in the buffers.
		const scalar_type &type = merge_bufs.front()->get_type();
		scalar_variable::ptr min_val
			= dispatcher.get_anchors().get_min_frontier();

		// Breaks the local buffers into two parts. The first part is to
		// merge with others; we have to keep the second part for further
		// merging.
		std::vector<std::pair<const char *, const char *> > merge_data;
		std::vector<std::pair<const char *, const char *> > leftovers;
		std::vector<size_t> merge_sizes(merge_bufs.size());
		size_t leftover_size = 0;
		size_t merge_size = 0;
		local_buf_vec_store::const_ptr prev_leftover
			= dispatcher.get_prev_leftover(0);
		// We go through all the buffers to be merged and merge elements
		// that are smaller than `min_val' and keep all elements in the `leftover'
		// buffer, which have been read from the disks but are larger than
		// `min_val'.
		for (size_t i = 0; i < merge_bufs.size(); i++) {
			size_t entry_size = merge_bufs[i]->get_entry_size();
			const size_t tot_len = merge_bufs[i]->get_length();
			const char *start = merge_bufs[i]->get_raw_arr();
			const char *end = merge_bufs[i]->get_raw_arr()
				+ tot_len * entry_size;
			off_t leftover_start;
			if (min_val != NULL) {
				leftover_start = type.get_stl_algs().lower_bound(
						start, end, min_val->get_raw());
				// lower_bound finds the location so that all elements before
				// the location have values smaller than `min_val'. Actually,
				// we can also merge all elements whose value is equal to
				// `min_val'.
				if ((size_t) leftover_start < tot_len && min_val->equals(start
							+ leftover_start * entry_size)) {
					size_t rel_loc;
					type.get_agg_ops().get_find_next().run(
							tot_len - leftover_start,
							start + leftover_start * entry_size, &rel_loc);
					// There is at least one element with the same value as
					// `min_val'.
					assert(rel_loc > 0 && rel_loc <= tot_len - leftover_start);
					leftover_start += rel_loc;
				}
				assert((size_t) leftover_start <= tot_len);
			}
			else
				leftover_start = tot_len;
			assert((size_t) leftover_start <= tot_len);
			merge_sizes[i] = leftover_start;
			merge_size += leftover_start;
			leftover_size += (tot_len - leftover_start);
			if (leftover_start > 0)
				merge_data.push_back(std::pair<const char *, const char *>(
							merge_bufs[i]->get(0),
							merge_bufs[i]->get(leftover_start)));
			if (tot_len - leftover_start)
				leftovers.push_back(std::pair<const char *, const char *>(
							merge_bufs[i]->get(leftover_start),
							merge_bufs[i]->get(tot_len)));
		}

		// Here we rely on OpenMP to merge the data in the buffer in parallel.
		local_buf_vec_store::ptr merge_res(new local_buf_vec_store(-1,
					merge_size, type, -1));
		std::vector<std::pair<int, off_t> > merge_index(merge_size);
		type.get_sorter().merge_with_index(merge_data,
				merge_res->get_raw_arr(), merge_size, merge_index);
		// Write the merge result to disks.
		dispatcher.get_merge_writer(0).append(merge_res);
		merge_res = NULL;

		std::vector<std::pair<int, off_t> > leftover_merge_index(leftover_size);
		std::vector<local_buf_vec_store::ptr> leftover_bufs(stores.size());
		if (leftover_size > 0) {
			// Keep the leftover and merge them into a single buffer.
			local_buf_vec_store::ptr leftover_buf = local_buf_vec_store::ptr(
					new local_buf_vec_store(-1, leftover_size, type, -1));
			type.get_sorter().merge_with_index(leftovers,
					leftover_buf->get_raw_arr(), leftover_size,
					leftover_merge_index);
			leftover_bufs[0] = leftover_buf;
		}

		// Merge the remaining vectors accordingly.
		for (size_t i = 1; i < stores.size(); i++) {
			std::vector<std::pair<const char *, const char *> > merge_data;
			std::vector<std::pair<const char *, const char *> > leftovers;

			merge_set_t &set = stores[i];
			assert(set.size() == merge_bufs.size());
			for (size_t i = 0; i < set.size(); i++) {
				off_t leftover_start = merge_sizes[i];
				assert(set[i]->get_length() == merge_bufs[i]->get_length());
				if (leftover_start > 0)
					merge_data.push_back(std::pair<const char *, const char *>(
								set[i]->get(0), set[i]->get(leftover_start)));
				if (set[i]->get_length() - leftover_start)
					leftovers.push_back(std::pair<const char *, const char *>(
								set[i]->get(leftover_start),
								set[i]->get(set[i]->get_length())));
			}

			// Merge the part that can be merged.
			const scalar_type &type = set.front()->get_type();
			merge_res = local_buf_vec_store::ptr(new local_buf_vec_store(-1,
						merge_size, type, -1));
			type.get_sorter().merge(merge_data, merge_index,
					merge_res->get_raw_arr(), merge_size);
			dispatcher.get_merge_writer(i).append(merge_res);

			if (leftover_size > 0) {
				// Keep the leftover and merge them into a single buffer.
				local_buf_vec_store::ptr leftover_buf = local_buf_vec_store::ptr(
						new local_buf_vec_store(-1, leftover_size, type, -1));
				type.get_sorter().merge(leftovers, leftover_merge_index,
						leftover_buf->get_raw_arr(), leftover_size);
				leftover_bufs[i] = leftover_buf;
			}
		}

		dispatcher.set_prev_leftovers(leftover_bufs);
	}
}

/*
 * The two functions compute the sort buffer size and the anchor gap size
 * based on the sort buffer size and min I/O size provided by the user.
 * One version works for a single vector and the other version works for
 * multiple vectors.
 *
 * The anchor gap size should be multiple of the entry sizes of all
 * vectors as well as multiple of the min I/O size.
 * The sort buffer size should be multiple of the anchor gap size.
 */

std::pair<size_t, size_t> cal_sort_buf_size(const scalar_type &type,
		size_t num_eles)
{
	size_t min_anchor_gap_bytes = 1;		// in the number of bytes.
	min_anchor_gap_bytes = boost::math::lcm(min_anchor_gap_bytes,
			type.get_size());
	min_anchor_gap_bytes = boost::math::lcm(min_anchor_gap_bytes,
			(size_t) PAGE_SIZE);
	assert(min_anchor_gap_bytes % type.get_size() == 0);
	assert(min_anchor_gap_bytes % PAGE_SIZE == 0);

	// The number of elements between two anchors.
	size_t min_anchor_gap_size = min_anchor_gap_bytes / type.get_size();
	size_t num_anchors
		= matrix_conf.get_sort_buf_size() / type.get_size() / min_anchor_gap_size;
	size_t sort_buf_size = num_anchors * min_anchor_gap_size;
	assert((sort_buf_size * type.get_size()) % PAGE_SIZE == 0);

	// Find the maximal size of anchor gap size allowed for the sort buffer
	// size and the vector length.
	assert(sort_buf_size % min_anchor_gap_size == 0);
	size_t num_sort_bufs = ceil(((double) num_eles) / sort_buf_size);
	size_t num_min_anchors = sort_buf_size / min_anchor_gap_size;
	size_t factor = sort_buf_size / num_sort_bufs / min_anchor_gap_size;
	for (; factor > 0; factor--)
		if (num_min_anchors % factor == 0)
			break;
	assert(factor != 0);
	return std::pair<size_t, size_t>(sort_buf_size, min_anchor_gap_size * factor);
}

std::pair<size_t, size_t> cal_sort_buf_size(
		const std::vector<const scalar_type *> &types, size_t num_eles)
{
	size_t tot_entry_size = 0;
	size_t min_anchor_gap_bytes = 1;		// in the number of bytes.
	for (size_t i = 0; i < types.size(); i++) {
		tot_entry_size += types[i]->get_size();
		min_anchor_gap_bytes = boost::math::lcm(min_anchor_gap_bytes,
				types[i]->get_size());
	}
	min_anchor_gap_bytes = boost::math::lcm(min_anchor_gap_bytes,
			(size_t) PAGE_SIZE);
	for (size_t i = 0; i < types.size(); i++) {
		assert(min_anchor_gap_bytes % types[i]->get_size() == 0);
	}
	assert(min_anchor_gap_bytes % PAGE_SIZE == 0);

	// The number of elements between two anchors.
	size_t min_anchor_gap_size = min_anchor_gap_bytes / types[0]->get_size();
	size_t num_anchors
		= matrix_conf.get_sort_buf_size() / tot_entry_size / min_anchor_gap_size;
	size_t sort_buf_size = num_anchors * min_anchor_gap_size;
	for (size_t i = 0; i < types.size(); i++)
		assert((sort_buf_size * types[i]->get_size()) % PAGE_SIZE == 0);

	// Find the maximal size of anchor gap size allowed for the sort buffer
	// size and the vector length.
	assert(sort_buf_size % min_anchor_gap_size == 0);
	size_t num_sort_bufs = ceil(((double) num_eles) / sort_buf_size);
	size_t num_min_anchors = sort_buf_size / min_anchor_gap_size;
	size_t factor = sort_buf_size / num_sort_bufs / min_anchor_gap_size;
	for (; factor > 0; factor--)
		if (num_min_anchors % factor == 0)
			break;
	assert(factor != 0);
	return std::pair<size_t, size_t>(sort_buf_size,
			min_anchor_gap_size * factor);
}

}

std::vector<EM_vec_store::ptr> sort(
		const std::vector<EM_vec_store::const_ptr> &vecs)
{
	assert(vecs.size() > 0);
	for (size_t i = 1; i < vecs.size(); i++) {
		if (vecs[i]->get_length() != vecs[0]->get_length()) {
			BOOST_LOG_TRIVIAL(error) << "Not all vectors have the same length";
			return std::vector<EM_vec_store::ptr>();
		}
	}

	std::vector<const scalar_type *> types(vecs.size());
	for (size_t i = 0; i < vecs.size(); i++)
		types[i] = &vecs[i]->get_type();
	std::pair<size_t, size_t> sizes = EM_sort_detail::cal_sort_buf_size(types,
			vecs.front()->get_length());
	size_t sort_buf_size = sizes.first;
	size_t anchor_gap_size = sizes.second;
	printf("sort buf size: %ld, anchor gap size: %ld\n", sort_buf_size,
			anchor_gap_size);
	for (size_t i = 0; i < vecs.size(); i++) {
		size_t num_sort_bufs
			= ceil(((double) vecs[i]->get_length()) / sort_buf_size);
		// We have to make sure the sort buffer can contain a anchor portion
		// from each partially sorted buffer.
		assert(num_sort_bufs * anchor_gap_size <= sort_buf_size);
	}

	/*
	 * Divide the vector into multiple large parts and sort each part in parallel.
	 */
	std::vector<EM_vec_store::ptr> tmp_vecs(vecs.size());
	for (size_t i = 0; i < vecs.size(); i++)
		tmp_vecs[i] = EM_vec_store::create(vecs[i]->get_length(),
				vecs[i]->get_type());
	EM_sort_detail::EM_vec_sort_dispatcher::ptr sort_dispatcher(
			new EM_sort_detail::EM_vec_sort_dispatcher(vecs, tmp_vecs,
				sort_buf_size, anchor_gap_size));
	io_worker_task sort_worker(sort_dispatcher, 1);
	for (size_t i = 0; i < vecs.size(); i++) {
		sort_worker.register_EM_obj(const_cast<EM_vec_store *>(vecs[i].get()));
		sort_worker.register_EM_obj(tmp_vecs[i].get());
	}
	sort_worker.run();

	/* Merge all parts.
	 * Here we assume that one level of merging is enough and we rely on
	 * OpenMP to parallelize merging.
	 */
	std::vector<EM_vec_store::ptr> out_vecs(vecs.size());
	for (size_t i = 0; i < vecs.size(); i++)
		out_vecs[i] = EM_vec_store::create(vecs[i]->get_length(),
				vecs[i]->get_type());
	std::vector<EM_vec_store::const_ptr> tmp_vecs1(tmp_vecs.begin(),
			tmp_vecs.end());
	EM_sort_detail::EM_vec_merge_dispatcher::ptr merge_dispatcher(
			new EM_sort_detail::EM_vec_merge_dispatcher(tmp_vecs1, out_vecs,
				sort_dispatcher->get_sort_summary().get_prio_queue(),
				sort_buf_size));
	// TODO let's not use asynchornous I/O for now.
	io_worker_task merge_worker(merge_dispatcher, 0);
	for (size_t i = 0; i < vecs.size(); i++) {
		merge_worker.register_EM_obj(tmp_vecs[i].get());
		merge_worker.register_EM_obj(out_vecs[i].get());
	}
	merge_worker.run();
	return out_vecs;
}

void EM_vec_store::sort()
{
	std::pair<size_t, size_t> sizes
		= EM_sort_detail::cal_sort_buf_size(get_type(), get_length());
	size_t sort_buf_size = sizes.first;
	size_t anchor_gap_size = sizes.second;
	size_t num_sort_bufs = ceil(((double) get_length()) / sort_buf_size);
	printf("sort buf size: %ld, anchor gap size: %ld\n", sort_buf_size,
			anchor_gap_size);
	// We have to make sure the sort buffer can contain a anchor portion
	// from each partially sorted buffer.
	assert(num_sort_bufs * anchor_gap_size <= sort_buf_size);

	/*
	 * Divide the vector into multiple large parts and sort each part in parallel.
	 */
	std::vector<EM_vec_store::const_ptr> in_vecs(1);
	std::vector<EM_vec_store::ptr> out_vecs(1);
	in_vecs[0] = EM_vec_store::const_ptr(this, empty_free());
	out_vecs[0] = EM_vec_store::ptr(this, empty_free());
	EM_sort_detail::EM_vec_sort_dispatcher::ptr sort_dispatcher(
			new EM_sort_detail::EM_vec_sort_dispatcher(in_vecs, out_vecs,
				sort_buf_size, anchor_gap_size));
	io_worker_task sort_worker(sort_dispatcher, 1);
	sort_worker.register_EM_obj(this);
	sort_worker.run();

	/* Merge all parts.
	 * Here we assume that one level of merging is enough and we rely on
	 * OpenMP to parallelize merging.
	 */
	EM_vec_store::ptr tmp = EM_vec_store::create(get_length(), get_type());
	in_vecs[0] = EM_vec_store::const_ptr(this, empty_free());
	out_vecs[0] = tmp;
	EM_sort_detail::EM_vec_merge_dispatcher::ptr merge_dispatcher(
			new EM_sort_detail::EM_vec_merge_dispatcher(in_vecs, out_vecs,
				sort_dispatcher->get_sort_summary().get_prio_queue(),
				sort_buf_size));
	// TODO let's not use asynchornous I/O for now.
	io_worker_task merge_worker(merge_dispatcher, 0);
	merge_worker.register_EM_obj(this);
	merge_worker.register_EM_obj(tmp.get());
	merge_worker.run();

	// In the end, we points to the new file.
	holder = tmp->holder;
	ios = tmp->ios;
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
		ends = local_vec_store::ptr(new local_buf_vec_store(-1,
					num_portions * 2, vec.get_type(), -1));
		issorted.resize(num_portions);
		portion_size = vec.get_portion_size();
	}

	void set_portion_result(local_buf_vec_store::const_ptr store) {
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
	local_buf_vec_store::const_ptr store;
	issorted_summary &summary;
public:
	EM_vec_issorted_compute(issorted_summary &_summary): summary(_summary) {
	}

	void set_buf(local_buf_vec_store::const_ptr store) {
		this->store = store;
	}

	virtual void run(char *buf, size_t size) {
		assert(store->get_raw_arr() == buf);
		summary.set_portion_result(store);
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

	virtual void create_vec_task(off_t global_start,
			size_t length, const EM_vec_store &from_vec);
};

void EM_vec_issorted_dispatcher::create_vec_task(off_t global_start,
			size_t orig_length, const EM_vec_store &from_vec)
{
	// The vector doesn't necessary have the number of elements to make
	// the end of the vector aligned with the page size. If it's not,
	// we should align it. We end up writing more data to the underlying
	// file.
	size_t entry_size = from_vec.get_type().get_size();
	assert(round_ele(global_start, PAGE_SIZE, entry_size) == global_start);
	size_t length = roundup_ele(orig_length, PAGE_SIZE, entry_size);
	EM_vec_issorted_compute *compute = new EM_vec_issorted_compute(summary);
	local_vec_store::const_ptr portion = from_vec.get_portion_async(global_start,
			length, portion_compute::ptr(compute));
	if (length != orig_length)
		const_cast<local_vec_store &>(*portion).expose_sub_vec(0, orig_length);
	compute->set_buf(portion);
}

bool EM_vec_store::is_sorted() const
{
	mem_thread_pool::ptr threads = mem_thread_pool::get_global_mem_threads();
	EM_vec_issorted_dispatcher::ptr dispatcher(
			new EM_vec_issorted_dispatcher(*this));
	for (size_t i = 0; i < threads->get_num_threads(); i++) {
		io_worker_task *task = new io_worker_task(dispatcher);
		task->register_EM_obj(const_cast<EM_vec_store *>(this));
		threads->process_task(i % threads->get_num_nodes(), task);
	}
	threads->wait4complete();
	return dispatcher->get_summary().is_sorted();
}

////////////////////////// Set data of the vector ////////////////////////////

namespace
{

class EM_vec_setdata_dispatcher: public EM_vec_dispatcher
{
	const set_vec_operate &op;
	EM_vec_store &to_vec;
public:
	EM_vec_setdata_dispatcher(EM_vec_store &store,
			const set_vec_operate &_op): EM_vec_dispatcher(store), op(
				_op), to_vec(store) {
	}

	virtual void create_vec_task(off_t global_start,
			size_t length, const EM_vec_store &from_vec);
};

void EM_vec_setdata_dispatcher::create_vec_task(off_t global_start,
		size_t orig_length, const EM_vec_store &from_vec)
{
	// The vector doesn't necessary have the number of elements to make
	// the end of the vector aligned with the page size. If it's not,
	// we should align it. We end up writing more data to the underlying
	// file.
	size_t entry_size = from_vec.get_type().get_size();
	assert(round_ele(global_start, PAGE_SIZE, entry_size) == global_start);
	size_t length = roundup_ele(orig_length, PAGE_SIZE, entry_size);
	local_buf_vec_store::ptr buf(new local_buf_vec_store(
				global_start, length, to_vec.get_type(), -1));
	if (length != orig_length)
		buf->expose_sub_vec(0, orig_length);
	buf->set_data(op);
	if (length != orig_length)
		buf->reset_expose();
	to_vec.write_portion_async(buf);
}

}

void EM_vec_store::set_data(const set_vec_operate &op)
{
	mem_thread_pool::ptr threads = mem_thread_pool::get_global_mem_threads();
	EM_vec_setdata_dispatcher::ptr dispatcher(
			new EM_vec_setdata_dispatcher(*this, op));
	for (size_t i = 0; i < threads->get_num_threads(); i++) {
		io_worker_task *task = new io_worker_task(dispatcher);
		task->register_EM_obj(this);
		threads->process_task(i % threads->get_num_nodes(), task);
	}
	threads->wait4complete();
}

matrix_store::const_ptr EM_vec_store::conv2mat(size_t nrow, size_t ncol,
			bool byrow) const
{
	BOOST_LOG_TRIVIAL(error)
		<< "can't convert an EM vector to a matrix";
	return matrix_store::ptr();
}

bool EM_vec_store::set_persistent(const std::string &name)
{
	if (!holder->set_persistent(name))
		return false;
	// TODO we have to make sure no other threads are accessing the data
	// in the vector. How can we do that?
	safs::file_io_factory::shared_ptr factory = safs::create_io_factory(
			holder->get_name(), safs::REMOTE_ACCESS);
	ios = io_set::ptr(new io_set(factory));
	return true;
}

}

}
