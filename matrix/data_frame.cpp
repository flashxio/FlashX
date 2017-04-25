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

#include <boost/format.hpp>

#include "log.h"

#include "matrix_config.h"
#include "data_frame.h"
#include "bulk_operate.h"
#include "mem_vec_store.h"
#include "EM_vector.h"
#include "mem_vv_store.h"
#include "vector_vector.h"
#include "EM_vv_store.h"
#include "bulk_operate_ext.h"

namespace fm
{

data_frame::data_frame(const std::vector<named_vec_t> &named_vecs)
{
	assert(!named_vecs.empty());
	this->named_vecs = named_vecs;
	for (auto it = named_vecs.begin(); it != named_vecs.end(); it++)
		vec_map.insert(*it);
}

bool data_frame::is_in_mem() const
{
	// If one of the vector is on disk, this data frame is considered
	// external-memory.
	bool in_mem = true;
	for (size_t i = 0; i < named_vecs.size(); i++)
		in_mem &= named_vecs[i].second->is_in_mem();
	return in_mem;
}

data_frame::const_ptr data_frame::shuffle_vecs(
		const std::vector<off_t> &vec_idxs) const
{
	std::vector<named_vec_t> sub_vecs(vec_idxs.size());
	for (size_t i = 0; i < sub_vecs.size(); i++) {
		size_t idx = vec_idxs[i];
		if (idx >= get_num_vecs()) {
			BOOST_LOG_TRIVIAL(error) << "vec idx is out of bound";
			return data_frame::const_ptr();
		}
		sub_vecs[i] = named_vecs[idx];
	}
	return data_frame::create(sub_vecs);
}

bool data_frame::append(std::vector<data_frame::ptr>::const_iterator begin,
		std::vector<data_frame::ptr>::const_iterator end)
{
	std::unordered_map<std::string, std::vector<detail::vec_store::const_ptr> > vecs;
	for (size_t i = 0; i < named_vecs.size(); i++) {
		std::vector<detail::vec_store::const_ptr> _vecs;
		_vecs.push_back(get_vec(named_vecs[i].first));
		vecs.insert(std::pair<std::string, std::vector<detail::vec_store::const_ptr> >(
					named_vecs[i].first, _vecs));
	}

	for (auto it = begin; it != end; it++) {
		data_frame::ptr df = *it;
		if (df->get_num_vecs() != get_num_vecs()) {
			BOOST_LOG_TRIVIAL(error)
				<< "The data frames have different numbers of vectors";
			return false;
		}
		for (auto vec_it = vecs.begin(); vec_it != vecs.end(); vec_it++) {
			detail::vec_store::const_ptr vec = df->get_vec(vec_it->first);
			if (vec == NULL) {
				BOOST_LOG_TRIVIAL(error)
					<< "The data frames have different names on the vectors";
				return false;
			}
			if (vec->get_type() != vec_it->second.front()->get_type()) {
				BOOST_LOG_TRIVIAL(error)
					<< "The data frames have different types for the vectors with the same name";
				return false;
			}
			vec_it->second.push_back(vec);
		}
	}

	for (auto it = vecs.begin(); it != vecs.end(); it++) {
		detail::vec_store::ptr vec = get_vec(it->first);
		vec->append(it->second.begin() + 1, it->second.end());
	}
	return true;
}

bool data_frame::append(data_frame::ptr df)
{
	for (auto it = named_vecs.begin(); it != named_vecs.end(); it++) {
		if (df->get_vec(it->first) == NULL) {
			BOOST_LOG_TRIVIAL(error)
				<< boost::format("The new data frame doesn't have column %1%")
				% it->first;
			return false;
		}
	}

	for (auto it = named_vecs.begin(); it != named_vecs.end(); it++)
		it->second->append(*df->get_vec(it->first));
	return true;
}

data_frame::const_ptr data_frame::sort(const std::string &col_name) const
{
	detail::vec_store::const_ptr sorted_col = get_vec(col_name);
	if (sorted_col == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"The column %1% doesn't exist") % col_name;
		return data_frame::const_ptr();
	}
	if (sorted_col->is_sorted())
		return this->shallow_copy();

	data_frame::ptr ret(new data_frame());
	if (is_in_mem()) {
		detail::vec_store::ptr copy_col = sorted_col->deep_copy();
		detail::smp_vec_store::ptr idxs = detail::smp_vec_store::cast(
				copy_col->sort_with_index());
		for (size_t i = 0; i < get_num_vecs(); i++) {
			detail::smp_vec_store::const_ptr mem_vec
				= detail::smp_vec_store::cast(get_vec(i));
			if (mem_vec == sorted_col) {
				ret->add_vec(col_name, copy_col);
			}
			else {
				detail::mem_vec_store::ptr tmp = mem_vec->get(*idxs);
				ret->add_vec(get_vec_name(i), tmp);
			}
		}
	}
	else {
		std::vector<std::string> names;
		std::vector<detail::EM_vec_store::const_ptr> vecs;
		names.push_back(col_name);
		assert(!sorted_col->is_in_mem());
		vecs.push_back(detail::EM_vec_store::cast(sorted_col));
		for (size_t i = 0; i < named_vecs.size(); i++) {
			if (named_vecs[i].second != sorted_col) {
				assert(!named_vecs[i].second->is_in_mem());
				vecs.push_back(detail::EM_vec_store::cast(named_vecs[i].second));
				names.push_back(named_vecs[i].first);
			}
		}
		std::vector<detail::EM_vec_store::ptr> sorted = detail::sort(vecs);
		// We should reshuffle the columns so that the columns in the returned
		// data frame have the same order as the current data frame.
		size_t j = 1;
		for (size_t i = 0; i < named_vecs.size(); i++) {
			if (named_vecs[i].first == col_name)
				ret->add_vec(col_name, sorted[0]);
			else {
				assert(names[j] == named_vecs[i].first);
				ret->add_vec(names[j], sorted[j]);
				j++;
			}
		}
	}
	return ret;
}

bool data_frame::is_sorted(const std::string &col_name) const
{
	return get_vec_ref(col_name).is_sorted();
}

data_frame::const_ptr data_frame::shallow_copy() const
{
	return data_frame::const_ptr(new data_frame(*this));
}

bool data_frame::add_vec(const std::string &name, detail::vec_store::ptr vec)
{
	if (get_num_vecs() > 0) {
		if (vec->get_length() != get_num_entries()) {
			BOOST_LOG_TRIVIAL(error)
				<< "Add a vector with different number of entries from the data frame";
			return false;
		}
	}
	named_vecs.push_back(named_vec_t(name, vec));
	vec_map.insert(named_vec_t(name, vec));
	return true;
}

data_frame::ptr merge_data_frame(const std::vector<data_frame::const_ptr> &dfs,
		bool in_mem)
{
	assert(!dfs.empty());
	size_t num_vecs = dfs[0]->get_num_vecs();
	for (size_t i = 1; i < dfs.size(); i++) {
		if (num_vecs != dfs[i]->get_num_vecs()) {
			BOOST_LOG_TRIVIAL(error)
				<< "The data frames have different numbers of vectors";
			return data_frame::ptr();
		}
	}

	data_frame::ptr df = data_frame::create();
	for (size_t vec_idx = 0; vec_idx < num_vecs; vec_idx++) {
		std::string vec_name = dfs[0]->named_vecs[vec_idx].first;
		const scalar_type &vec_type
			= dfs[0]->named_vecs[vec_idx].second->get_type();
		detail::vec_store::ptr vec = dfs[0]->named_vecs[vec_idx].second->deep_copy();

		// Here we collect the same column from all the data frame
		// except the first one.
		std::vector<detail::vec_store::const_ptr> vecs;
		for (size_t df_idx = 1; df_idx < dfs.size(); df_idx++) {
			if (vec_name != dfs[df_idx]->named_vecs[vec_idx].first
					|| vec_type != dfs[df_idx]->named_vecs[vec_idx].second->get_type()) {
				BOOST_LOG_TRIVIAL(error)
					<< "The data frames have different vectors";
				return data_frame::ptr();
			}
			vecs.push_back(dfs[df_idx]->named_vecs[vec_idx].second);
		}
		vec->append(vecs.begin(), vecs.end());
		df->add_vec(vec_name, vec);
	}
	return df;
}

std::vector<off_t> partition_vector(const detail::mem_vec_store &sorted_vec,
		int num_parts);

namespace
{

class vv_index_append;

/*
 * This generates a task on a partition of a data frame.
 */
class groupby_task_queue
{
	// The global offset of the partition on the data frame.
	off_t global_off;
	// The local offsets in the partition.
	std::vector<off_t> par_starts;
	// This indicates the current task.
	std::atomic<size_t> curr_idx;
	size_t first_task_idx;
public:
	typedef std::shared_ptr<groupby_task_queue> ptr;

	class task {
		// The start in the global vector.
		off_t global_start;
		// The start of the partition.
		off_t lstart;
		// The end of the partition.
		off_t lend;
		// The index of the partition in the vector.
		size_t idx;
	public:
		task(off_t global_start, off_t lstart, off_t lend, size_t idx) {
			this->global_start = global_start;
			this->lstart = lstart;
			this->lend = lend;
			this->idx = idx;
		}

		bool is_valid() const {
			return lstart >= 0 && lend > 0;
		}

		off_t get_global_start() const {
			return global_start;
		}

		off_t get_local_start() const {
			return lstart;
		}

		size_t get_length() const {
			return lend - lstart;
		}

		size_t get_idx() const {
			return idx;
		}
	};

	static ptr create(off_t global_off,
			detail::mem_vec_store::const_ptr sorted_vec,
			bool is_vec, std::shared_ptr<vv_index_append> append);

	groupby_task_queue(off_t global_off,
			detail::mem_vec_store::const_ptr sorted_vec,
			size_t num_parts, size_t first_task_idx) {
		this->global_off = global_off;
		par_starts = partition_vector(*sorted_vec, num_parts);
		// It's possible that two partitions end up having the same start location
		// because the vector is small or a partition has only one value.
		assert(std::is_sorted(par_starts.begin(), par_starts.end()));
		auto end_par_starts = std::unique(par_starts.begin(), par_starts.end());
		par_starts.resize(end_par_starts - par_starts.begin());
		this->first_task_idx = first_task_idx;
		curr_idx = 0;
	}

	size_t get_num_parts() const {
		// The last one indicates the end of the vector for groupby.
		return par_starts.size() - 1;
	}

	task get_task() {
		size_t curr = curr_idx.fetch_add(1);
		if (curr >= get_num_parts())
			return task(-1, -1, -1, -1);
		else
			return task(global_off + par_starts[curr], par_starts[curr],
					par_starts[curr + 1], curr + first_task_idx);
	}
};

/*
 * This helps to append data from multiple threads to a single vv store.
 * Each append is associated with a number. Data needs to be appended
 * according to the specified order.
 * The challenges here are:
 * * each append has an arbitrary size,
 * * the data may come in an arbitrary order,
 * * we need to maximize parallelization.
 * When the data is written to an EM vv store,
 * * the data has to be written to buffers with aligned memory.
 */
class vv_index_append
{
	class filled_chunk {
		local_vec_store::ptr data;
		std::atomic<size_t> filled_size;
	public:
		typedef std::shared_ptr<filled_chunk> ptr;

		filled_chunk(local_vec_store::ptr data) {
			this->data = data;
			filled_size = 0;
		}

		bool is_empty() const {
			return data == NULL;
		}

		bool fully_filled(bool last) const {
			if (data == NULL)
				return false;
			else if (last)
				return data->get_length() == filled_size;
			else
				return data->get_reserved_size() == filled_size;
		}

		local_vec_store::ptr get_vec() {
			return data;
		}

		void set_portion(local_vec_store::const_ptr buf, off_t off) {
			assert(data);
			data->set_portion(buf, off);
			filled_size += buf->get_length();
		}

	};
	typedef std::shared_ptr<std::vector<off_t> > off_vec_ptr;
	// The offsets of vectors in the vv store.
	// The key is the append index, the value is the offsets of
	// the vectors from the append in the global vv store.
	std::map<off_t, off_vec_ptr> vv_off_map;
	// When we write data to the global vector in the external memory,
	// we need to make sure data is stored in aligned memory and is written
	// to the aligned offset in a file.
	// We also need to make sure data is stored in the global order.
	// Therefore, we partition the data in the global vector into many
	// fixed-size chunks. Whenever a chunk is full, we write it to the
	// external-memory vector and reuse the memory for another chunk.
	std::deque<filled_chunk::ptr> global_chunks;
	// This is the global vector where we'll write data to.
	detail::vec_store::ptr global_vec;
	size_t chunk_size;

	std::atomic<size_t> tot_num_appends;
	std::atomic<size_t> num_appends;
	std::atomic<size_t> num_bytes_append;
	// The last append to the buffer.
	off_t last_append_idx;
	// The buffer contains the vectors from the threads.
	std::map<off_t, detail::vv_store::const_ptr> buf;
	spin_lock lock;

	const scalar_type &ele_type;

	void fill_chunks(const detail::vec_store &data,
			std::deque<filled_chunk::ptr> &chunks,
			std::deque<std::pair<off_t, size_t> > &off_sizes);
	std::vector<off_t> get_vv_offs();
public:
	typedef std::shared_ptr<vv_index_append> ptr;

	vv_index_append(detail::vec_store::ptr vec): ele_type(vec->get_type()) {
		this->global_vec = vec;
		this->chunk_size = matrix_conf.get_write_io_buf_size() / ele_type.get_size();
		last_append_idx = -1;
		this->tot_num_appends = 0;
		this->num_appends = 0;
		// The vector might already have some data.
		this->num_bytes_append = vec->get_length() * vec->get_type().get_size();
	}

	std::vector<safs::io_interface::ptr> create_ios() const {
		if (global_vec->is_in_mem())
			return std::vector<safs::io_interface::ptr>();
		else {
			detail::EM_object *obj
				= dynamic_cast<detail::EM_object *>(global_vec.get());
			assert(obj);
			return obj->create_ios();
		}
	}

	// This is called in the worker threads.
	void append(off_t idx, detail::vv_store::const_ptr vv,
			size_t max_num_bytes_taken);
	void append(off_t idx, detail::vv_store::const_ptr vv);

	// The methods below are called in the main thread.

	// This method returns the complete vv store.
	void flush() {
		append(-1, NULL, std::numeric_limits<size_t>::max());
	}
	detail::vv_store::ptr get_vv();

	void inc_tot_appends(size_t new_appends) {
		tot_num_appends += new_appends;
	}
	size_t get_tot_appends() const {
		return tot_num_appends.load();
	}

	bool has_data() {
		lock.lock();
		bool ret;
		if (num_appends < tot_num_appends)
			ret = true;
		else
			ret = !global_chunks.empty();
		lock.unlock();
		return ret;
	}
};

std::vector<off_t> vv_index_append::get_vv_offs()
{
	lock.lock();
	assert(num_appends == tot_num_appends);

	assert(buf.empty());
	assert(!vv_off_map.empty());
	auto it = vv_off_map.begin();
	std::vector<off_t> vv_offs = *it->second;
	for (it++; it != vv_off_map.end(); it++) {
		assert(vv_offs.back() == it->second->front());
		vv_offs.insert(vv_offs.end(), it->second->begin() + 1,
				it->second->end());
	}

	lock.unlock();
	return vv_offs;
}

/*
 * This method is only called to construct in-memory vv store.
 */
detail::vv_store::ptr vv_index_append::get_vv()
{
	flush();
	auto offs = get_vv_offs();
	for (size_t i = 0; i < global_chunks.size(); i++) {
		assert(global_chunks[i]->fully_filled(true));
		auto vec = global_chunks[i]->get_vec();
		assert(global_vec->get_length() + vec->get_length()
				<= global_vec->get_reserved_size());
		global_vec->append(*vec);
	}
	return detail::vv_store::create(offs, global_vec);
}

/*
 * Here we try to write data to the chunks in the specified locations.
 */
void vv_index_append::fill_chunks(const detail::vec_store &data,
		std::deque<filled_chunk::ptr> &chunks,
		std::deque<std::pair<off_t, size_t> > &off_sizes)
{
	size_t data_off = 0;
	size_t data_size = data.get_length();
	while (data_off < data.get_length()) {
		off_t off = off_sizes.front().first;
		size_t size = std::min(off_sizes.front().second, data_size);
		local_vec_store::const_ptr buf = data.get_portion(data_off, size);
		chunks.front()->set_portion(buf, off);
		data_off += size;
		data_size -= size;
		off_sizes.front().first += size;
		off_sizes.front().second -= size;
		if (off_sizes.front().second == 0) {
			chunks.pop_front();
			off_sizes.pop_front();
		}
	}

}

void vv_index_append::append(off_t idx, detail::vv_store::const_ptr vv)
{
	// We take out at most some number of bytes from the buffer so that writing
	// data to chunks can be parallelized among multiple threads.
	size_t max_num_bytes_taken;
	if (global_vec->is_in_mem())
		max_num_bytes_taken = std::numeric_limits<size_t>::max();
	else
		max_num_bytes_taken = std::max(vv->get_num_bytes() * 2,
				chunk_size * ele_type.get_size() * 2);
	append(idx, vv, max_num_bytes_taken);
}

void vv_index_append::append(off_t idx, detail::vv_store::const_ptr vv,
		size_t max_num_bytes_taken)
{
	off_t first_idx = -1;
	std::vector<detail::vv_store::const_ptr> contig;
	std::vector<off_vec_ptr> vv_off_vec;
	// The data chunks where we'll write data to.
	std::deque<filled_chunk::ptr> chunks;
	// The locations and sizes we'll write data to the chunks above.
	std::deque<std::pair<off_t, size_t> > off_sizes;
	// The last offset in bytes in the global vec..
	off_t last_off_bytes = 0;
	lock.lock();
	// Add the vv store.
	if (vv)
		buf.insert(std::pair<off_t, detail::vv_store::const_ptr>(idx, vv));

	// Get contiguous vv.
	auto it = buf.begin();
	size_t num_bytes = 0;
	while (it != buf.end() && last_append_idx + 1 == it->first
			&& num_bytes < max_num_bytes_taken) {
		if (first_idx < 0)
			first_idx = it->first;
		contig.push_back(it->second);
		num_bytes += contig.back()->get_num_bytes();
		last_append_idx++;
		it++;
	}
	if (it != buf.begin())
		buf.erase(buf.begin(), it);
	size_t max_num_chunks_to_write
		= num_bytes / (chunk_size * ele_type.get_size()) + 2;

	// If we have contiguous vv, we insert them to the global data structure.
	// It takes two steps:
	// * determine the location in the global vector where vvs are written to.
	// * write data to the determined location.
	// We can do this because we have reserved space in advance, so we don't
	// need to reallocate space.

	// Here is to determine the location.
	size_t num_new_vvs = 0;
	for (size_t i = 0; i < contig.size(); i++) {
		num_new_vvs += contig[i]->get_num_vecs();
		auto p = off_vec_ptr(new std::vector<off_t>());
		vv_off_vec.push_back(p);
		vv_off_map.insert(std::pair<off_t, off_vec_ptr>(first_idx + i, p));
	}
	last_off_bytes = num_bytes_append.load();
	num_bytes_append += num_bytes;

	size_t fill_size = num_bytes / ele_type.get_size();
	// We assume the preallocated space is large enough to store all
	// appended data.
	assert(global_vec->get_length() + fill_size
			<= global_vec->get_reserved_size());
	off_t global_off = global_vec->get_length();
	size_t global_size = fill_size;
	// For in-memory vector, we resize the vector and record the location
	// where we should write data to.
	if (fill_size > 0 && global_vec->is_in_mem()) {
		global_vec->resize(global_vec->get_length() + fill_size);
		fill_size = 0;
	}
	// For external-memory vector, we first write data to aligned memory
	// buffers and record the locations where the aligned buffers should be
	// written to.
	else if (fill_size > 0 && !global_vec->is_in_mem()) {
		local_vec_store::ptr last_vec;
		off_t portion_global_start;
		if (global_chunks.empty())
			portion_global_start = global_off;
		else {
			last_vec = global_chunks.back()->get_vec();
			// The next local vector starts at the ending point of this local vector
			portion_global_start = last_vec->get_global_start() + chunk_size;
		}

		if (last_vec && last_vec->get_length() < last_vec->get_reserved_size()) {
			chunks.push_back(global_chunks.back());
			size_t num_eles = std::min(fill_size,
					last_vec->get_reserved_size() - last_vec->get_length());
			off_sizes.emplace_back(last_vec->get_length(), num_eles);
			last_vec->resize(last_vec->get_length() + num_eles);
			fill_size -= num_eles;
		}
		while (fill_size > 0) {
			size_t num_eles = std::min(fill_size, chunk_size);
			// We store data in the local vectors.
			local_vec_store::ptr buf(new local_buf_vec_store(
						portion_global_start, chunk_size, ele_type, -1, false));
			portion_global_start += chunk_size;
			buf->resize(num_eles);
			global_chunks.emplace_back(new filled_chunk(buf));
			fill_size -= num_eles;
			chunks.push_back(global_chunks.back());
			off_sizes.emplace_back(0, num_eles);
		}
	}
	assert(fill_size == 0);

	lock.unlock();

	// Here is to write data to the right location.
	for (size_t i = 0; i < contig.size(); i++) {
		detail::vv_store::const_ptr vv = contig[i];
		const detail::vec_store &data = vv->get_data();
		if (global_vec->is_in_mem()) {
			local_vec_store::const_ptr buf = data.get_portion(0,
					data.get_length());
			global_vec->set_portion(buf, global_off);
			global_off += data.get_length();
			global_size -= data.get_length();
		}
		else
			fill_chunks(data, chunks, off_sizes);

		// Fill the offsets for each vector in the vv store.
		std::vector<off_t> &vv_offs = *vv_off_vec[i];
		vv_offs.push_back(last_off_bytes);
		for (size_t j = 0; j < vv->get_num_vecs(); j++)
			vv_offs.push_back(vv_offs.back() + vv->get_num_bytes(j));
		last_off_bytes += vv->get_num_bytes();
	}
	assert(chunks.empty());
	assert(off_sizes.empty());

	std::vector<local_vec_store::const_ptr> filled_chunks;
	lock.lock();
	num_appends += contig.size();
	// Here is to count the number of chunks that have been completely filled,
	// so we can write them to the global vector.
	size_t to_global_size = 0;
	while (!global_chunks.empty() && global_chunks.front()->fully_filled(false)
			&& filled_chunks.size() <= max_num_chunks_to_write) {
		auto vec = global_chunks.front()->get_vec();
		to_global_size += vec->get_length();
		assert(global_vec->get_length() + to_global_size
				>= vec->get_global_start() + vec->get_length());
		filled_chunks.push_back(vec);
		global_chunks.pop_front();
	}
	if (to_global_size > 0)
		global_vec->resize(global_vec->get_length() + to_global_size);
	size_t global_len = global_vec->get_length();
	lock.unlock();

	if (!filled_chunks.empty()) {
		assert(!global_vec->is_in_mem());
		detail::EM_vec_store::ptr em_vec
			= std::static_pointer_cast<detail::EM_vec_store>(global_vec);
		for (size_t i = 0; i < filled_chunks.size(); i++) {
			auto vec = filled_chunks[i];
			assert(vec->get_global_start() + vec->get_length() <= global_len);
			em_vec->write_portion_async(vec, vec->get_global_start());
		}
	}
}

class portion_groupby_complete;

class local_groupby_task: public thread_task
{
	groupby_task_queue::ptr q;
	// If the data frame is stored in EM, sorted_col is a local vector.
	// We have this field to avoid reading the data from disks again.
	detail::mem_vec_store::const_ptr sorted_col;
	size_t sorted_col_idx;
	data_frame::const_ptr df;
	vv_index_append::ptr append;
	const gr_apply_operate<sub_data_frame> &op;
	bool last_task;
public:
	local_groupby_task(groupby_task_queue::ptr q, size_t sorted_col_idx,
			detail::mem_vec_store::const_ptr sorted_col, data_frame::const_ptr df,
			const gr_apply_operate<sub_data_frame> &_op,
			vv_index_append::ptr append, bool last_task): op(_op) {
		this->q = q;
		this->sorted_col_idx = sorted_col_idx;
		this->sorted_col = sorted_col;
		this->df = df;
		this->append = append;
		this->last_task = last_task;
	}

	std::shared_ptr<portion_groupby_complete> get_part(
			const groupby_task_queue::task &t, local_buf_vec_store::ptr row_buf,
			bool last_task);

	void run();
};

class portion_groupby_complete: public detail::portion_compute
{
	// The local buf to store the temporary result.
	// This buffer is shared by many portion_groupby_complete tasks.
	// However, all of these tasks run in the same thread, it's fine
	// to share the same buffer. There is no concurrent access to it.
	local_buf_vec_store::ptr row_buf;
	// A partition in the data frame for computation in this task.
	sub_data_frame df;
	size_t sorted_col_idx;
	int num_EM;
	const gr_apply_operate<sub_data_frame> &op;
	vv_index_append::ptr append;
	groupby_task_queue::task t;
	bool last_task;
public:
	typedef std::shared_ptr<portion_groupby_complete> ptr;

	portion_groupby_complete(vv_index_append::ptr append,
			local_buf_vec_store::ptr row_buf,
			const gr_apply_operate<sub_data_frame> &_op,
			const groupby_task_queue::task &_t,
			size_t sorted_col_idx, bool last_task): op(_op), t(_t) {
		this->row_buf = row_buf;
		this->append = append;
		this->num_EM = 0;
		this->sorted_col_idx = sorted_col_idx;
		this->last_task = last_task;
	}

	bool is_complete() const {
		return num_EM == 0;
	}

	void run_complete();

	virtual void run(char *buf, size_t size) {
		assert(num_EM > 0);
		num_EM--;
		if (num_EM == 0)
			run_complete();
	}

	void set_data(const sub_data_frame &df, int num_EM) {
		this->df = df;;
		this->num_EM = num_EM;
	}
};

static void expose_portion(sub_data_frame &sub_df, off_t loc, size_t length)
{
	for (size_t i = 0; i < sub_df.size(); i++)
		const_cast<local_vec_store *>(sub_df[i].get())->expose_sub_vec(loc,
				length);
}

void portion_groupby_complete::run_complete()
{
	// The local vector from EM vector might have been resized before.
	// so we should get a new local vector whose local offset is 0.
	sub_data_frame sub_df = df;
	for (size_t i = 0; i < sub_df.size(); i++) {
		if (sub_df[i]->get_local_start() > 0)
			sub_df[i] = sub_df[i]->get_portion(0, sub_df[i]->get_length());
	}

	auto sorted_col = sub_df[sorted_col_idx];
	detail::mem_vv_store::ptr ret = detail::mem_vv_store::create(
			op.get_output_type());
	agg_operate::const_ptr find_next
		= sorted_col->get_type().get_agg_ops().get_find_next();
	size_t loc = 0;
	size_t col_len = sorted_col->get_length();
	// We can't search a vv store.
	assert(!detail::vv_store::is_vector_vector(*sorted_col));
	const char *start = sorted_col->get_raw_arr();
	size_t entry_size = sorted_col->get_type().get_size();
	while (loc < col_len) {
		size_t curr_length = col_len - loc;
		const char *curr_ptr = start + entry_size * loc;
		size_t rel_end;
		find_next->runAgg(curr_length, curr_ptr, &rel_end);
		// This expose a portion of the data frame.
		expose_portion(sub_df, loc, rel_end);
		// The first argument is the key and the second one is the value
		// (a data frame)
		op.run(curr_ptr, sub_df, *row_buf);
		if (row_buf->get_length() > 0)
			ret->append(*row_buf);
		loc += rel_end;
	}
	append->append(t.get_idx(), ret);
	if (last_task)
		append->flush();
}

portion_groupby_complete::ptr local_groupby_task::get_part(
		const groupby_task_queue::task &t, local_buf_vec_store::ptr row_buf,
		bool last_task)
{
	sub_data_frame sub_df(df->get_num_vecs());

	// Read the part from the data frame.
	size_t num_EM = 0;
	portion_groupby_complete::ptr compute(new portion_groupby_complete(
				append, row_buf, op, t, sorted_col_idx, last_task));

	for (size_t i = 0; i < df->get_num_vecs(); i++) {
		if (i == sorted_col_idx)
			sub_df[i] = sorted_col->get_portion(t.get_local_start(),
					t.get_length());
		else {
			detail::vec_store::const_ptr vec = df->get_vec(i);
			local_vec_store::const_ptr col;
			if (vec->is_in_mem())
				col = vec->get_portion(t.get_global_start(), t.get_length());
			// TODO we need to unify the async get_portion interface of vec store.
			else if (detail::vv_store::is_vector_vector(*vec)) {
				detail::EM_vv_store::const_ptr EM_vv
					= std::dynamic_pointer_cast<const detail::EM_vv_store>(vec);
				assert(EM_vv);
				col = EM_vv->get_portion_async(t.get_global_start(),
						t.get_length(), compute);
				num_EM++;
				assert(detail::vv_store::is_vector_vector(*col));
			}
			else {
				detail::EM_vec_store::const_ptr EM_vec
					= std::dynamic_pointer_cast<const detail::EM_vec_store>(vec);
				assert(EM_vec);
				col = EM_vec->get_portion_async(t.get_global_start(),
						t.get_length(), compute);
				num_EM++;
			}
			sub_df[i] = col;
		}
	}

	compute->set_data(sub_df, num_EM);
	return compute;
}

void local_groupby_task::run()
{
	// Prepare I/O.
	size_t max_pending_ios = 1;
	std::vector<safs::io_interface::ptr> ios;
	for (size_t i = 0; i < df->get_num_vecs(); i++) {
		detail::vec_store::const_ptr vec = df->get_vec(i);
		const detail::EM_object *obj
			= dynamic_cast<const detail::EM_object *>(vec.get());
		if (obj) {
			std::vector<safs::io_interface::ptr> tmp = obj->create_ios();
			ios.insert(ios.end(), tmp.begin(), tmp.end());
		}
	}
	auto tmp = append->create_ios();
	ios.insert(ios.end(), tmp.begin(), tmp.end());
	safs::io_select::ptr select;
	// If we need to access I/O.
	if (!ios.empty())
		select = safs::create_io_select(ios);

	size_t out_size;
	// If the user can predict the number of output elements, we can create
	// a buffer of the expected size.
	if (op.get_num_out_eles() > 0)
		out_size = op.get_num_out_eles();
	else
		// If the user can't, we create a small buffer.
		out_size = 16;
	local_buf_vec_store::ptr row_buf(new local_buf_vec_store(0,
				out_size, op.get_output_type(), -1));
	while (true) {
		// Get the a partition from the queue.
		auto t = q->get_task();
		if (!t.is_valid())
			break;

		portion_groupby_complete::ptr compute = get_part(t, row_buf, last_task);
		if (compute->is_complete())
			compute->run_complete();
		if (select)
			safs::wait4ios(select, max_pending_ios);
	}

	// If this is the last task, we need to wait for all I/O to complete.
	if (last_task) {
		size_t num_pending = 0;
		if (select)
			num_pending = safs::wait4ios(select, 0);
		assert(num_pending == 0);

		for (size_t i = 0; i < ios.size(); i++) {
			detail::portion_callback &cb = static_cast<detail::portion_callback &>(
					ios[i]->get_callback());
			assert(!cb.has_callback());
		}
	}
}

groupby_task_queue::ptr groupby_task_queue::create(off_t global_off,
		detail::mem_vec_store::const_ptr sorted_col, bool is_vec,
		vv_index_append::ptr append)
{
	// If all fields in the data frame are vectors, we use a larger
	// partition size. Otherwise, we use a small partition size.
	size_t part_size = is_vec ? 1024 * 1024 : 16 * 1024;
	return groupby_task_queue::ptr(new groupby_task_queue(global_off,
				sorted_col, std::max(sorted_col->get_length() / part_size, 1UL),
				append->get_tot_appends()));
}

}

static bool is_all_vec(const data_frame &df)
{
	for (size_t i = 0; i < df.get_num_vecs(); i++) {
		auto vec = df.get_vec(i);
		if (vec->get_num_bytes()
				!= vec->get_length() * vec->get_type().get_size())
			return false;
	}
	return true;
}

static void parallel_groupby_async(data_frame::const_ptr sorted_df,
		off_t sorted_col_idx, const gr_apply_operate<sub_data_frame> &op,
		vv_index_append::ptr append)
{
	detail::mem_vec_store::const_ptr sorted_col
		= std::dynamic_pointer_cast<const detail::mem_vec_store>(
				sorted_df->get_vec(sorted_col_idx));
	groupby_task_queue::ptr q = groupby_task_queue::create(0,
			sorted_col, is_all_vec(*sorted_df), append);
	append->inc_tot_appends(q->get_num_parts());

	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	int num_threads = mem_threads->get_num_threads();
	for (int i = 0; i < num_threads; i++) {
		// It's difficult to localize computation.
		// TODO can we localize computation?
		int node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id, new local_groupby_task(q,
					sorted_col_idx, sorted_col, sorted_df, op, append, true));
	}
}

static vector_vector::ptr in_mem_groupby(
		data_frame::const_ptr sorted_df, const std::string &col_name,
		const gr_apply_operate<sub_data_frame> &op,
		detail::mem_vec_store::ptr out_vec)
{
	off_t sorted_col_idx = -1;
	for (size_t i = 0; i < sorted_df->get_num_vecs(); i++) {
		if (sorted_df->get_vec_name(i) == col_name)
			sorted_col_idx = i;
	}
	assert(sorted_col_idx >= 0);

	vv_index_append::ptr append(new vv_index_append(out_vec));

	parallel_groupby_async(sorted_df, sorted_col_idx, op, append);
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	mem_threads->wait4complete();
	auto vv = append->get_vv();
	assert(vv);
	return vector_vector::create(vv);
}

namespace
{

class EM_df_groupby_dispatcher: public detail::task_dispatcher
{
	off_t ele_idx;
	size_t portion_size;
	off_t sorted_col_idx;
	detail::vec_store::const_ptr groupby_col;
	data_frame::const_ptr df;
	detail::EM_vec_store::ptr out_vec;
	const gr_apply_operate<sub_data_frame> &op;
	vv_index_append::ptr append;
public:
	typedef std::shared_ptr<EM_df_groupby_dispatcher> ptr;

	EM_df_groupby_dispatcher(data_frame::const_ptr df, const std::string &col_name,
			detail::EM_vec_store::ptr out_vec, size_t portion_size,
			const gr_apply_operate<sub_data_frame> &_op): op(_op) {
		ele_idx = 0;
		this->portion_size = portion_size;

		sorted_col_idx = -1;
		for (size_t i = 0; i < df->get_num_vecs(); i++) {
			if (df->get_vec_name(i) == col_name)
				sorted_col_idx = i;
		}
		assert(sorted_col_idx >= 0);

		this->sorted_col_idx = sorted_col_idx;
		this->groupby_col = df->get_vec(col_name);
		this->df = df;
		this->out_vec = out_vec;
		append = vv_index_append::ptr(new vv_index_append(out_vec));
	}

	virtual bool issue_task();

	detail::vv_store::ptr get_vv() {
		return append->get_vv();
	}
};

bool EM_df_groupby_dispatcher::issue_task()
{
	if ((size_t) ele_idx >= groupby_col->get_length())
		return false;

	// Indicate whether this portion is the last portion in the vector.
	bool last_part;
	// The number of elements we need to read from the sorted vector.
	size_t read_len;
	// The number of elements we need to expose on the vectors.
	size_t real_len;
	if (portion_size >= groupby_col->get_length() - ele_idx) {
		last_part = true;
		read_len = groupby_col->get_length() - ele_idx;
	}
	else {
		last_part = false;
		read_len = portion_size;
	}
	real_len = read_len;

	local_vec_store::const_ptr vec = groupby_col->get_portion(ele_idx, read_len);
	// If this isn't the last portion, we should expose the part of
	// the sorted vector that all elements with the same values in the vector
	// are guaranteed in the part.
	if (!last_part) {
		agg_operate::const_ptr find_prev
			= vec->get_type().get_agg_ops().get_find_prev();
		size_t off;
		size_t entry_size = vec->get_type().get_size();
		find_prev->runAgg(read_len, vec->get_raw_arr() + read_len * entry_size,
				&off);
		// All keys are the same in this portion of the vector. Then we don't
		// know if we have got all values for the key.
		// If we know we can ignore the key, we can ignore the entire portion.
		if (off == read_len && op.ignore_key(vec->get_raw_arr()))
			off = 0;
		// Otherwise, we have to make sure we have seen all values for the keys
		// we are going to run UDFs.
		else
			assert(off < read_len);
		real_len = read_len - off;
		// The local buffer may already be a sub vector.
		const_cast<local_vec_store &>(*vec).expose_sub_vec(vec->get_local_start(),
				real_len);
	}

	groupby_task_queue::ptr groupby_q = groupby_task_queue::create(ele_idx,
			vec, is_all_vec(*df), append);
	append->inc_tot_appends(groupby_q->get_num_parts());

	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	int num_threads = mem_threads->get_num_threads();
	for (int i = 0; i < num_threads; i++) {
		// It's difficult to localize computation.
		// TODO can we localize computation?
		int node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id, new local_groupby_task(
					groupby_q, sorted_col_idx, vec, df, op, append, last_part));
	}

	ele_idx += real_len;
	return true;
}

}

static vector_vector::ptr EM_groupby(
		data_frame::const_ptr sorted_df, const std::string &col_name,
		const gr_apply_operate<sub_data_frame> &op,
		detail::EM_vec_store::ptr out_vec)
{
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	int num_threads = mem_threads->get_num_threads();

	size_t portion_size = num_threads * 1024 * 1024 * 5;
	// If one of the columns in the data frame is a vector vector,
	// we should use a smaller portion size.
	for (size_t i = 0; i < sorted_df->get_num_vecs(); i++)
		if (sorted_df->get_vec_ref(i).get_entry_size() == 0) {
			portion_size = num_threads * 16 * 1024 * 5;
			break;
		}

	EM_df_groupby_dispatcher::ptr groupby_dispatcher(
			new EM_df_groupby_dispatcher(sorted_df, col_name, out_vec,
				portion_size, op));
	detail::io_worker_task worker(groupby_dispatcher, 1);
	for (size_t i = 0; i < sorted_df->get_num_vecs(); i++) {
		detail::vec_store::const_ptr col = sorted_df->get_vec(i);
		if (col->is_in_mem())
			continue;

		const detail::EM_object *obj
			= dynamic_cast<const detail::EM_object *>(col.get());
		worker.register_EM_obj(const_cast<detail::EM_object *>(obj));
	}
	worker.register_EM_obj(out_vec.get());
	worker.run();
	mem_threads->wait4complete();
	return vector_vector::create(groupby_dispatcher->get_vv());
}

vector_vector::ptr data_frame::groupby(const std::string &col_name,
		const gr_apply_operate<sub_data_frame> &op,
		detail::vec_store::ptr store) const
{
	const scalar_type &type = op.get_output_type();
	if (store == NULL) {
		// We use the storage size of the data frame to approximate the storage size
		// for the groupby result.
		size_t num_bytes = 0;
		for (size_t i = 0; i < get_num_vecs(); i++)
			num_bytes += get_vec(i)->get_num_bytes();
		store = detail::vec_store::create(0, type, -1, is_in_mem());
		store->reserve(num_bytes / type.get_size());
	}

	data_frame::const_ptr sorted_df = sort(col_name);
	if (is_in_mem()) {
		assert(store->is_in_mem());
		return in_mem_groupby(sorted_df, col_name, op,
				std::dynamic_pointer_cast<detail::mem_vec_store>(store));
	}
	else {
		assert(!store->is_in_mem());
		return EM_groupby(sorted_df, col_name, op,
				std::dynamic_pointer_cast<detail::EM_vec_store>(store));
	}
}

}
