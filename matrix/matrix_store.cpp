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

#include "matrix_store.h"
#include "local_matrix_store.h"
#include "mem_matrix_store.h"
#include "EM_dense_matrix.h"
#include "EM_object.h"
#include "matrix_config.h"
#include "local_mem_buffer.h"

namespace fm
{

namespace detail
{

std::atomic<size_t> matrix_store::mat_counter;

matrix_store::ptr matrix_store::create(size_t nrow, size_t ncol,
		matrix_layout_t layout, const scalar_type &type, int num_nodes,
		bool in_mem, safs::safs_file_group::ptr group)
{
	if (in_mem)
		return mem_matrix_store::create(nrow, ncol, layout, type, num_nodes);
	else
		return EM_matrix_store::create(nrow, ncol, layout, type, group);
}

matrix_store::matrix_store(size_t nrow, size_t ncol, bool in_mem,
		const scalar_type &_type): type(_type)
{
	this->nrow = nrow;
	this->ncol = ncol;
	this->in_mem = in_mem;
	this->entry_size = type.get_size();
	this->cache_portion = true;
}

size_t matrix_store::get_num_portions() const
{
	std::pair<size_t, size_t> chunk_size = get_portion_size();
	if (is_wide())
		return ceil(((double) get_num_cols()) / chunk_size.second);
	else
		return ceil(((double) get_num_rows()) / chunk_size.first);
}

local_matrix_store::ptr matrix_store::get_portion(size_t id)
{
	size_t start_row;
	size_t start_col;
	size_t num_rows;
	size_t num_cols;
	std::pair<size_t, size_t> chunk_size = get_portion_size();
	if (is_wide()) {
		start_row = 0;
		start_col = chunk_size.second * id;
		num_rows = get_num_rows();
		num_cols = std::min(chunk_size.second, get_num_cols() - start_col);
	}
	else {
		start_row = chunk_size.first * id;
		start_col = 0;
		num_rows = std::min(chunk_size.first, get_num_rows() - start_row);
		num_cols = get_num_cols();
	}
	return get_portion(start_row, start_col, num_rows, num_cols);
}

local_matrix_store::const_ptr matrix_store::get_portion(size_t id) const
{
	size_t start_row;
	size_t start_col;
	size_t num_rows;
	size_t num_cols;
	std::pair<size_t, size_t> chunk_size = get_portion_size();
	if (is_wide()) {
		start_row = 0;
		start_col = chunk_size.second * id;
		num_rows = get_num_rows();
		num_cols = std::min(chunk_size.second, get_num_cols() - start_col);
	}
	else {
		start_row = chunk_size.first * id;
		start_col = 0;
		num_rows = std::min(chunk_size.first, get_num_rows() - start_row);
		num_cols = get_num_cols();
	}
	return get_portion(start_row, start_col, num_rows, num_cols);
}

namespace
{

class reset_op: public set_operate
{
	const scalar_type &type;
	size_t entry_size;
public:
	reset_op(const scalar_type &_type): type(_type) {
		this->entry_size = type.get_size();
	}

	virtual void set(void *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		memset(arr, 0, num_eles * entry_size);
	}
	virtual const scalar_type &get_type() const {
		return type;
	}
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

class set_task: public thread_task
{
	detail::local_matrix_store::ptr local_store;
	const set_operate &op;
public:
	set_task(detail::local_matrix_store::ptr local_store,
			const set_operate &_op): op(_op) {
		this->local_store = local_store;
	}

	void run() {
		local_store->set_data(op);
	}
};

/*
 * These two functions define the length and portion size for 1D partitioning
 * on a matrix.
 */

static inline size_t get_tot_len(const matrix_store &mat)
{
	    return mat.is_wide() ? mat.get_num_cols() : mat.get_num_rows();
}

static inline size_t get_portion_size(const matrix_store &mat)
{
	    return mat.is_wide() ? mat.get_portion_size().second : mat.get_portion_size().first;
}

class EM_mat_setdata_dispatcher: public EM_portion_dispatcher
{
	const set_operate &op;
	matrix_store &to_mat;
public:
	EM_mat_setdata_dispatcher(matrix_store &store, const set_operate &_op);

	virtual void create_task(off_t global_start, size_t length);
};

EM_mat_setdata_dispatcher::EM_mat_setdata_dispatcher(matrix_store &store,
		const set_operate &_op): EM_portion_dispatcher(get_tot_len(store),
			fm::detail::get_portion_size(store)), op(_op), to_mat(store)
{
}

void EM_mat_setdata_dispatcher::create_task(off_t global_start,
		size_t length)
{
	size_t global_start_row, global_start_col;
	size_t num_rows, num_cols;
	if (to_mat.is_wide()) {
		global_start_row = 0;
		global_start_col = global_start;
		num_rows = to_mat.get_num_rows();
		num_cols = length;
	}
	else {
		global_start_row = global_start;
		global_start_col = 0;
		num_rows = length;
		num_cols = to_mat.get_num_cols();
	}
	local_matrix_store::ptr buf;
	if (to_mat.store_layout() == matrix_layout_t::L_COL)
		buf = local_matrix_store::ptr(new local_buf_col_matrix_store(
					global_start_row, global_start_col, num_rows, num_cols,
					to_mat.get_type(), -1));
	else
		buf = local_matrix_store::ptr(new local_buf_row_matrix_store(
					global_start_row, global_start_col, num_rows, num_cols,
					to_mat.get_type(), -1));
	buf->set_data(op);
	to_mat.write_portion_async(buf, global_start_row, global_start_col);
}

}

void matrix_store::reset_data()
{
	set_data(reset_op(get_type()));
}

void matrix_store::set_data(const set_operate &op)
{
	size_t num_chunks = get_num_portions();
	if (is_in_mem() && num_chunks == 1) {
		local_matrix_store::ptr buf;
		if (store_layout() == matrix_layout_t::L_ROW)
			buf = local_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
						get_num_rows(), get_num_cols(), get_type(), -1));
		else
			buf = local_matrix_store::ptr(new local_buf_col_matrix_store(0, 0,
						get_num_rows(), get_num_cols(), get_type(), -1));
		buf->set_data(op);
		write_portion_async(buf, 0, 0);
		// After computation, some matrices buffer local portions in the thread,
		// we should try to clean these local portions. These local portions
		// may contain pointers to some matrices that don't exist any more.
		// We also need to clean them to reduce memory consumption.
		// We might want to keep the memory buffer for I/O on dense matrices.
		if (matrix_conf.is_keep_mem_buf())
			detail::local_mem_buffer::clear_bufs(
					detail::local_mem_buffer::MAT_PORTION);
		else
			detail::local_mem_buffer::clear_bufs();
	}
	else if (is_in_mem()) {
		detail::mem_thread_pool::ptr mem_threads
			= detail::mem_thread_pool::get_global_mem_threads();
		for (size_t i = 0; i < num_chunks; i++) {
			detail::local_matrix_store::ptr local_store = get_portion(i);

			int node_id = local_store->get_node_id();
			// If the local matrix portion is not assigned to any node, 
			// assign the tasks in round robin fashion.
			if (node_id < 0)
				node_id = i % mem_threads->get_num_nodes();
			mem_threads->process_task(node_id, new set_task(local_store, op));
		}
		mem_threads->wait4complete();
	}
	else {
		mem_thread_pool::ptr threads = mem_thread_pool::get_global_mem_threads();
		EM_mat_setdata_dispatcher::ptr dispatcher(
				new EM_mat_setdata_dispatcher(*this, op));
		EM_matrix_store *em_this = dynamic_cast<EM_matrix_store *>(this);
		assert(em_this);
		em_this->start_stream();
		for (size_t i = 0; i < threads->get_num_threads(); i++) {
			io_worker_task *task = new io_worker_task(dispatcher);
			const EM_object *obj = dynamic_cast<const EM_object *>(this);
			task->register_EM_obj(const_cast<EM_object *>(obj));
			threads->process_task(i % threads->get_num_nodes(), task);
		}
		threads->wait4complete();
		em_this->end_stream();
	}
}

matrix_stream::ptr matrix_stream::create(matrix_store::ptr store)
{
	if (store->is_in_mem()) {
		mem_matrix_store::ptr mem_store
			= std::dynamic_pointer_cast<mem_matrix_store>(store);
		if (mem_store == NULL) {
			BOOST_LOG_TRIVIAL(error)
				<< "The in-mem matrix store isn't writable";
			return matrix_stream::ptr();
		}
		else
			return mem_matrix_stream::create(mem_store);
	}
	else {
		EM_matrix_store::ptr em_store
			= std::dynamic_pointer_cast<EM_matrix_store>(store);
		if (em_store == NULL) {
			BOOST_LOG_TRIVIAL(error)
				<< "The ext-mem matrix store isn't writable";
			return matrix_stream::ptr();
		}
		else
			return EM_matrix_stream::create(em_store);
	}
}

matrix_store::const_ptr matrix_store::get_cols(
			const std::vector<off_t> &idxs) const
{
	matrix_store::const_ptr tm = transpose();
	matrix_store::const_ptr rows = tm->get_rows(idxs);
	if (rows == NULL)
		return matrix_store::const_ptr();
	else
		return rows->transpose();
}

matrix_store::const_ptr matrix_store::get_cols(off_t start, off_t end) const
{
	if (start < 0 || end < 0 || end - start < 0) {
		BOOST_LOG_TRIVIAL(error) << "invalid range for selecting columns";
		return matrix_store::const_ptr();
	}

	std::vector<off_t> idxs(end - start);
	for (size_t i = 0; i < idxs.size(); i++)
		idxs[i] = start + i;
	return get_cols(idxs);
}

matrix_store::const_ptr matrix_store::get_rows(off_t start, off_t end) const
{
	if (start < 0 || end < 0 || end - start < 0) {
		BOOST_LOG_TRIVIAL(error) << "invalid range for selecting rows";
		return matrix_store::const_ptr();
	}

	std::vector<off_t> idxs(end - start);
	for (size_t i = 0; i < idxs.size(); i++)
		idxs[i] = start + i;
	return get_rows(idxs);
}

bool matrix_store::share_data(const matrix_store &store) const
{
	// By default, we can use data id to determine if two matrices have
	// the same data.
	return get_data_id() == store.get_data_id()
		&& get_data_id() != INVALID_MAT_ID;
}

matrix_append::matrix_append(matrix_store::ptr store)
{
	this->res = store;
	q.resize(1000);
	last_append = -1;
	written_eles = 0;
	empty_portion = local_matrix_store::const_ptr(new local_buf_row_matrix_store(
				0, 0, 0, 0, store->get_type(), -1, false));
}

void matrix_append::write_async(local_matrix_store::const_ptr portion,
		off_t seq_id)
{
	if (seq_id <= last_append) {
		BOOST_LOG_TRIVIAL(error) << "Append a repeated portion";
		return;
	}

	if (portion == NULL)
		portion = empty_portion;
	std::vector<local_matrix_store::const_ptr> data;
	lock.lock();
	// Add the new portion to the queue. If the queue is too small,
	// we should resize the queue first.
	off_t loc = seq_id - last_append - 1;
	assert(loc >= 0);
	if ((size_t) loc >= q.size())
		q.resize(q.size() * 2);
	q[loc] = portion;

	off_t start_loc = -1;
	if (q.front())
		start_loc = written_eles;
	// Get the portions from the queue.
	while (q.front()) {
		auto mat = q.front();
		// If the portion isn't empty.
		if (mat->get_num_rows() > 0 && mat->get_num_cols() > 0)
			data.push_back(mat);
		q.pop_front();
		q.push_back(local_matrix_store::const_ptr());
		last_append++;
		written_eles += mat->get_num_rows() * mat->get_num_cols();
	}
	lock.unlock();

	for (size_t i = 0; i < data.size(); i++) {
		assert(start_loc >= 0);
		// TODO this works if the result matrix is stored in memory.
		if (res->is_wide()) {
			off_t start_row = 0;
			off_t start_col = start_loc / res->get_num_rows();
			res->write_portion_async(data[i], start_row, start_col);
		}
		else {
			off_t start_row = start_loc / res->get_num_cols();
			off_t start_col = 0;
			res->write_portion_async(data[i], start_row, start_col);
		}
		start_loc += data[i]->get_num_rows() * data[i]->get_num_cols();
	}
}

matrix_append::~matrix_append()
{
	for (size_t i = 0; i < q.size(); i++)
		assert(q[i] == NULL);
}

void matrix_append::flush()
{
	for (size_t i = 0; i < q.size(); i++)
		assert(q[i] == NULL);
}

}

}
