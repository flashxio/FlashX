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
#include "rand_gen.h"

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
		detail::mem_thread_pool::disable_thread_pool();
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
		detail::mem_thread_pool::enable_thread_pool();
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
		for (size_t i = 0; i < threads->get_num_threads(); i++) {
			io_worker_task *task = new io_worker_task(dispatcher);
			const EM_object *obj = dynamic_cast<const EM_object *>(this);
			task->register_EM_obj(const_cast<EM_object *>(obj));
			threads->process_task(i % threads->get_num_nodes(), task);
		}
		threads->wait4complete();
	}
}

namespace
{

/*
 * This class set elements in a container randomly.
 * set_operate can't change its own state and has to be thread-safe when
 * running on multiple threads. However, random generators aren't
 * thread-safe, so we have to create a random generator for each thread.
 */
class rand_init: public set_operate
{
public:
	enum rand_dist_type {
		NORM,
		UNIF,
		MAX_NUM,
	};
private:
	class rand_gen_wrapper {
		rand_gen::ptr gen;
	public:
		rand_gen_wrapper(rand_gen::ptr gen) {
			this->gen = gen;
		}

		rand_gen &get_gen() {
			return *gen;
		}
	};

	pthread_key_t gen_key;
	const scalar_type &type;
	const scalar_variable &var1;
	const scalar_variable &var2;
	rand_dist_type rand_dist;

	rand_gen &get_rand_gen() const {
		void *addr = pthread_getspecific(gen_key);
		if (addr == NULL) {
			if (rand_dist == rand_dist_type::NORM)
				addr = new rand_gen_wrapper(type.create_randn_gen(var1, var2));
			else if (rand_dist == rand_dist_type::UNIF)
				addr = new rand_gen_wrapper(type.create_randu_gen(var1, var2));
			else
				assert(0);
			int ret = pthread_setspecific(gen_key, addr);
			assert(ret == 0);
		}
		rand_gen_wrapper *wrapper = (rand_gen_wrapper *) addr;
		return wrapper->get_gen();
	}

	static void destroy_rand_gen(void *gen) {
		rand_gen_wrapper *wrapper = (rand_gen_wrapper *) gen;
		delete wrapper;
		printf("destroy rand gen\n");
	}
public:
	rand_init(const scalar_variable &_var1, const scalar_variable &_var2,
			rand_dist_type rand_dist): type(_var1.get_type()), var1(
				_var1), var2(_var2) {
		int ret = pthread_key_create(&gen_key, destroy_rand_gen);
		this->rand_dist = rand_dist;
		assert(ret == 0);
	}

	~rand_init() {
		pthread_key_delete(gen_key);
	}

	virtual void set(void *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		get_rand_gen().gen(arr, num_eles);
	}
	virtual const scalar_type &get_type() const {
		return get_rand_gen().get_type();
	}
};

}

void matrix_store::init_randu(const scalar_variable &min,
		const scalar_variable &max)
{
	assert(min.get_type() == max.get_type());
	set_data(rand_init(min, max, rand_init::rand_dist_type::UNIF));
}

void matrix_store::init_randn(const scalar_variable &mean,
		const scalar_variable &var)
{
	assert(mean.get_type() == var.get_type());
	set_data(rand_init(mean, var, rand_init::rand_dist_type::NORM));
}

}

}
