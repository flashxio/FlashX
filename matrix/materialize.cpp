/*
 * Copyright 2016 Open Connectome Project (http://openconnecto.me)
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

#include <unordered_set>
#include <unordered_map>

#include "common.h"
#include "thread.h"

#include "local_mem_buffer.h"
#include "matrix_config.h"
#include "matrix_stats.h"
#include "local_matrix_store.h"
#include "materialize.h"
#include "dense_matrix.h"
#include "EM_object.h"
#include "EM_dense_matrix.h"
#include "mapply_matrix_store.h"
#include "sink_matrix.h"

typedef std::vector<size_t> mat_id_set;

namespace std
{

template<>
struct hash<mat_id_set>
{
	size_t operator() (const mat_id_set &set) const {
		size_t id = 0;
		for (size_t i = 0; i < set.size(); i++)
			id = id * 10 + set[i];
		return id;
	}
};

template<>
struct equal_to<mat_id_set>
{
	bool operator()(const mat_id_set &s1, const mat_id_set &s2) const {
		if (s1.size() != s2.size())
			return false;

		for (size_t i = 0; i < s1.size(); i++)
			if (s1[i] != s2[i])
				return false;
		return true;
	}
};

}

namespace fm
{

namespace detail
{

namespace
{

class mapply_task: public thread_task
{
	const std::vector<matrix_store::const_ptr> &mats;
	const std::vector<matrix_store::ptr> &out_mats;
	size_t portion_idx;
	const portion_mapply_op &op;
	bool one_portion;
public:
	mapply_task(const std::vector<matrix_store::const_ptr> &_mats,
			size_t portion_idx, const portion_mapply_op &_op,
			const std::vector<matrix_store::ptr> &_out_mats): mats(
				_mats), out_mats(_out_mats), op(_op) {
		this->portion_idx = portion_idx;
		this->one_portion = mats.front()->get_num_portions() == 1;
	}

	void run();
};

void mapply_task::run()
{
	std::vector<detail::local_matrix_store::const_ptr> local_stores(
			mats.size());
	std::vector<detail::local_matrix_store::ptr> local_out_stores(
			out_mats.size());
	int node_id = thread::get_curr_thread()->get_node_id();
	for (size_t j = 0; j < mats.size(); j++) {
		local_stores[j] = mats[j]->get_portion(portion_idx);
#ifdef MATRIX_DEBUG
		if (local_stores[j] && local_stores[j]->get_node_id() >= 0
				&& !one_portion)
			assert(node_id == local_stores[j]->get_node_id());
#endif
	}
	for (size_t j = 0; j < out_mats.size(); j++) {
		local_out_stores[j] = out_mats[j]->get_portion(portion_idx);
#ifdef MATRIX_DEBUG
		if (local_out_stores[j] && local_out_stores[j]->get_node_id() >= 0
				&& !one_portion)
			assert(node_id == local_out_stores[j]->get_node_id());
#endif
	}

	if (local_out_stores.empty())
		op.run(local_stores);
	else if (local_out_stores.size() == 1)
		op.run(local_stores, *local_out_stores[0]);
	else
		op.run(local_stores, local_out_stores);
}

class mem_task_queue
{
public:
	typedef std::shared_ptr<mem_task_queue> ptr;

	static const size_t INVALID_TASK = std::numeric_limits<size_t>::max();

	virtual size_t get_portion_idx() = 0;
	virtual int get_node_id() const = 0;
};

/*
 * This dispatches the computation tasks for all threads in a machine,
 * so it is shared by all threads.
 */
class smp_task_queue: public mem_task_queue
{
	std::atomic<size_t> curr_portion;
	size_t num_portions;
public:
	smp_task_queue(size_t num_portions) {
		this->num_portions = num_portions;
		this->curr_portion = 0;
	}

	virtual size_t get_portion_idx() {
		size_t portion_idx = curr_portion.fetch_add(1);
		if (portion_idx >= num_portions)
			return INVALID_TASK;
		else
			return portion_idx;
	}

	virtual int get_node_id() const {
		return -1;
	}
};

/*
 * This dispatches the computation tasks for the portions in a specific
 * NUMA node. All threads in a NUMA node share the same task queue.
 */
class numa_task_queue: public mem_task_queue
{
	matrix_store::const_ptr numa_mat;
	size_t num_portions;
	std::atomic<size_t> curr_portion;
	int curr_node_id;
public:
	numa_task_queue(matrix_store::const_ptr mat, int node_id) {
		this->numa_mat = mat;
		this->num_portions = mat->get_num_portions();
		this->curr_portion = 0;
		this->curr_node_id = node_id;
	}

	virtual size_t get_portion_idx() {
		while (true) {
			// TODO this doesn't scale to a large number of NUMA nodes.
			size_t portion_idx = curr_portion.fetch_add(1);
			// If we have run out of portions, we need to indicate that
			// we have processed all tasks.
			if (portion_idx >= num_portions)
				return INVALID_TASK;
			// We only return the portion in the current node.
			int node_id = numa_mat->get_portion_node_id(portion_idx);
			if (node_id == curr_node_id)
				return portion_idx;
		}
	}

	virtual int get_node_id() const {
		return curr_node_id;
	}
};

class mem_worker_task: public thread_task
{
	std::vector<matrix_store::const_ptr> mats;
	std::vector<matrix_store::ptr> out_mats;
	const portion_mapply_op &op;
	mem_task_queue::ptr task_queue;
public:
	mem_worker_task(const std::vector<matrix_store::const_ptr> mats,
			const portion_mapply_op &_op,
			const std::vector<matrix_store::ptr> out_mats,
			mem_task_queue::ptr task_queue): op(_op) {
		this->mats = mats;
		this->out_mats = out_mats;
		this->task_queue = task_queue;
	}

	void run() {
		size_t portion_idx;
		while ((portion_idx = task_queue->get_portion_idx())
				!= mem_task_queue::INVALID_TASK) {
			mapply_task task(mats, portion_idx, op, out_mats);
			task.run();
		}
	}
};

static size_t cal_min_portion_size(
		const std::vector<matrix_store::const_ptr> &mats1,
		const std::vector<matrix_store::ptr> &mats2)
{
	assert(mats1.size() > 0);
	if (mats1[0]->is_wide()) {
		size_t min_portion_size = mats1[0]->get_portion_size().second;
		for (size_t i = 1; i < mats1.size(); i++)
			min_portion_size = std::min(min_portion_size,
					mats1[i]->get_portion_size().second);
		for (size_t i = 0; i < mats2.size(); i++) {
			min_portion_size = std::min(min_portion_size,
					mats2[i]->get_portion_size().second);
			assert(mats2[i]->get_portion_size().second % min_portion_size == 0);
		}
		for (size_t i = 0; i < mats1.size(); i++)
			assert(mats1[i]->get_portion_size().second % min_portion_size == 0);
		return min_portion_size;
	}
	else {
		size_t min_portion_size = mats1[0]->get_portion_size().first;
		for (size_t i = 1; i < mats1.size(); i++)
			min_portion_size = std::min(min_portion_size,
					mats1[i]->get_portion_size().first);
		for (size_t i = 0; i < mats2.size(); i++) {
			min_portion_size = std::min(min_portion_size,
					mats2[i]->get_portion_size().first);
			assert(mats2[i]->get_portion_size().first % min_portion_size == 0);
		}
		for (size_t i = 0; i < mats1.size(); i++)
			assert(mats1[i]->get_portion_size().first % min_portion_size == 0);
		return min_portion_size;
	}
}

/*
 * This dispatcher issues I/O to access the same portion of all dense matrices
 * simultaneously. This may have good I/O performance, but may consume a lot
 * of memory when it runs for a large group of dense matrices.
 */
class EM_mat_mapply_par_dispatcher: public detail::EM_portion_dispatcher
{
	std::vector<matrix_store::const_ptr> mats;
	std::vector<matrix_store::ptr> res_mats;
	std::vector<matrix_stream::ptr> res_mat_streams;
	std::vector<matrix_store::const_ptr> EM_mats;
	portion_mapply_op::const_ptr op;
	size_t min_portion_size;
public:
	EM_mat_mapply_par_dispatcher(
			const std::vector<matrix_store::const_ptr> &mats,
			const std::vector<matrix_store::ptr> &res_mats,
			portion_mapply_op::const_ptr op, size_t tot_len,
			size_t portion_size);
	~EM_mat_mapply_par_dispatcher();

	virtual void create_task(off_t global_start, size_t length);
};

/*
 * This collects all the portions in a partition that are required by
 * an operation and are ready in memory.
 */
class collected_portions
{
	// The output matrices.
	std::vector<matrix_stream::ptr> res_mat_streams;
	// The portions with partial result.
	// This is only used when the operation is aggregation.
	local_matrix_store::ptr res_portion;
	// The portions with data ready.
	// This is only used when the operation isn't aggregation.
	std::vector<local_matrix_store::const_ptr> ready_portions;

	// The number of portions that are ready.
	size_t num_ready;
	// The number of portions that are required for the computation.
	size_t num_required;

	// The location of the portions.
	off_t global_start;
	size_t length;

	portion_mapply_op::const_ptr op;
public:
	typedef std::shared_ptr<collected_portions> ptr;

	collected_portions(const std::vector<matrix_stream::ptr> &res_mat_streams,
			portion_mapply_op::const_ptr op, size_t num_required,
			off_t global_start, size_t length) {
		this->res_mat_streams = res_mat_streams;
		this->global_start = global_start;
		this->length = length;
		this->op = op;
		this->num_required = num_required;
		this->num_ready = 0;
	}

	off_t get_global_start() const {
		return global_start;
	}

	size_t get_length() const {
		return length;
	}

	void run_all_portions();

	void add_ready_portion(local_matrix_store::const_ptr portion);

	bool is_complete() const {
		return num_required == num_ready;
	}
};

/*
 * This dispatcher accesses one portion of a dense matrix at a time in a thread,
 * although I/O is still performed asynchronously. This is particularly useful
 * when we compute on a large group of dense matrices. Users can split the large
 * group and build a hierarchy to compute on the dense matrices. In this hierarchy,
 * we can use this dispatcher to limit the matrices that are being accessed
 * at the same time to reduce memory consumption.
 */
class EM_mat_mapply_serial_dispatcher: public detail::EM_portion_dispatcher
{
	/*
	 * This contains the data structure for fetching data in a partition
	 * across all dense matrices in an operation.
	 * The first member contains the matrices whose portions haven't been
	 * fetched; the second member contains the portions that have been
	 * successfully fetched.
	 */
	typedef std::pair<std::deque<matrix_store::const_ptr>,
			collected_portions::ptr> part_state_t;

	std::vector<matrix_store::const_ptr> mats;
	std::vector<matrix_store::ptr> res_mats;
	std::vector<matrix_stream::ptr> res_mat_streams;

	// These maintain per-thread states.
	std::vector<part_state_t> part_states;
	portion_mapply_op::const_ptr op;
	size_t min_portion_size;
public:
	EM_mat_mapply_serial_dispatcher(
			const std::vector<matrix_store::const_ptr> &mats,
			const std::vector<matrix_store::ptr> &res_mats,
			portion_mapply_op::const_ptr op, size_t tot_len,
			size_t portion_size): detail::EM_portion_dispatcher(tot_len,
				portion_size) {
		this->mats = mats;
		this->res_mats = res_mats;
		this->op = op;
		this->min_portion_size = cal_min_portion_size(mats, res_mats);
		part_states.resize(detail::mem_thread_pool::get_global_num_threads());
		res_mat_streams.resize(res_mats.size());
		for (size_t i = 0; i < res_mats.size(); i++)
			res_mat_streams[i] = matrix_stream::create(res_mats[i]);
	}

	virtual void create_task(off_t global_start, size_t length);
	virtual bool issue_task();
};

class mapply_portion_compute: public portion_compute
{
	std::vector<detail::local_matrix_store::const_ptr> local_stores;
	std::vector<matrix_stream::ptr> mat_streams;
	size_t num_required_reads;
	size_t num_reads;
	const portion_mapply_op &op;
public:
	mapply_portion_compute(const std::vector<matrix_stream::ptr> mat_streams,
			const portion_mapply_op &_op): op(_op) {
		this->mat_streams = mat_streams;
		this->num_required_reads = 0;
		this->num_reads = 0;
	}

	void set_buf(
			const std::vector<detail::local_matrix_store::const_ptr> &stores) {
		this->local_stores = stores;
	}
	void set_EM_parts(size_t num_EM_parts) {
		this->num_required_reads = num_EM_parts;
	}

	virtual void run(char *buf, size_t size);

	void run_complete();
};

static inline local_matrix_store::ptr create_local_buf_matrix(
		const matrix_store &to_mat, off_t start_row, off_t start_col,
		size_t num_rows, size_t num_cols)
{
	if (to_mat.store_layout() == matrix_layout_t::L_COL)
		return local_matrix_store::ptr(new local_buf_col_matrix_store(
					start_row, start_col, num_rows, num_cols,
					to_mat.get_type(), -1));
	else
		return local_matrix_store::ptr(new local_buf_row_matrix_store(
					start_row, start_col, num_rows, num_cols,
					to_mat.get_type(), -1));
}

void mapply_portion_compute::run_complete()
{
	assert(!local_stores.empty());
	const detail::local_matrix_store &first_mat = *local_stores.front();
	std::vector<local_matrix_store::ptr> local_res(mat_streams.size());
	for (size_t i = 0; i < mat_streams.size(); i++) {
		size_t res_num_rows;
		size_t res_num_cols;
		if (mat_streams[i]->get_mat().is_wide()) {
			res_num_rows = mat_streams[i]->get_mat().get_num_rows();
			res_num_cols = first_mat.get_num_cols();
		}
		else {
			res_num_rows = first_mat.get_num_rows();
			res_num_cols = mat_streams[i]->get_mat().get_num_cols();
		}
		local_res[i] = create_local_buf_matrix(mat_streams[i]->get_mat(),
				first_mat.get_global_start_row(), first_mat.get_global_start_col(),
				res_num_rows, res_num_cols);
	}
	if (local_res.empty())
		op.run(local_stores);
	else if (local_res.size() == 1)
		op.run(local_stores, *local_res[0]);
	else
		op.run(local_stores, local_res);
	for (size_t i = 0; i < mat_streams.size(); i++)
		mat_streams[i]->write_async(local_res[i],
				local_res[i]->get_global_start_row(),
				local_res[i]->get_global_start_col());
}

void mapply_portion_compute::run(char *buf, size_t size)
{
	assert(!local_stores.empty());
	num_reads++;
	if (num_required_reads == num_reads)
		run_complete();
}

static size_t cal_task_size(const std::vector<matrix_store::const_ptr> &mats,
		const std::vector<matrix_store::ptr> &res_mats)
{
	size_t max_num_bytes = 0;
	for (size_t i = 0; i < mats.size(); i++) {
		auto psize = mats[i]->get_portion_size();
		max_num_bytes = std::max(max_num_bytes,
				psize.first * psize.second * mats[i]->get_type().get_size());
	}
	for (size_t i = 0; i < res_mats.size(); i++) {
		auto psize = res_mats[i]->get_portion_size();
		max_num_bytes = std::max(max_num_bytes,
				psize.first * psize.second * res_mats[i]->get_type().get_size());
	}
	size_t num_portions = div_ceil(
			(size_t) safs::params.get_RAID_block_size() * PAGE_SIZE,
			max_num_bytes);
	size_t ret = 1;
	for (; ret < num_portions; ret *= 2);
	return ret;
}

static inline size_t get_reserved_portions()
{
	return mem_thread_pool::get_global_num_threads() * 10;
}

EM_mat_mapply_par_dispatcher::EM_mat_mapply_par_dispatcher(
		const std::vector<matrix_store::const_ptr> &mats,
		const std::vector<matrix_store::ptr> &res_mats,
		portion_mapply_op::const_ptr op, size_t tot_len,
		size_t portion_size): detail::EM_portion_dispatcher(tot_len,
			portion_size,
			// This is the number of portions we reserve to be processed
			// one at a time.
			get_reserved_portions(),
			// This is the number of portions in a regular task.
			cal_task_size(mats, res_mats))
{
	this->mats = mats;
	this->res_mats = res_mats;
	this->op = op;
	this->min_portion_size = cal_min_portion_size(mats, res_mats);

	// Set prefetch on each EM matrix.
	size_t num_reserved = get_reserved_portions();
	for (size_t i = 0; i < mats.size(); i++) {
		if (!mats[i]->is_in_mem()) {
			size_t prefetch_end = 0;
			if (mats[i]->get_num_portions() > num_reserved) {
				prefetch_end = mats[i]->get_num_portions() - num_reserved;
				// We want the prefetch range to be aligned with the prefetch
				// size.
				prefetch_end
					= (prefetch_end / get_task_size()) * get_task_size();
			}
			const_cast<matrix_store &>(*mats[i]).set_prefetches(get_task_size(),
					std::pair<size_t, size_t>(0, prefetch_end));
			EM_mats.push_back(mats[i]);
		}
	}
	res_mat_streams.resize(res_mats.size());
	for (size_t i = 0; i < res_mats.size(); i++) {
		if (!res_mats[i]->is_in_mem()) {
			size_t prefetch_end = 0;
			if (mats[i]->get_num_portions() > num_reserved) {
				prefetch_end = mats[i]->get_num_portions() - num_reserved;
				// We want the prefetch range to be aligned with the prefetch
				// size.
				prefetch_end
					= (prefetch_end / get_task_size()) * get_task_size();
			}
			const_cast<matrix_store &>(*res_mats[i]).set_prefetches(
					get_task_size(), std::pair<size_t, size_t>(0, prefetch_end));
			EM_mats.push_back(res_mats[i]);
		}
		res_mat_streams[i] = matrix_stream::create(res_mats[i]);
		assert(res_mat_streams[i]);
	}
}

EM_mat_mapply_par_dispatcher::~EM_mat_mapply_par_dispatcher()
{
	for (size_t i = 0; i < EM_mats.size(); i++)
		const_cast<matrix_store &>(*EM_mats[i]).set_prefetches(1,
				std::pair<size_t, size_t>(0, 0));
}

/*
 * This method is invoked in each worker thread.
 * It reads data portions from each matrix and perform computation.
 * To be able to adapt to matrices with different shapes, a task may read
 * multiple portions from a matrix.
 */
void EM_mat_mapply_par_dispatcher::create_task(off_t global_start,
		size_t length)
{
	// We fetch the portions using the minimum portion size among
	// the matrices. The idea is to reduce the amount of data cached
	// in virtual matrices.
	for (size_t local_start = 0; local_start < length;
			local_start += min_portion_size) {
		size_t local_length = std::min(min_portion_size, length - local_start);
		std::vector<detail::local_matrix_store::const_ptr> local_stores(
				mats.size());
		mapply_portion_compute *mapply_compute = new mapply_portion_compute(
				res_mat_streams, *op);
		mapply_portion_compute::ptr compute(mapply_compute);
		size_t num_EM_parts = 0;
		for (size_t j = 0; j < local_stores.size(); j++) {
			size_t global_start_row, global_start_col, num_rows, num_cols;
			if (mats[j]->is_wide()) {
				global_start_row = 0;
				global_start_col = global_start + local_start;
				num_rows = mats[j]->get_num_rows();
				num_cols = local_length;
			}
			else {
				global_start_row = global_start + local_start;
				global_start_col = 0;
				num_rows = local_length;
				num_cols = mats[j]->get_num_cols();
			}
			async_cres_t res = mats[j]->get_portion_async(global_start_row,
					global_start_col, num_rows, num_cols, compute);
			if (!res.first)
				num_EM_parts++;
			local_stores[j] = res.second;
		}
		mapply_compute->set_buf(local_stores);
		mapply_compute->set_EM_parts(num_EM_parts);
		// When all input parts are in memory or have been cached, we need to
		// run the portion compute manually by ourselves.
		if (num_EM_parts == 0)
			mapply_compute->run_complete();
	}
}

void collected_portions::add_ready_portion(local_matrix_store::const_ptr portion)
{
	num_ready++;
	assert(num_ready <= num_required);

	if (op->is_agg() && res_portion) {
		std::vector<local_matrix_store::const_ptr> local_stores(2);
		local_stores[0] = res_portion;
		local_stores[1] = portion;
		// We store the partial result in a single portion, regardless of
		// the number of output matrices we want to generate eventually.
		op->run(local_stores, *res_portion);
	}
	else if (op->is_agg()) {
		if (op->get_out_num_rows() > 0 && op->get_out_num_cols() > 0) {
			// Regardless of the number of output matrices we want to generate
			// eventually, we only use one matrix portion to store
			// the intermediate aggregation result.
			size_t start_row, start_col, num_rows, num_cols;
			if (res_mat_streams.front()->get_mat().is_wide()) {
				start_row = 0;
				start_col = global_start;
				num_rows = op->get_out_num_rows();
				num_cols = length;
			}
			else {
				start_row = global_start;
				start_col = 0;
				num_rows = length;
				num_cols = op->get_out_num_cols();
			}

			if (res_mat_streams.front()->get_mat().store_layout()
					== matrix_layout_t::L_COL)
				res_portion = local_matrix_store::ptr(
						new local_buf_col_matrix_store(start_row, start_col,
							num_rows, num_cols, op->get_output_type(),
							portion->get_node_id()));
			else
				res_portion = local_matrix_store::ptr(
						new local_buf_row_matrix_store(start_row, start_col,
							num_rows, num_cols, op->get_output_type(),
							portion->get_node_id()));

			std::vector<local_matrix_store::const_ptr> local_stores(1);
			local_stores[0] = portion;
			// We rely on the user-defined function to copy the first portion
			// to the partial result portion.
			op->run(local_stores, *res_portion);
		}
		else {
			std::vector<local_matrix_store::const_ptr> local_stores(1);
			local_stores[0] = portion;
			// We rely on the user-defined function to copy the first portion
			// to the partial result portion.
			op->run(local_stores);
		}
	}
	else {
		// This forces the portion of a mapply matrix to materialize and
		// release the data in the underlying portion.
		portion->materialize_self();
		ready_portions.push_back(portion);
	}
}

void collected_portions::run_all_portions()
{
	assert(is_complete());
	// If this is an aggregation operation and no result portion was generated,
	// return now.
	if (op->is_agg() && res_portion == NULL)
		return;

	mapply_portion_compute compute(res_mat_streams, *op);
	if (res_portion) {
		std::vector<detail::local_matrix_store::const_ptr> stores(1);
		stores[0] = res_portion;
		compute.set_buf(stores);
		compute.run_complete();
		res_portion = NULL;
	}
	else {
		compute.set_buf(ready_portions);
		compute.run_complete();

		// We don't need these portions now. They should be free'd.
		ready_portions.clear();
	}
}

/*
 * This class reads the EM portions one at a time. After it reads all portions,
 * it invokes mapply_portion_compute.
 */
class serial_read_portion_compute: public portion_compute
{
	collected_portions::ptr collected;

	// The portion that is being read.
	local_matrix_store::const_ptr pending_portion;
public:
	serial_read_portion_compute(collected_portions::ptr collected) {
		this->collected = collected;
	}

	virtual void run(char *buf, size_t size) {
		collected->add_ready_portion(pending_portion);
		pending_portion = NULL;
		if (collected->is_complete())
			collected->run_all_portions();
		collected = NULL;
	}

	static void fetch_portion(std::deque<matrix_store::const_ptr> &mats,
			collected_portions::ptr collected);
};

// TODO I need to copy the vectors multiple times. Is this problematic?
void serial_read_portion_compute::fetch_portion(
		std::deque<matrix_store::const_ptr> &mats,
		collected_portions::ptr collected)
{
	serial_read_portion_compute *_compute
		= new serial_read_portion_compute(collected);
	serial_read_portion_compute::ptr compute(_compute);
	size_t global_start_row, global_start_col, num_rows, num_cols;
	while (!mats.empty()) {
		matrix_store::const_ptr mat = mats.front();
		mats.pop_front();
		if (mat->is_wide()) {
			global_start_row = 0;
			global_start_col = collected->get_global_start();
			num_rows = mat->get_num_rows();
			num_cols = collected->get_length();
		}
		else {
			global_start_row = collected->get_global_start();
			global_start_col = 0;
			num_rows = collected->get_length();
			num_cols = mat->get_num_cols();
		}
		async_cres_t res = mat->get_portion_async(global_start_row,
				global_start_col, num_rows, num_cols, compute);
		if (!res.first) {
			_compute->pending_portion = res.second;
			break;
		}
		else
			collected->add_ready_portion(res.second);
	}
	if (collected->is_complete()) {
		// All portions are ready, we can perform computation now.
		collected->run_all_portions();
	}
}

/*
 * This method is invoked in each worker thread.
 */
void EM_mat_mapply_serial_dispatcher::create_task(off_t global_start,
		size_t length)
{
	int thread_id = mem_thread_pool::get_curr_thread_id();
	assert(part_states[thread_id].first.empty());
	assert(part_states[thread_id].second == NULL);

	// Here we don't use the minimum portion size.
	collected_portions::ptr collected(new collected_portions(res_mat_streams,
				op, mats.size(), global_start, length));
	part_states[thread_id].second = collected;
	part_states[thread_id].first.insert(part_states[thread_id].first.end(),
			mats.begin(), mats.end());
	serial_read_portion_compute::fetch_portion(part_states[thread_id].first,
			collected);
}

bool EM_mat_mapply_serial_dispatcher::issue_task()
{
	int thread_id = mem_thread_pool::get_curr_thread_id();
	// If the data of the matrices in the current partition hasn't been
	// fetched, we fetch them first.
	if (!part_states[thread_id].first.empty()) {
		assert(part_states[thread_id].second != NULL);
		serial_read_portion_compute::fetch_portion(
				part_states[thread_id].first, part_states[thread_id].second);
		return true;
	}

	// If we are ready to access the next partition, we reset the pointer
	// to the collection of the portions in the current partition.
	part_states[thread_id].second = NULL;
	return EM_portion_dispatcher::issue_task();
}

}

matrix_store::ptr __mapply_portion(
		const std::vector<matrix_store::const_ptr> &mats,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout,
		bool par_access)
{
	// As long as one of the input matrices is in external memory, we output
	// an EM matrix.
	bool out_in_mem = mats.front()->is_in_mem();
	for (size_t i = 1; i < mats.size(); i++)
		out_in_mem = out_in_mem && mats[i]->is_in_mem();

	int num_nodes = -1;
	if (out_in_mem) {
		num_nodes = mats[0]->get_num_nodes();
		for (size_t i = 1; i < mats.size(); i++)
			num_nodes = std::max(num_nodes, mats[i]->get_num_nodes());
	}
	std::vector<matrix_store::ptr> out_mats;
	if (op->get_out_num_rows() > 0 && op->get_out_num_cols() > 0) {
		detail::matrix_store::ptr res = detail::matrix_store::create(
				op->get_out_num_rows(), op->get_out_num_cols(),
				out_layout, op->get_output_type(), num_nodes, out_in_mem);
		if (res == NULL)
			return matrix_store::ptr();

		out_mats.push_back(res);
	}
	bool ret = __mapply_portion(mats, op, out_mats, par_access);
	if (ret && out_mats.size() == 1)
		return out_mats[0];
	else
		return matrix_store::ptr();
}

matrix_store::ptr __mapply_portion(
		const std::vector<matrix_store::const_ptr> &mats,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout,
		bool out_in_mem, int out_num_nodes, bool par_access)
{
	std::vector<matrix_store::ptr> out_mats;
	if (op->get_out_num_rows() > 0 && op->get_out_num_cols() > 0) {
		detail::matrix_store::ptr res = detail::matrix_store::create(
				op->get_out_num_rows(), op->get_out_num_cols(),
				out_layout, op->get_output_type(), out_num_nodes, out_in_mem);
		if (res == NULL)
			return matrix_store::ptr();

		out_mats.push_back(res);
	}
	bool ret = __mapply_portion(mats, op, out_mats, par_access);
	if (ret && out_mats.size() == 1)
		return out_mats[0];
	else
		return matrix_store::ptr();
}

// This function computes the result of mapply.
// We allow to not output a matrix.
bool __mapply_portion(
		const std::vector<matrix_store::const_ptr> &mats,
		portion_mapply_op::const_ptr op,
		const std::vector<matrix_store::ptr> &out_mats, bool par_access)
{
	std::vector<matrix_store::const_ptr> numa_mats;
	bool out_in_mem;
	int out_num_nodes;
	if (out_mats.empty()) {
		out_in_mem = true;
		out_num_nodes = -1;
	}
	else {
		out_in_mem = out_mats.front()->is_in_mem();
		out_num_nodes = out_mats.front()->get_num_nodes();
		if (out_num_nodes > 0)
			numa_mats.push_back(out_mats.front());
		detail::matrix_stats.inc_write_bytes(
				out_mats[0]->get_num_rows() * out_mats[0]->get_num_cols()
				* out_mats[0]->get_entry_size(), out_in_mem);
	}
	for (size_t i = 1; i < out_mats.size(); i++) {
		assert(out_in_mem == out_mats[i]->is_in_mem());
		assert(out_num_nodes == out_mats[i]->get_num_nodes());
		detail::matrix_stats.inc_write_bytes(
				out_mats[i]->get_num_rows() * out_mats[i]->get_num_cols()
				* out_mats[i]->get_entry_size(), out_in_mem);
		if (out_mats[i]->get_num_nodes() > 0)
			numa_mats.push_back(out_mats[i]);
	}
	assert(mats.size() >= 1);
	assert(op->is_success());

	// Collect all NUMA matrices in the input matrices.
	for (size_t i = 0; i < mats.size(); i++) {
		if (mats[i]->get_num_nodes() > 0)
			numa_mats.push_back(mats[i]);
	}

	bool all_in_mem = mats[0]->is_in_mem();
	size_t num_chunks = mats.front()->get_num_portions();
	std::pair<size_t, size_t> first_size = mats.front()->get_portion_size();
	size_t tot_len;
	size_t portion_size;
	bool in_square = mats.front()->get_num_rows() == mats.front()->get_num_cols();
	// The input matrix is wide.
	if (mats.front()->is_wide()
			// The input matrix is square and the output matrix is wide.
			|| (in_square && op->is_wide())) {
		tot_len = mats.front()->get_num_cols();
		portion_size = first_size.second;
		if (op->get_out_num_cols() > 0)
			assert(op->get_out_num_cols() == mats.front()->get_num_cols());
		for (size_t i = 1; i < mats.size(); i++) {
			portion_size = std::max(portion_size,
					mats[i]->get_portion_size().second);
			assert(mats[i]->get_num_cols() == tot_len);
			all_in_mem = all_in_mem && mats[i]->is_in_mem();
		}
	}
	// There are two cases that we are here.
	// The input matrix is tall.
	// The input matrix is square and the output matrix isn't wide.
	else {
		tot_len = mats.front()->get_num_rows();
		portion_size = first_size.first;
		if (op->get_out_num_rows() > 0)
			assert(op->get_out_num_rows() == mats.front()->get_num_rows());
		for (size_t i = 1; i < mats.size(); i++) {
			portion_size = std::max(portion_size,
					mats[i]->get_portion_size().first);
			assert(mats[i]->get_num_rows() == tot_len);
			all_in_mem = all_in_mem && mats[i]->is_in_mem();
		}
	}
	all_in_mem = all_in_mem && out_in_mem;

	// We need to test if all NUMA matrices distribute data to NUMA nodes
	// with the same mapping. It's kind of difficult to test it.
	// Right now the system has only one mapper. As long as a NUMA matrix
	// is distributed to the same number of NUMA nodes, they should be
	// distributed with the same mapping.
	for (size_t i = 1; i < numa_mats.size(); i++)
		assert(numa_mats[i]->get_num_nodes() == numa_mats[0]->get_num_nodes());

	// If all matrices are in memory and all matrices have only one portion,
	// we should perform computation in the local thread.
	if (all_in_mem && num_chunks == 1) {
		mapply_task task(mats, 0, *op, out_mats);
		task.run();
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
	else if (all_in_mem) {
		detail::mem_thread_pool::ptr threads
			= detail::mem_thread_pool::get_global_mem_threads();
		std::vector<mem_task_queue::ptr> task_queues;
		// If there aren't NUMA matrices, there will be only one queue.
		if (numa_mats.empty()) {
			task_queues.resize(1);
			task_queues[0] = mem_task_queue::ptr(
					new smp_task_queue(mats[0]->get_num_portions()));
		}
		// Otherwise, each NUMA node has a task queue.
		else {
			int num_nodes = numa_mats.front()->get_num_nodes();
			assert(num_nodes > 0);
			task_queues.resize(num_nodes);
			for (int i = 0; i < num_nodes; i++)
				task_queues[i] = mem_task_queue::ptr(
						new numa_task_queue(numa_mats.front(), i));
		}
		for (size_t i = 0; i < threads->get_num_threads(); i++) {
			int node_id = i % threads->get_num_nodes();
			mem_task_queue::ptr queue;
			if (task_queues.size() == 1)
				queue = task_queues[0];
			else {
				assert((size_t) node_id < task_queues.size());
				queue = task_queues[node_id];
				assert(queue->get_node_id() == node_id);
			}
			threads->process_task(node_id, new mem_worker_task(mats, *op,
						out_mats, queue));
		}
		threads->wait4complete();
	}
	else {
		mem_thread_pool::ptr threads = mem_thread_pool::get_global_mem_threads();
		detail::EM_portion_dispatcher::ptr dispatcher;

		// Start streams for external-memory result matrices.
		std::vector<EM_matrix_store::ptr> em_outs;
		for (size_t i = 0; i < out_mats.size(); i++) {
			if (out_mats[i]->is_in_mem())
				continue;
			EM_matrix_store::ptr em_out
				= std::dynamic_pointer_cast<EM_matrix_store>(out_mats[i]);
			assert(em_out);
			em_outs.push_back(em_out);
			em_out->start_stream();
		}

		if (par_access)
			dispatcher = detail::EM_portion_dispatcher::ptr(
					new EM_mat_mapply_par_dispatcher(mats, out_mats, op,
						tot_len, EM_matrix_store::CHUNK_SIZE));
		else
			dispatcher = detail::EM_portion_dispatcher::ptr(
					new EM_mat_mapply_serial_dispatcher(mats, out_mats, op,
						tot_len, EM_matrix_store::CHUNK_SIZE));
		for (size_t i = 0; i < threads->get_num_threads(); i++) {
			io_worker_task *task = new io_worker_task(dispatcher);
			for (size_t j = 0; j < mats.size(); j++) {
				if (!mats[j]->is_in_mem()) {
					const EM_object *obj
						= dynamic_cast<const EM_object *>(mats[j].get());
					assert(obj);
					task->register_EM_obj(const_cast<EM_object *>(obj));
				}
			}
			for (size_t j = 0; j < out_mats.size(); j++) {
				if (!out_mats[j]->is_in_mem())
					task->register_EM_obj(
							dynamic_cast<EM_object *>(out_mats[j].get()));
			}
			threads->process_task(i % threads->get_num_nodes(), task);
		}
		threads->wait4complete();

		// End the streams for external-memory result matrices.
		for (size_t i = 0; i < em_outs.size(); i++)
			em_outs[i]->end_stream();
	}
	return op->is_success();
}

matrix_store::ptr __mapply_portion_virtual(
		const std::vector<matrix_store::const_ptr> &stores,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout,
		bool par_access)
{
	mapply_matrix_store *store = new mapply_matrix_store(stores, op, out_layout);
	store->set_par_access(par_access);
	return matrix_store::ptr(store);
}

dense_matrix::ptr mapply_portion(
		const std::vector<dense_matrix::const_ptr> &mats,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout,
		bool par_access)
{
	std::vector<matrix_store::const_ptr> stores(mats.size());
	for (size_t i = 0; i < stores.size(); i++)
		stores[i] = mats[i]->get_raw_store();
	mapply_matrix_store *store = new mapply_matrix_store(stores, op, out_layout);
	store->set_par_access(par_access);
	matrix_store::const_ptr ret(store);
	return dense_matrix::create(ret);
}

}

namespace
{

class materialize_mapply_op: public detail::portion_mapply_op
{
	bool is_wide;
	bool _agg;
public:
	// The type doesn't really matter.
	materialize_mapply_op(const scalar_type &type,
			bool is_wide, bool _agg): detail::portion_mapply_op(0, 0, type) {
		this->is_wide = is_wide;
		this->_agg = _agg;
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		throw unsupported_exception(
				"Don't support transpose of groupby_row_mapply_op");
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const {
		if (is_wide)
			detail::materialize_wide(ins);
		else
			detail::materialize_tall(ins);
	}

	virtual bool is_agg() const {
		return _agg;
	}

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		throw unsupported_exception(
				"Don't support to_string of groupby_row_mapply_op");
	}
};

/*
 * This class helps us to identify the underlying matrices that a virtual
 * matrix relies on. We'll keep two virtual matrices together if
 * the underlying matrices of a virtual matrix are the subset of the other
 * virtual matrices. Eventually, a collection of virtual matrices in this
 * class will be materialized together to reduce I/O.
 */
class underlying_mat_set
{
	std::unordered_set<size_t> under_mat_set;
	mat_id_set under_mats;
	std::unordered_map<size_t, detail::virtual_matrix_store::const_ptr> vmats;

	underlying_mat_set(detail::virtual_matrix_store::const_ptr mat) {
		// TODO If this is a set of mixed IM and EM matrices, we might want to
		// consider only EM matrices.
		auto under_map = mat->get_underlying_mats();
		for (auto it = under_map.begin(); it != under_map.end(); it++) {
			// Some of the virtual matrices (e.g., one_val_matrix_store and
			// set_data_matrix_store) don't store data physically. We should
			// ignore these matrices.
			if (it->second > 0)
				this->under_mats.push_back(it->first);
		}
		std::sort(under_mats.begin(), under_mats.end());
		this->under_mat_set.insert(under_mats.begin(), under_mats.end());
		vmats.insert(std::pair<size_t, detail::virtual_matrix_store::const_ptr>(
					mat->get_data_id(), mat));
	}
public:
	typedef std::shared_ptr<underlying_mat_set> ptr;

	static ptr create(detail::virtual_matrix_store::const_ptr mat) {
		return ptr(new underlying_mat_set(mat));
	}

	mat_id_set get_underlying() const {
		return under_mats;
	}

	size_t get_num_underlying() const {
		return under_mat_set.size();
	}

	bool is_subset_of(const underlying_mat_set &set) const {
		for (size_t i = 0; i < under_mats.size(); i++) {
			auto it = set.under_mat_set.find(under_mats[i]);
			// We can't find the underlying matrix in the other set.
			if (it == set.under_mat_set.end())
				return false;
		}
		return true;
	}

	bool merge(const underlying_mat_set &set) {
		// We don't merge two sets if the underlying matrices aren't
		// the subset of the other.
		if (!set.is_subset_of(*this))
			return false;
		vmats.insert(set.vmats.begin(), set.vmats.end());
		return true;
	}

	void materialize(bool par_access);
};

void underlying_mat_set::materialize(bool par_access)
{
	// TODO we need to deal with the case that some matrices are tall and
	// some are wide. It's better to materialize tall and wide matrices together.
	std::vector<detail::matrix_store::const_ptr> wide_vmats;
	std::vector<detail::matrix_store::const_ptr> tall_vmats;
	for (auto it = vmats.begin(); it != vmats.end(); it++) {
		auto vmat = it->second;
		if (vmat->is_wide()) {
			wide_vmats.push_back(vmat);
			assert(wide_vmats.front()->get_num_cols()
					== wide_vmats.back()->get_num_cols());
		}
		else {
			tall_vmats.push_back(vmat);
			assert(tall_vmats.front()->get_num_rows()
					== tall_vmats.back()->get_num_rows());
		}
	}

	// TODO we might want to materialize with matrices from other underlying set
	// together.
	if (!wide_vmats.empty()) {
		detail::portion_mapply_op::const_ptr materialize_op(
				new materialize_mapply_op(wide_vmats[0]->get_type(),
					true, !par_access));
		__mapply_portion(wide_vmats, materialize_op, matrix_layout_t::L_ROW,
				par_access);
	}
	if (!tall_vmats.empty()) {
		detail::portion_mapply_op::const_ptr materialize_op(
				new materialize_mapply_op(tall_vmats[0]->get_type(),
					false, !par_access));
		__mapply_portion(tall_vmats, materialize_op, matrix_layout_t::L_ROW,
				par_access);
	}
}

class vmat_level
{
	// The two containers have the same matrices. The hashtable is used for
	// merging the matrices in the same level and the vector is used for
	// the rest cases.
	std::unordered_map<mat_id_set, underlying_mat_set::ptr> map;
	std::vector<underlying_mat_set::ptr> vec;

	// Merge the matrix to one of the matrices in the this level.
	bool merge(const underlying_mat_set &set);
	// The level of materialization.
	// It is equal to the number of underlying matrices of a virtual
	// matrix.
	size_t level_id;
public:
	typedef std::shared_ptr<vmat_level> ptr;

	vmat_level(size_t level_id) {
		this->level_id = level_id;
	}

	bool is_empty() const {
		return vec.empty();
	}

	size_t get_num_matrices() const {
		return vec.size();
	}

	size_t get_num_underlying() const {
		return level_id;
	}

	// Add the matrix to this level. It automatically merges the new matrix
	// to an existing one if possible.
	void add(underlying_mat_set::ptr set) {
		assert(set->get_num_underlying() == level_id);
		auto ret = map.insert(
				std::pair<mat_id_set, underlying_mat_set::ptr>(
					set->get_underlying(), set));
		// If the underlying matrices already exist, we merge them.
		if (!ret.second) {
			bool success = ret.first->second->merge(*set);
			assert(success);
		}
		else
			vec.push_back(set);
	}

	// Merge the matrices in this level to upper levels
	// If some of them can't be merged, we create a new one to contain
	// the unmerged ones.
	vmat_level::ptr merge_to(
			std::vector<vmat_level::ptr>::const_iterator ul_begin,
			std::vector<vmat_level::ptr>::const_iterator ul_end) const;

	void materialize(bool par_access) {
		for (size_t i = 0; i < vec.size(); i++)
			vec[i]->materialize(par_access);
	}
};

bool vmat_level::merge(const underlying_mat_set &set)
{
	std::vector<underlying_mat_set::ptr> mergable;
	for (size_t i = 0; i < vec.size(); i++)
		if (set.is_subset_of(*vec[i]))
			mergable.push_back(vec[i]);
	// If we can't merge to any of the matrices.
	if (mergable.empty())
		return false;
	if (mergable.size() == 1)
		mergable[0]->merge(set);
	else
		// We just randomly pick one to merge to.
		mergable[random() % mergable.size()]->merge(set);
	return true;
}

vmat_level::ptr vmat_level::merge_to(
		std::vector<vmat_level::ptr>::const_iterator ul_begin,
		std::vector<vmat_level::ptr>::const_iterator ul_end) const
{
	std::vector<underlying_mat_set::ptr> unmerged;
	for (auto it = map.begin(); it != map.end(); it++) {
		auto set = it->second;
		auto ul_it = ul_begin;
		for (; ul_it != ul_end; ul_it++) {
			assert((*ul_it)->get_num_underlying() > get_num_underlying());
			// If we can merge to a level, we stop here.
			if ((*ul_it)->merge(*set))
				break;
		}
		// If we can't merge to any upper level, we add to the unmerged set.
		if (ul_it == ul_end)
			unmerged.push_back(set);
	}
	if (unmerged.empty())
		return vmat_level::ptr();
	else {
		vmat_level::ptr level(new vmat_level(this->level_id));
		// We don't need to construct the map any more because we won't add
		// another matrix to this level.
		level->vec = unmerged;
		return level;
	}
}

class vmat_levels
{
	std::vector<vmat_level::ptr> underlying_levels;

	void expand(size_t num_levels);
public:
	typedef std::shared_ptr<vmat_levels> ptr;

	void add(underlying_mat_set::ptr set) {
		size_t level_id = set->get_num_underlying();
		// If there isn't a level to store the new matrix, we expend the vector.
		expand(level_id + 1);
		assert(underlying_levels.size() > level_id);
		underlying_levels[level_id]->add(set);
	}

	bool is_empty() const {
		for (size_t i = 0; i < underlying_levels.size(); i++)
			if (!underlying_levels[i]->is_empty())
				return false;
		return true;
	}

	void materialize(bool par_access);
};

void vmat_levels::expand(size_t num_levels)
{
	if (underlying_levels.size() < num_levels) {
		size_t num_orig = underlying_levels.size();
		underlying_levels.resize(num_levels);
		for (size_t i = num_orig; i < num_levels; i++) {
			if (underlying_levels[i] == NULL)
				underlying_levels[i]
					= vmat_level::ptr(new vmat_level(i));
		}
	}
}

void vmat_levels::materialize(bool par_access)
{
	// We merge them first.
	for (auto it = underlying_levels.begin();
			it != underlying_levels.end() - 1; it++) {
		// we merge the current level to one of the upper levels.
		// If some matrices in this level can't be merged, we keep
		// them in this level.
		vmat_level::ptr unmerged = (*it)->merge_to(it + 1,
				underlying_levels.end());
		*it = unmerged;
	}

	for (size_t i = 0; i < underlying_levels.size(); i++)
		if (underlying_levels[i])
			underlying_levels[i]->materialize(par_access);
}

}

bool materialize(std::vector<dense_matrix::ptr> &mats, bool par_access,
		bool mater_self)
{
	if (mats.empty())
		return true;

	for (size_t i = 0; i < mats.size(); i++)
		const_cast<detail::matrix_store &>(mats[i]->get_data()).reset_dag_ref();
	for (size_t i = 0; i < mats.size(); i++)
		const_cast<detail::matrix_store &>(mats[i]->get_data()).inc_dag_ref(
				detail::INVALID_MAT_ID);

	bool ret = true;
	try {
		vmat_levels::ptr levels(new vmat_levels());
		for (size_t i = 0; i < mats.size(); i++) {
			// If this isn't a virtual matrix, skip it.
			if (!mats[i]->is_virtual())
				continue;

			// If the virtual matrix is a TAS matrix, we want it to save
			// the materialized result.
			mats[i]->set_materialize_level(materialize_level::MATER_FULL);

			auto vmats = mats[i]->get_compute_matrices();
			for (size_t j = 0; j < vmats.size(); j++)
				levels->add(underlying_mat_set::create(vmats[j]));
		}
		if (!levels->is_empty()) {
			levels->materialize(par_access);
			// Now all virtual matrices contain the materialized results.
			if (mater_self) {
				for (size_t i = 0; i < mats.size(); i++)
					ret = ret && mats[i]->materialize_self();
			}
		}
	} catch (std::exception &e) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"fail to materialize multiple matrices: %1%") % e.what();
		ret = false;
	}
	for (size_t i = 0; i < mats.size(); i++)
		const_cast<detail::matrix_store &>(mats[i]->get_data()).reset_dag_ref();
	return ret;
}

}
