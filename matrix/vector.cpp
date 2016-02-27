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

#include "common.h"
#include "thread.h"

#include "vector.h"
#include "data_frame.h"
#include "local_vec_store.h"
#include "mem_worker_thread.h"
#include "dense_matrix.h"
#include "bulk_operate_ext.h"
#include "generic_hashtable.h"
#include "local_matrix_store.h"
#include "dense_matrix.h"

namespace fm
{

vector::ptr vector::create(size_t length, const scalar_type &type,
		int num_nodes, bool in_mem, const set_vec_operate &op)
{
	detail::vec_store::ptr vec = detail::vec_store::create(length, type,
			num_nodes, in_mem);
	vec->set_data(op);
	return ptr(new vector(vec));
}

bool vector::verify_groupby(
		const gr_apply_operate<local_vec_store> &op) const
{
	if (op.get_key_type() != get_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "the operator's key type is incompatible with the vector";
		return false;
	}

	return true;
}

namespace
{

class local_groupby_task: public thread_task
{
	local_vec_store::const_ptr sub_vec;
	const gr_apply_operate<local_vec_store> &op;
	std::vector<data_frame::ptr> &sub_results;
	bool with_val;
	off_t idx;
public:
	local_groupby_task(local_vec_store::const_ptr sub_vec,
			const gr_apply_operate<local_vec_store> &_op,
			bool with_val, std::vector<data_frame::ptr> &_sub_results,
			off_t idx): op(_op), sub_results(_sub_results) {
		this->sub_vec = sub_vec;
		this->with_val = with_val;
		this->idx = idx;
	}

	void run() {
		sub_results[idx] = sub_vec->groupby(op, with_val);
	}
};

}

/*
 * Partition the vector into multiple partitions of relatively equal size.
 * It is guaranteed that the same values are all split into the same partition.
 * It returns a vector with `num_parts' + 1 elements. Each element except
 * the last indicates the start location of a partition in the vector.
 * The last element indicates the end of the vector.
 */
std::vector<off_t> partition_vector(const detail::mem_vec_store &sorted_vec,
		int num_parts)
{
	agg_operate::const_ptr find_next
		= sorted_vec.get_type().get_agg_ops().get_find_next();
	std::vector<off_t> par_starts(num_parts + 1);
	for (int i = 0; i < num_parts; i++) {
		off_t start = sorted_vec.get_length() / num_parts * i;
		// This returns the relative start location of the next value.
		find_next->runAgg(sorted_vec.get_length() - start,
				sorted_vec.get_raw_arr() + sorted_vec.get_entry_size() * start,
				&par_starts[i]);
		// This is the absolute start location of this partition.
		par_starts[i] += start;
	}
	par_starts[0] = 0;
	par_starts[num_parts] = sorted_vec.get_length();
	return par_starts;
}

data_frame::ptr vector::groupby(
		const gr_apply_operate<local_vec_store> &op, bool with_val) const
{
	if (!verify_groupby(op))
		return data_frame::ptr();

	// If the vector hasn't been sorted, we need to sort it.
	detail::mem_vec_store::const_ptr sorted_vec;
	if (!this->is_sorted()) {
		// We don't want groupby changes the original vector.
		detail::vec_store::ptr tmp = get_data().deep_copy();
		tmp->sort();
		sorted_vec = detail::mem_vec_store::cast(tmp);
	}
	else
		sorted_vec = detail::mem_vec_store::cast(this->get_raw_store());

	// We need to find the start location for each thread.
	// The start location is where the value in the sorted array
	// first appears.
	// TODO this only works for vectors stored contiguously in memory.
	// It doesn't work for NUMA vector.
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	int num_threads = mem_threads->get_num_threads();
	std::vector<off_t> par_starts = partition_vector(*sorted_vec, num_threads);

	// It's possible that two partitions end up having the same start location
	// because the vector is small or a partition has only one value.
	assert(std::is_sorted(par_starts.begin(), par_starts.end()));
	auto end_par_starts = std::unique(par_starts.begin(), par_starts.end());
	int num_parts = end_par_starts - par_starts.begin() - 1;
	std::vector<data_frame::ptr> sub_results(num_parts);
	for (int i = 0; i < num_parts; i++) {
		off_t start = par_starts[i];
		off_t end = par_starts[i + 1];
		local_vec_store::const_ptr sub_vec = sorted_vec->get_portion(
				start, end - start);

		int node_id = sub_vec->get_node_id();
		if (node_id < 0)
			node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id, new local_groupby_task(sub_vec, op,
					with_val, sub_results, i));
	}
	mem_threads->wait4complete();

	if (num_parts == 1)
		return sub_results[0];
	else {
		std::vector<data_frame::const_ptr> const_sub_results(
				sub_results.begin(), sub_results.end());
		return merge_data_frame(const_sub_results, true);
	}
}

namespace
{

class agg_vec_portion_op: public detail::portion_mapply_op
{
	agg_operate::const_ptr find_next;
	agg_operate::const_ptr agg_op;
	std::vector<generic_hashtable::ptr> tables;
	std::vector<local_vec_store::ptr> lvec_bufs;
	std::vector<local_vec_store::ptr> lkey_bufs;
	std::vector<local_vec_store::ptr> lagg_bufs;
public:
	agg_vec_portion_op(agg_operate::const_ptr agg_op): detail::portion_mapply_op(
			0, 0, agg_op->get_output_type()) {
		this->find_next = agg_op->get_input_type().get_agg_ops().get_find_next();
		this->agg_op = agg_op;
		size_t nthreads = detail::mem_thread_pool::get_global_num_threads();
		tables.resize(nthreads);
		lvec_bufs.resize(nthreads);
		lkey_bufs.resize(nthreads);
		lagg_bufs.resize(nthreads);
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		return detail::portion_mapply_op::const_ptr();
	}

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return "";
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const;

	generic_hashtable::ptr get_agg() const;
};

void agg_vec_portion_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins) const
{
	int thread_id = detail::mem_thread_pool::get_curr_thread_id();
	agg_vec_portion_op *mutable_this = const_cast<agg_vec_portion_op *>(this);
	generic_hashtable::ptr ltable;
	local_vec_store::ptr lvec;
	local_vec_store::ptr lkeys;
	local_vec_store::ptr laggs;
	if (tables[thread_id] == NULL) {
		mutable_this->tables[thread_id] = ins[0]->get_type().create_hashtable(
				agg_op->get_output_type());
		mutable_this->lvec_bufs[thread_id] = local_vec_store::ptr(
				new local_buf_vec_store(0, ins[0]->get_num_rows(),
					ins[0]->get_type(), -1));
		mutable_this->lkey_bufs[thread_id] = local_vec_store::ptr(
				new local_buf_vec_store(0, ins[0]->get_num_rows(),
					ins[0]->get_type(), -1));
		mutable_this->lagg_bufs[thread_id] = local_vec_store::ptr(
				new local_buf_vec_store(0, ins[0]->get_num_rows(),
					agg_op->get_output_type(), -1));
	}
	ltable = tables[thread_id];
	lvec = lvec_bufs[thread_id];
	lkeys = lkey_bufs[thread_id];
	laggs = lagg_bufs[thread_id];
	assert(ltable);
	assert(lvec);
	assert(lkeys);
	assert(laggs);
	size_t llength = ins[0]->get_num_rows();
	assert(ins[0]->get_raw_arr());
	assert(lvec->get_length() >= llength);
	memcpy(lvec->get_raw_arr(), ins[0]->get_raw_arr(),
			lvec->get_entry_size() * llength);
	lvec->get_type().get_sorter().serial_sort(lvec->get_raw_arr(), llength,
			false);

	// Start to aggregate on the values.
	const char *arr = lvec->get_raw_arr();
	size_t key_idx = 0;
	while (llength > 0) {
		size_t num_same = 0;
		find_next->runAgg(llength, arr, &num_same);
		assert(num_same <= llength);
		memcpy(lkeys->get(key_idx), arr, lvec->get_entry_size());
		agg_op->runAgg(num_same, arr, laggs->get(key_idx));

		key_idx++;
		arr += num_same * lvec->get_entry_size();
		llength -= num_same;
	}
	ltable->insert(key_idx, lkeys->get_raw_arr(), laggs->get_raw_arr(), *agg_op);
}

generic_hashtable::ptr agg_vec_portion_op::get_agg() const
{
	generic_hashtable::ptr ret;
	size_t i;
	// Find a local table.
	for (i = 0; i < tables.size(); i++) {
		if (tables[i]) {
			ret = tables[i];
			break;
		}
	}
	// We need to move to the next local table.
	i++;
	// Merge with other local tables if they exist.
	for (; i < tables.size(); i++)
		if (tables[i])
			ret->merge(*tables[i], *agg_op);
	return ret;
}

}

data_frame::ptr vector::groupby(agg_operate::const_ptr op, bool with_val) const
{
	dense_matrix::ptr mat = conv2mat(get_length(), 1, true);
	if (mat == NULL)
		return data_frame::ptr();

	std::vector<detail::matrix_store::const_ptr> stores(1);
	stores[0] = mat->get_raw_store();
	agg_vec_portion_op *_portion_op = new agg_vec_portion_op(op);
	detail::portion_mapply_op::const_ptr portion_op(_portion_op);
	detail::__mapply_portion(stores, portion_op, matrix_layout_t::L_COL);
	generic_hashtable::ptr agg_res = _portion_op->get_agg();

	// The key-value pairs got from the hashtable aren't in any order.
	// We need to sort them before returning them.
	data_frame::ptr df = agg_res->conv2df();
	vector::ptr keys = vector::create(df->get_vec(0));
	data_frame::ptr sorted_keys = keys->sort_with_index();
	detail::smp_vec_store::const_ptr idx_store
		= std::dynamic_pointer_cast<const detail::smp_vec_store>(
				sorted_keys->get_vec("idx"));
	detail::smp_vec_store::const_ptr agg_store
		= std::dynamic_pointer_cast<const detail::smp_vec_store>(df->get_vec(1));
	data_frame::ptr ret = data_frame::create();
	if (with_val)
		ret->add_vec("val", sorted_keys->get_vec("val"));
	ret->add_vec("agg", agg_store->get(*idx_store));
	return ret;
}

bool vector::equals(const vector &vec) const
{
	if (vec.get_length() != this->get_length())
		return false;
	else if (vec.get_type() != this->get_type())
		return false;
	else {
		assert(is_in_mem());
		return memcmp(
				dynamic_cast<const detail::mem_vec_store &>(get_data()).get_raw_arr(),
				dynamic_cast<const detail::mem_vec_store &>(vec.get_data()).get_raw_arr(),
				get_length() * get_entry_size()) == 0;
	}
}

namespace
{

class local_agg_task: public thread_task
{
	local_vec_store::const_ptr sub_vec;
	const bulk_operate &op;
	char *agg_res;
public:
	local_agg_task(local_vec_store::const_ptr sub_vec,
			const bulk_operate &_op, char *agg_res): op(_op) {
		this->sub_vec = sub_vec;
		this->agg_res = agg_res;
	}

	void run() {
		op.runAgg(sub_vec->get_length(), sub_vec->get_raw_arr(), agg_res);
	}
};

}

scalar_variable::ptr vector::aggregate(const bulk_operate &op) const
{
	scalar_variable::ptr res = op.get_output_type().create_scalar();
	size_t num_portions = get_data().get_num_portions();
	size_t portion_size = get_data().get_portion_size();
	// TODO we don't need to allocate a slot for every portion.
	std::unique_ptr<char[]> raw_res(new char[res->get_size() * num_portions]);

	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	for (size_t i = 0; i < num_portions; i++) {
		off_t start = i * portion_size;
		size_t length = std::min(portion_size, get_length() - start);
		local_vec_store::const_ptr sub_vec
			= static_cast<const detail::mem_vec_store &>(
					get_data()).get_portion(start, length);

		int node_id = sub_vec->get_node_id();
		if (node_id < 0)
			node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id, new local_agg_task(sub_vec, op,
					raw_res.get() + i * res->get_size()));
	}
	mem_threads->wait4complete();
	char final_res[res->get_size()];
	op.runAgg(num_portions, raw_res.get(), final_res);
	res->set_raw(final_res, res->get_size());
	return res;
}

scalar_variable::ptr vector::dot_prod(const vector &vec) const
{
	if (get_type() != vec.get_type()) {
		BOOST_LOG_TRIVIAL(error) << "The type isn't compatible";
		return scalar_variable::ptr();
	}

	if (get_length() != vec.get_length()) {
		BOOST_LOG_TRIVIAL(error) << "The two vectors have different lengths";
		return scalar_variable::ptr();
	}

#if 0
	const bulk_operate &multiply = get_type().get_basic_ops().get_multiply();
	const bulk_operate &add = get_type().get_basic_ops().get_add();
	std::vector<size_t> local_lens = mapper.cal_local_lengths(get_length());
	std::unique_ptr<char[]> agg_buf(new char[data.size() * get_entry_size()]);
	// TODO this is a very inefficient implementation.
	for (size_t i = 0; i < data.size(); i++) {
		std::unique_ptr<char[]> multiply_buf(
				new char[local_lens[i] * get_entry_size()]);
		multiply.runAA(local_lens[i], data[i].get_raw(), vec.data[i].get_raw(),
				multiply_buf.get());
		add.runA(local_lens[i], multiply_buf.get(),
				agg_buf.get() + get_entry_size() * i);
	}
	char final_agg[get_entry_size()];
	add.runA(data.size(), agg_buf.get(), final_agg);
#endif

	scalar_variable::ptr res = get_type().create_scalar();
#if 0
	res->set_raw(final_agg, get_entry_size());
#endif

	return res;
}

vector::ptr vector::sort() const
{
	// TODO it's unnecessary to copy first and sort.
	detail::vec_store::ptr vec = get_data().deep_copy();
	vec->sort();
	return vector::create(vec);
}

data_frame::ptr vector::sort_with_index() const
{
	detail::vec_store::ptr vec = get_data().deep_copy();
	detail::vec_store::ptr idx = vec->sort_with_index();
	data_frame::ptr df = data_frame::create();
	df->add_vec("idx", idx);
	df->add_vec("val", vec);
	return df;
}

dense_matrix::ptr vector::conv2mat(size_t nrow, size_t ncol,
		bool byrow) const
{
	detail::matrix_store::const_ptr mat = get_data().conv2mat(nrow, ncol, byrow);
	if (mat == NULL)
		return dense_matrix::ptr();
	return dense_matrix::create(mat);
}

bool vector::export2(FILE *f) const
{
	if (!is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "Doesn't support to write an EM vector to a file";
		return false;
	}

	size_t ret = fwrite(
			dynamic_cast<const detail::mem_vec_store &>(get_data()).get_raw_arr(),
			get_length() * get_type().get_size(), 1, f);
	if (ret == 0) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"can't write, error: %1%") % strerror(errno);
		return false;
	}
	return true;
}

vector::ptr vector::sapply(bulk_uoperate::const_ptr op) const
{
	dense_matrix::ptr tmp = conv2mat(get_length(), 1, false);
	if (tmp == NULL)
		return vector::ptr();
	dense_matrix::ptr tmp1 = tmp->sapply(op);
	return tmp1->get_col(0);
}

vector::ptr vector::get(const std::vector<off_t> &idxs) const
{
	detail::matrix_store::const_ptr mat = get_data().conv2mat(get_length(),
			1, false);
	if (mat == NULL)
		return vector::ptr();
	mat = mat->get_rows(idxs);
	if (mat == NULL)
		return vector::ptr();
	detail::vec_store::const_ptr vec = mat->get_col_vec(0);
	if (vec)
		return vector::create(vec);
	else
		return vector::ptr();
}

}
