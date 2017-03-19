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

#include "col_vec.h"
#include "vector.h"
#include "mem_worker_thread.h"
#include "generic_hashtable.h"

namespace fm
{

col_vec::ptr col_vec::create(detail::matrix_store::ptr store)
{
	if (store->get_num_cols() > 1 && store->get_num_rows() > 1) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't convert a matrix store with multiple cols&rows to a vector";
		assert(0);
		return ptr();
	}
	dense_matrix::ptr mat = dense_matrix::create(store);
	if (mat->get_num_cols() > 1)
		mat = mat->transpose();
	if (mat->get_data().store_layout() == matrix_layout_t::L_ROW)
		mat = mat->conv2(matrix_layout_t::L_COL);
	return ptr(new col_vec(mat->get_raw_store()));
}

col_vec::ptr col_vec::create(dense_matrix::ptr mat)
{
	if (mat->get_num_cols() > 1 && mat->get_num_rows() > 1) {
		printf("the input matrix has %ld rows and %ld cols\n",
				mat->get_num_rows(), mat->get_num_cols());
		BOOST_LOG_TRIVIAL(error)
			<< "can't convert a matrix with multiple cols&rows to a vector";
		assert(0);
		return ptr();
	}
	if (mat->get_num_cols() > 1)
		mat = mat->transpose();

	if (mat->get_data().store_layout() == matrix_layout_t::L_ROW)
		mat = mat->conv2(matrix_layout_t::L_COL);
	return ptr(new col_vec(mat->get_raw_store()));
}

col_vec::ptr col_vec::create(vector::const_ptr vec)
{
	dense_matrix::ptr mat = vec->conv2mat(vec->get_length(), 1, false);
	assert(mat->store_layout() == matrix_layout_t::L_COL);
	if (mat == NULL)
		return col_vec::ptr();
	else
		return col_vec::create(mat);
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

data_frame::ptr col_vec::groupby(agg_operate::const_ptr op, bool with_val,
		bool sorted) const
{
	std::vector<detail::matrix_store::const_ptr> stores(1);
	stores[0] = get_raw_store();
	agg_vec_portion_op *_portion_op = new agg_vec_portion_op(op);
	detail::portion_mapply_op::const_ptr portion_op(_portion_op);
	detail::__mapply_portion(stores, portion_op, matrix_layout_t::L_COL);
	generic_hashtable::ptr agg_res = _portion_op->get_agg();

	// The key-value pairs got from the hashtable aren't in any order.
	// We need to sort them before returning them.
	data_frame::ptr ret = data_frame::create();
	data_frame::ptr df = agg_res->conv2df();
	if (sorted) {
		vector::ptr keys = vector::create(df->get_vec(0));
		data_frame::ptr sorted_keys = keys->sort_with_index();
		detail::smp_vec_store::const_ptr idx_store
			= std::dynamic_pointer_cast<const detail::smp_vec_store>(
					sorted_keys->get_vec("idx"));
		detail::smp_vec_store::const_ptr agg_store
			= std::dynamic_pointer_cast<const detail::smp_vec_store>(df->get_vec(1));
		if (with_val)
			ret->add_vec("val", sorted_keys->get_vec("val"));
		ret->add_vec("agg", agg_store->get(*idx_store));
	}
	else {
		if (with_val)
			ret->add_vec("val", df->get_vec(0));
		ret->add_vec("agg", df->get_vec(1));
	}
	return ret;
}

}
