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

#include "common.h"
#include "thread.h"

#include "mem_vector.h"
#include "mem_data_frame.h"
#include "mem_vector_vector.h"
#include "local_vec_store.h"
#include "mem_worker_thread.h"

namespace fm
{

mem_vector::ptr mem_vector::create(size_t length, const scalar_type &type,
		const set_operate &op)
{
	// TODO I should allow users to create a NUMA vector store as well.
	detail::mem_vec_store::ptr vec = detail::smp_vec_store::create(length, type);
	vec->set_data(op);
	return ptr(new mem_vector(vec));
}

mem_vector::ptr mem_vector::cast(vector::ptr vec)
{
	if (!vec->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't cast a non-in-mem vector to in-mem vector";
		return mem_vector::ptr();
	}
	return std::static_pointer_cast<mem_vector>(vec);
}

bool mem_vector::verify_groupby(
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

data_frame::ptr mem_vector::groupby(
		const gr_apply_operate<local_vec_store> &op, bool with_val) const
{
	const agg_operate &find_next = get_type().get_agg_ops().get_find_next();
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
	off_t par_starts[num_threads + 1];
	for (int i = 0; i < num_threads; i++) {
		off_t start = sorted_vec->get_length() / num_threads * i;
		// This returns the relative start location of the next value.
		find_next.run(sorted_vec->get_length() - start,
				sorted_vec->get_raw_arr() + get_entry_size() * start,
				par_starts + i);
		// This is the absolute start location of this partition.
		par_starts[i] += start;
	}
	par_starts[0] = 0;
	par_starts[num_threads] = sorted_vec->get_length();

	// It's possible that two partitions end up having the same start location
	// because the vector is small or a partition has only one value.
	assert(std::is_sorted(par_starts, par_starts + num_threads + 1));
	off_t *end_par_starts = std::unique(par_starts, par_starts + num_threads + 1);
	int num_parts = end_par_starts - par_starts - 1;
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
		return merge_data_frame(const_sub_results);
	}
}

bool mem_vector::equals(const mem_vector &vec) const
{
	if (vec.get_length() != this->get_length())
		return false;
	else if (vec.get_type() != this->get_type())
		return false;
	else
		return memcmp(this->get_raw_arr(), vec.get_raw_arr(),
				get_length() * get_entry_size()) == 0;
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
		op.runA(sub_vec->get_length(), sub_vec->get_raw_arr(), agg_res);
	}
};

}

scalar_variable::ptr mem_vector::aggregate(const bulk_operate &op) const
{
	scalar_variable::ptr res = op.get_output_type().create_scalar();
	size_t num_portions = get_data().get_num_portions();
	size_t portion_size = get_data().get_portion_size();
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
	op.runA(num_portions, raw_res.get(), final_res);
	res->set_raw(final_res, res->get_size());
	return res;
}

scalar_variable::ptr mem_vector::dot_prod(const vector &vec) const
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

vector::ptr mem_vector::sort() const
{
	detail::mem_vec_store::ptr vec = detail::mem_vec_store::cast(
			get_data().deep_copy());
	vec->sort();
	return mem_vector::create(vec);
}

data_frame::ptr mem_vector::sort_with_index() const
{
	detail::vec_store::ptr vec = get_data().deep_copy();
	detail::vec_store::ptr idx = vec->sort_with_index();
	mem_data_frame::ptr df = mem_data_frame::create();
	df->add_vec("idx", idx);
	df->add_vec("val", vec);
	return df;
}

bool mem_vector::export2(FILE *f) const
{
	size_t ret = fwrite(get_raw_arr(),
			get_length() * get_type().get_size(), 1, f);
	if (ret == 0) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"can't write, error: %1%") % strerror(errno);
		return false;
	}
	return true;
}

}
