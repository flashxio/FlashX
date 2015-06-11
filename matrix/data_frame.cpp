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

#include "data_frame.h"
#include "bulk_operate.h"
#include "mem_vec_store.h"
#include "EM_vector.h"
#include "mem_vv_store.h"
#include "mem_vector_vector.h"

namespace fm
{

data_frame::data_frame(const std::vector<named_vec_t> &named_vecs)
{
	assert(!named_vecs.empty());
	this->named_vecs = named_vecs;
	bool in_mem = named_vecs.front().second->is_in_mem();
	for (auto it = named_vecs.begin(); it != named_vecs.end(); it++) {
		assert(in_mem == it->second->is_in_mem());
		vec_map.insert(*it);
	}
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
	if (sorted_col->is_in_mem()) {
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
		vecs.push_back(detail::EM_vec_store::cast(sorted_col));
		for (size_t i = 0; i < named_vecs.size(); i++) {
			if (named_vecs[i].second != sorted_col) {
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
		if (vec->is_in_mem() != named_vecs.front().second->is_in_mem()) {
			BOOST_LOG_TRIVIAL(error)
				<< "Add a vector in different storage from the data frame.";
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

namespace
{

void expose_portion(const std::vector<detail::mem_vec_store::const_ptr> &sorted_df,
		off_t loc, size_t length, sub_data_frame &sub_df)
{
	// TODO This is an very inefficient implementation.
	sub_df.resize(sorted_df.size());
	for (size_t i = 0; i < sub_df.size(); i++)
		sub_df[i] = sorted_df[i]->get_portion(loc, length);
}

class local_groupby_task: public thread_task
{
	detail::mem_vec_store::const_ptr sub_sorted_col;
	std::vector<detail::mem_vec_store::const_ptr> sub_dfs;
	const gr_apply_operate<sub_data_frame> &op;
	std::vector<detail::mem_vv_store::ptr> &sub_results;
	off_t idx;
public:
	local_groupby_task(detail::mem_vec_store::const_ptr sub_sorted_col,
			const std::vector<detail::mem_vec_store::const_ptr> &sub_dfs,
			const gr_apply_operate<sub_data_frame> &_op,
			std::vector<detail::mem_vv_store::ptr> &_sub_results,
			off_t idx): op(_op), sub_results(_sub_results) {
		this->sub_sorted_col = sub_sorted_col;
		this->sub_dfs = sub_dfs;
		this->idx = idx;
	}

	void run();
};

void local_groupby_task::run()
{
	sub_data_frame sub_df;
	size_t out_size;
	// If the user can predict the number of output elements, we can create
	// a buffer of the expected size.
	if (op.get_num_out_eles() > 0)
		out_size = op.get_num_out_eles();
	else
		// If the user can't, we create a small buffer.
		out_size = 16;
	local_buf_vec_store row(0, out_size, op.get_output_type(), -1);

	detail::mem_vv_store::ptr ret = detail::mem_vv_store::create(
			op.get_output_type());
	sub_results[idx] = ret;
	const agg_operate &find_next
		= sub_sorted_col->get_type().get_agg_ops().get_find_next();
	size_t loc = 0;
	size_t col_len = sub_sorted_col->get_length();
	const char *start = sub_sorted_col->get_raw_arr();
	size_t entry_size = sub_sorted_col->get_entry_size();
	while (loc < col_len) {
		size_t curr_length = col_len - loc;
		const char *curr_ptr = start + entry_size * loc;
		size_t rel_end;
		find_next.run(curr_length, curr_ptr, &rel_end);
		// This expose a portion of the data frame.
		expose_portion(sub_dfs, loc, rel_end, sub_df);
		// The first argument is the key and the second one is the value
		// (a data frame)
		op.run(curr_ptr, sub_df, row);
		if (row.get_length() > 0)
			ret->append(row);
		loc += rel_end;
	}
}

}

std::vector<off_t> partition_vector(const detail::mem_vec_store &sorted_vec,
		int num_parts);

vector_vector::ptr data_frame::groupby(const std::string &col_name,
		gr_apply_operate<sub_data_frame> &op) const
{
	data_frame::const_ptr sorted_df = sort(col_name);
	detail::mem_vec_store::const_ptr sorted_col
		= detail::mem_vec_store::cast(sorted_df->get_vec(col_name));

	// We need to find the start location for each thread.
	// The start location is where the value in the sorted array
	// first appears.
	// TODO this only works for vectors stored contiguously in memory.
	// It doesn't work for NUMA vector.
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	int num_threads = mem_threads->get_num_threads();
	std::vector<off_t> par_starts = partition_vector(*sorted_col, num_threads);

	// It's possible that two partitions end up having the same start location
	// because the vector is small or a partition has only one value.
	assert(std::is_sorted(par_starts.begin(), par_starts.end()));
	auto end_par_starts = std::unique(par_starts.begin(), par_starts.end());
	int num_parts = end_par_starts - par_starts.begin() - 1;
	std::vector<detail::mem_vv_store::ptr> sub_results(num_parts);
	for (int i = 0; i < num_parts; i++) {
		off_t start = par_starts[i];
		off_t end = par_starts[i + 1];
		std::vector<detail::mem_vec_store::const_ptr> sub_dfs(
				sorted_df->get_num_vecs());
		for (size_t i = 0; i < sorted_df->get_num_vecs(); i++) {
			detail::smp_vec_store::const_ptr vec = detail::smp_vec_store::cast(
					sorted_df->get_vec(i));
			detail::smp_vec_store::ptr sub_vec
				= detail::smp_vec_store::cast(const_cast<detail::smp_vec_store *>(
							vec.get())->shallow_copy());
			bool ret = sub_vec->expose_sub_vec(start, end - start);
			assert(ret);
			sub_dfs[i] = sub_vec;
		}
		detail::smp_vec_store::ptr sub_sorted_col = detail::smp_vec_store::cast(
				const_cast<detail::mem_vec_store *>(sorted_col.get())->shallow_copy());
		sub_sorted_col->expose_sub_vec(start, end - start);

		// It's difficult to localize computation.
		// TODO can we localize computation?
		int node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id, new local_groupby_task(
					sub_sorted_col, sub_dfs, op, sub_results, i));
	}
	mem_threads->wait4complete();

	if (num_parts == 1)
		return mem_vector_vector::create(sub_results[0]);
	else {
		detail::mem_vv_store::ptr res_vv = sub_results[0];
		std::vector<detail::vec_store::const_ptr> remain(
				sub_results.begin() + 1, sub_results.end());
		bool ret = res_vv->append(remain.begin(), remain.end());
		assert(ret);
		return mem_vector_vector::create(res_vv);
	}
}

}
