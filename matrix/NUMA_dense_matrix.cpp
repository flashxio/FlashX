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

#if 0

#include <boost/format.hpp>

#include "NUMA_dense_matrix.h"
#include "mem_worker_thread.h"

namespace fm
{

NUMA_row_tall_dense_matrix::NUMA_row_tall_dense_matrix(size_t nrow, size_t ncol,
		int num_nodes, const scalar_type &type): dense_matrix(nrow, ncol,
			type, true), mapper(num_nodes)
{
	data.resize(num_nodes);
	std::vector<size_t> local_lens = mapper.cal_local_lengths(nrow);
	for (int node_id = 0; node_id < num_nodes; node_id++)
		data[node_id] = detail::raw_data_array(
				local_lens[node_id] * ncol * get_entry_size(), node_id);
}

char *NUMA_row_tall_dense_matrix::get_row(off_t row_idx)
{
	auto phy_loc = mapper.map2physical(row_idx);
	return data[phy_loc.first].get_raw()
		+ phy_loc.second * get_num_cols() * get_entry_size();
}

const char *NUMA_row_tall_dense_matrix::get_row(off_t row_idx) const
{
	auto phy_loc = mapper.map2physical(row_idx);
	return data[phy_loc.first].get_raw()
		+ phy_loc.second * get_num_cols() * get_entry_size();
}

const char *NUMA_row_tall_dense_matrix::get_rows(off_t row_start,
		off_t row_end) const
{
	if (mapper.get_logical_range_id(row_start)
			!= mapper.get_logical_range_id(row_end - 1)) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"[%1%, %2%) isn't in the same range") % row_start % row_end;
		return NULL;
	}
	return get_row(row_start);
}

char *NUMA_row_tall_dense_matrix::get_rows(off_t row_start, off_t row_end)
{
	if (mapper.get_logical_range_id(row_start)
			!= mapper.get_logical_range_id(row_end - 1)) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"[%1%, %2%) isn't in the same range") % row_start % row_end;
		return NULL;
	}
	return (char *) get_row(row_start);
}

namespace
{

class set_row_operate: public detail::set_range_operate
{
	const set_operate &op;
	const detail::NUMA_mapper &mapper;
	size_t num_cols;
	size_t entry_size;
public:
	set_row_operate(const set_operate &_op, const detail::NUMA_mapper &_mapper,
			size_t num_cols, size_t entry_size): op(_op), mapper(_mapper) {
		this->num_cols = num_cols;
		this->entry_size = entry_size;
	}

	void set(char *buf, size_t size, off_t local_off, int node_id) const {
		size_t row_size = num_cols * entry_size;
		assert(local_off % row_size == 0);
		assert(size % row_size == 0);
		size_t num_rows = size / row_size;
		size_t local_row_idx = local_off / row_size;
		size_t global_row_idx = mapper.map2logical(node_id, local_row_idx);
		for (size_t i = 0; i < num_rows; i++)
			op.set(buf + i * row_size, num_cols, global_row_idx + i, 0);
	}
};

class set_col_operate: public set_operate
{
	const set_operate &op;
	size_t col_idx;
public:
	set_col_operate(const set_operate &_op, size_t col_idx): op(_op) {
		this->col_idx = col_idx;
	}

	virtual void set(void *arr, size_t num_eles, off_t row_idx, off_t) const {
		op.set(arr, num_eles, row_idx, this->col_idx);
	}

	virtual const scalar_type &get_type() const {
		return op.get_type();
	}
};

}

void NUMA_row_tall_dense_matrix::set_data(const set_operate &op)
{
	size_t row_size = get_num_cols() * get_entry_size();
	set_row_operate set_row(op, mapper, get_num_cols(), get_entry_size());
	set_array_ranges(mapper, get_num_rows(), row_size, set_row, data);
}

dense_matrix::ptr NUMA_row_tall_dense_matrix::deep_copy() const
{
	NUMA_row_tall_dense_matrix *ret = new NUMA_row_tall_dense_matrix(*this);
	for (size_t i = 0; i < ret->data.size(); i++) {
		// TODO we should do this in parallel.
		ret->data[i] = ret->data[i].deep_copy();
	}
	return dense_matrix::ptr(ret);
}

NUMA_col_tall_dense_matrix::NUMA_col_tall_dense_matrix(size_t nrow,
		size_t ncol, int num_nodes, const scalar_type &type): dense_matrix(
			nrow, ncol, type, true)
{
	data.resize(ncol);
	for (size_t i = 0; i < ncol; i++)
		data[i] = NUMA_vector::create(nrow, num_nodes, type);
}

void NUMA_col_tall_dense_matrix::set_data(const set_operate &op)
{
	for (size_t i = 0; i < data.size(); i++)
		data[i]->set_data(set_col_operate(op, i));
}

dense_matrix::ptr NUMA_col_tall_dense_matrix::deep_copy() const
{
	NUMA_col_tall_dense_matrix *ret = new NUMA_col_tall_dense_matrix(*this);
	for (size_t i = 0; i < ret->data.size(); i++)
		ret->data[i] = NUMA_vector::cast(ret->data[i]->deep_copy());
	return dense_matrix::ptr(ret);
}

dense_matrix::ptr NUMA_col_tall_dense_matrix::inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	// TODO
	assert(0);
}

dense_matrix::ptr NUMA_col_tall_dense_matrix::mapply2(const dense_matrix &m,
		const bulk_operate &op) const
{
	// TODO
	assert(0);
}

}

#endif
