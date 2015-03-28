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

const char *NUMA_row_tall_dense_matrix::get_row(off_t row_idx) const
{
	auto phy_loc = mapper.map2physical(row_idx);
	return data[phy_loc.first].get_raw()
		+ phy_loc.second * get_num_cols() * get_entry_size();
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

}

void NUMA_row_tall_dense_matrix::set_data(const set_operate &op)
{
	size_t row_size = get_num_cols() * get_entry_size();
	set_row_operate set_row(op, mapper, get_num_cols(), get_entry_size());
	set_array_ranges(mapper, get_num_rows(), row_size, set_row, data);
}

}
