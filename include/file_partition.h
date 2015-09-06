#ifndef __FILE_PARTITION_H__
#define __FILE_PARTITION_H__

/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
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

#include "file_mapper.h"

namespace safs
{

/**
 * This class represents a partition in a logical file.
 * A logical file is a data collection viewed by applications. It is
 * divided into partitions and stored on multiple disks physically.
 * When an application needs to access data, it only needs to give locations
 * in the logical file.
 */
class logical_file_partition
{
	file_mapper *mapper;
	std::vector<int> indices;
	// Map the file index in the mapper to the file index in the partition.
	std::vector<int> file_map;
public:
	logical_file_partition(std::vector<int> indices,
			file_mapper *mapper): file_map(mapper->get_num_files(), -1) {
		this->mapper = mapper;
		this->indices = indices;
		for (size_t i = 0; i < indices.size(); i++)
			file_map[indices[i]] = i;
	}

	logical_file_partition(std::vector<int> indices) {
		this->indices = indices;
		mapper = NULL;
	}

	logical_file_partition *create_file_partition(file_mapper *mapper) const {
		assert(mapper);
		return new logical_file_partition(indices, mapper);
	}

	int get_file_id() const {
		assert(mapper);
		return mapper->get_file_id();
	}

	bool is_active() const {
		return mapper != NULL;
	}

	const file_mapper *get_mapper() const {
		assert(mapper);
		return mapper;
	}

	int get_num_files() const {
		return indices.size();
	}

	const std::vector<int> &get_phy_file_indices() const {
		return indices;
	}

	const std::string &get_file_name(int idx) const {
		assert(mapper);
		return mapper->get_file_name(indices[idx]);
	}

	void map(off_t pg_off, block_identifier &bid) const {
		assert(mapper);
		mapper->map(pg_off, bid);
		// We have to make sure the offset does exist in the partition.
		assert(file_map[bid.idx] >= 0);
		bid.idx = file_map[bid.idx];
	}

	int map2file(off_t pg_off) const {
		assert(mapper);
		int idx = mapper->map2file(pg_off);
		assert(file_map[idx] >= 0);
		return file_map[idx];
	}

	int get_partition_id() const {
		return indices[0];
	}
};

}

#endif
