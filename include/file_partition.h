#ifndef __FILE_PARTITION_H__
#define __FILE_PARTITION_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "file_mapper.h"

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

	off_t map_backwards(int idx, off_t off_in_file) const {
		return mapper->map_backwards(indices[idx], off_in_file);
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

#endif
