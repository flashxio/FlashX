/*
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
#include "RAID_config.h"

namespace safs
{

file_mapper::ptr file_mapper::create(const safs_header &header,
		const std::vector<part_file_info> &files, const std::string &file_name)
{
	if (!header.is_valid())
		return file_mapper::ptr();

	int block_size = header.get_block_size();
	int mapping_option = header.get_mapping_option();
	switch (mapping_option) {
		case RAID0:
			return file_mapper::ptr(new RAID0_mapper(file_name, files,
						block_size));
		case RAID5:
			return file_mapper::ptr(new RAID5_mapper(file_name, files,
						block_size));
		case HASH:
			return file_mapper::ptr(new hash_mapper(file_name, files,
						block_size));
		default:
			fprintf(stderr, "wrong RAID mapping option\n");
			return file_mapper::ptr();
	}
}

int gen_RAID_rand_start(int num_files)
{
	int r = random();
	while (r == 0) {
		r = random();
	}
	return r % num_files;
}

int RAID0_mapper::rand_start;
int RAID5_mapper::rand_start;

atomic_integer file_mapper::file_id_gen;

std::vector<size_t> file_mapper::get_size_per_disk(size_t size) const
{
	std::vector<size_t> ret(get_num_files());
	for (size_t i = 0; i < ret.size() * 2; i++) {
		if (size < i * STRIPE_BLOCK_SIZE)
			break;
		// Get the offset of the last page in a block.
		off_t off = (size / STRIPE_BLOCK_SIZE - i) * STRIPE_BLOCK_SIZE
			+ STRIPE_BLOCK_SIZE - 1;
		struct block_identifier bid;
		map(off, bid);
		ret[bid.idx] = std::max((size_t) bid.off + 1, ret[bid.idx]);
	}
	return ret;
}

std::vector<size_t> hash_mapper::get_size_per_disk(size_t size) const
{
	std::vector<size_t> ret(get_num_files());
	// We should try the entire cycle the figure out the sizes of all files.
	for (int i = 0; i < CONST_P; i++) {
		if (size < (size_t) i * STRIPE_BLOCK_SIZE)
			break;
		off_t off = (size / STRIPE_BLOCK_SIZE - i) * STRIPE_BLOCK_SIZE
			+ STRIPE_BLOCK_SIZE - 1;
		struct block_identifier bid;
		map(off, bid);
		ret[bid.idx] = std::max((size_t) bid.off + 1, ret[bid.idx]);
	}
	return ret;
}

}
