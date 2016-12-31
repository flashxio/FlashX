#ifndef __SAFS_HEADER_H__
#define __SAFS_HEADER_H__

/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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

#include "io_request.h"

namespace safs
{

class safs_header
{
	static const int64_t MAGIC_NUMBER = 0x123456789FFFFFEL;
	static const int CURR_VERSION = 1;

	int64_t magic_number;
	int version_number;
	// In the number of pages.
	uint32_t block_size;
	uint32_t mapping_option;
	uint32_t writable;
	uint64_t num_bytes;
public:
	static size_t get_header_size() {
		return PAGE_SIZE;
	}

	safs_header() {
		this->magic_number = MAGIC_NUMBER;
		this->version_number = CURR_VERSION;
		this->block_size = 0;
		this->mapping_option = 0;
		this->writable = false;
		this->num_bytes = 0;
	}

	safs_header(int block_size, int mapping_option, bool writable,
			size_t file_size) {
		this->magic_number = MAGIC_NUMBER;
		this->version_number = CURR_VERSION;
		this->block_size = block_size;
		this->mapping_option = mapping_option;
		this->writable = writable;
		this->num_bytes = file_size;
	}

	int get_block_size() const {
		return block_size;
	}

	int get_mapping_option() const {
		return mapping_option;
	}

	bool is_writable() const {
		return writable;
	}

	bool is_safs_file() const {
		return magic_number == MAGIC_NUMBER;
	}

	bool is_right_version() const {
		return version_number == CURR_VERSION;
	}

	bool is_valid() const {
		return block_size != 0;
	}

	void resize(size_t new_size) {
		this->num_bytes = new_size;
	}

	size_t get_size() const {
		return num_bytes;
	}
};

}

#endif
