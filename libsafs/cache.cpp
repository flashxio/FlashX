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

#include "cache.h"
#include "io_interface.h"

void page_byte_array::memcpy(off_t rel_off, char buf[], size_t size) const
{
	if (size == 0)
		return;

	// The offset relative to the beginning of the page array.
	off_t off = get_offset_in_first_page() + rel_off;
	off_t end = off + size;
	assert((size_t) end <= get_size() + get_offset_in_first_page());

	off_t page_begin = ROUND(off, PAGE_SIZE);
	// If the element crosses the page boundary.
	if (end - page_begin > PAGE_SIZE) {
		off_t pg_idx = off / PAGE_SIZE;
		off_t last_pg_idx = end / PAGE_SIZE;
		off_t off_in_pg = off % PAGE_SIZE;
		size_t part1_size = PAGE_SIZE - off_in_pg;
		size_t copied = 0;
		::memcpy(buf, ((char *) get_page(pg_idx)->get_data()) + off_in_pg,
				part1_size);
		copied += part1_size;
		pg_idx++;
		while (pg_idx < last_pg_idx) {
			::memcpy(buf + copied, get_page(pg_idx)->get_data(), PAGE_SIZE);
			copied += PAGE_SIZE;
			pg_idx++;
		}
		size_t part2_size = size - copied;
		if (part2_size > 0)
			::memcpy(buf + copied, get_page(pg_idx)->get_data(), part2_size);
	}
	else {
		off_t pg_idx = off / PAGE_SIZE;
		off_t off_in_pg = off % PAGE_SIZE;
		::memcpy(buf, ((char *) get_page(pg_idx)->get_data()) + off_in_pg, size);
	}
}

void page::hit()
{
	int get_file_weight(file_id_t file_id);
	// Not all files are treated equally. We make the pages of the files with
	// higher weight stay in the page cache longer.
	hits += get_file_weight(file_id);
}
