/**
 * Copyright 2014 Da Zheng
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

#include "cache.h"

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
