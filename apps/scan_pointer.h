#ifndef __SCAN_POINTER_H__
#define __SCAN_POINTER_H__

/**
 * Copyright 2014 Da Zheng
 *
 * This file is part of SA-GraphLib.
 *
 * SA-GraphLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SA-GraphLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SA-GraphLib.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <stdlib.h>

/**
 * This class helps to scan a list from the beginning to the end
 * or from the end to the beginning.
 */
struct scan_pointer
{
	size_t range_start;
	size_t size;
	size_t idx;
	bool forward;
public:
	scan_pointer(size_t size, bool forward) {
		this->forward = forward;
		this->range_start = 0;
		this->size = size;
		if (forward)
			idx = range_start;
		else
			idx = size + range_start;
	}

	size_t get_num_remaining() const {
		if (forward)
			return size + range_start - idx;
		else
			return idx - range_start;
	}

	size_t get_curr_loc() const {
		return idx;
	}

	size_t move(size_t dist) {
		if (forward)
			idx += dist;
		else
			idx -= dist;
		return idx;
	}

	size_t get_range_start() const {
		return range_start;
	}

	size_t get_range_end() const {
		return range_start + size;
	}

	size_t get_range_size() const {
		return size;
	}

	void set_scan_dir(bool forward) {
		if (this->forward == forward)
			return;

		size_t new_range_start;
		size_t new_range_size;
		if (forward) {
			// So the original direction is backward.
			new_range_start = get_range_start();
			new_range_size
				= get_curr_loc() - get_range_start();
			idx = new_range_start;
		}
		else {
			// So the original direction is forward.
			new_range_start = get_curr_loc();
			new_range_size
				= get_range_end() - get_curr_loc();
			idx = new_range_start + new_range_size;

		}
		range_start = new_range_start;
		size = new_range_size;
		this->forward = forward;
	}
};

#endif
