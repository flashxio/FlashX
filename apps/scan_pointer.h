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
	size_t size;
	size_t idx;
	bool forward;
public:
	scan_pointer(size_t size, bool forward) {
		this->forward = forward;
		this->size = size;
		if (forward)
			idx = 0;
		else
			idx = size;
	}

	size_t get_num_remaining() const {
		if (forward)
			return size - idx;
		else
			return idx;
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
};

#endif
