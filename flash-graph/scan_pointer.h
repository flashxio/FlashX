#ifndef __SCAN_POINTER_H__
#define __SCAN_POINTER_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashGraph.
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
		if (forward) {
			idx += dist;
			idx = min(idx, size);
		}
		else if (idx >= dist)
			idx -= dist;
		else
			idx = 0;
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
