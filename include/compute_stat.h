#ifndef __COMPUTE_STAT_H__
#define __COMPUTE_STAT_H__

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

#include "concurrency.h"

/**
 * These provides a thread-safe way to compute some statistics:
 *	mean, max, min, etc on streaming data.
 */

/**
 * Compute the average value of elements in a vector.
 */
template<class T>
class stat_mean
{
	atomic_number<T> tot;
	atomic_number<T> num;
public:
	void add(T v) {
		tot.inc(v);
		num.inc(1);
	}

	double get() const {
		return ((double) tot.get()) / num.get();
	}
};

template<class T>
class stat_max
{
	atomic_number<T> max;
public:
	void add(T v) {
		while (true) {
			T orig_max = max.get();
			// The original max is smaller, try to switch to the new max.
			if (orig_max < v) {
				bool ret = max.CAS(orig_max, v);
				// If we succeed, we can stop now.
				if (ret)
					break;
				// Otherwise, we need to try again.
			}
			else
				// The original max is larger, we just stop the loop and
				// do nothing.
				break;
		}
	}

	T get() const {
		return max.get();
	}
};

#endif
