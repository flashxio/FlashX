#ifndef __COMPUTE_STAT_H__
#define __COMPUTE_STAT_H__

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

#include "concurrency.h"

/*
 * These provides a thread-safe way to compute some statistics:
 *	mean, max, min, etc on streaming data.
 */

/*
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
