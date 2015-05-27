#ifndef __STL_ALGS_H__
#define __STL_ALGS_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
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

#include <algorithm>

namespace fm
{

class stl_algs
{
public:
	virtual off_t lower_bound(const char *start, const char *end,
			const char *val) const = 0;
};

template<class T>
class stl_algs_impl: public stl_algs
{
public:
	virtual off_t lower_bound(const char *start, const char *end,
			const char *val) const {
		const T *t_start = (const T *) start;
		const T *t_end = (const T *) end;
		T t_val = *(const T *) val;
		const T *ret = std::lower_bound(t_start, t_end, t_val);
		return ret - t_start;
	}
};

}

#endif
