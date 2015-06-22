#ifndef __TYPE_CAST_H__
#define __TYPE_CAST_H__

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

#include <stdlib.h>

namespace fm
{

class scalar_type;

class type_cast
{
public:
	static bool require_cast(const scalar_type &t1, const scalar_type &t2);
	virtual void cast(size_t num, const void *in, void *out) const = 0;
};

template<class T1, class T2>
class type_cast_impl: public type_cast
{
public:
	virtual void cast(size_t num, const void *in, void *out) const {
		const T1 *t_in = (const T1 *) in;
		T2 *t_out = (T2 *) out;
		for (size_t i = 0; i < num; i++)
			t_out[i] = t_in[i];
	}
};

template<class T1, class T2>
const type_cast &get_type_cast()
{
	static type_cast_impl<T1, T2> cast;
	return cast;
}

}

#endif
