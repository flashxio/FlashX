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

#include "generic_type.h"
#include "type_cast.h"

namespace fm
{

bool type_cast::require_cast(const scalar_type &t1, const scalar_type &t2)
{
	if (t1 == t2)
		return false;
	// If the two types require different memory storage size, we definitely
	// need to cast them.
	if (t1.get_size() != t2.get_size())
		return true;

	if ((t1.get_type() == prim_type::P_SHORT
				&& t2.get_type() == prim_type::P_USHORT)
			|| (t1.get_type() == prim_type::P_INTEGER
				&& t2.get_type() == prim_type::P_UINT)
			|| (t1.get_type() == prim_type::P_LONG
				&& t2.get_type() == prim_type::P_ULONG))
		return false;

	if ((t2.get_type() == prim_type::P_SHORT
				&& t1.get_type() == prim_type::P_USHORT)
			|| (t2.get_type() == prim_type::P_INTEGER
				&& t1.get_type() == prim_type::P_UINT)
			|| (t2.get_type() == prim_type::P_LONG
				&& t1.get_type() == prim_type::P_ULONG))
		return false;

	// We need to cast the rest of the type pairs.
	return true;
}

}
