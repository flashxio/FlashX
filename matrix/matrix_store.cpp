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

#include "matrix_store.h"
#include "local_matrix_store.h"

namespace fm
{

namespace detail
{

size_t matrix_store::get_num_portions() const
{
	std::pair<size_t, size_t> chunk_size = get_portion_size();
	if (is_wide())
		return ceil(((double) get_num_cols()) / chunk_size.second);
	else
		return ceil(((double) get_num_rows()) / chunk_size.first);
}

}

}
