#ifndef __FM_UTILS_H__
#define __FM_UTILS_H__

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

#include "sparse_matrix_format.h"
#include "data_frame.h"

namespace fm
{

class sparse_matrix;

typedef uint32_t ele_idx_t;
const ele_idx_t INVALID_IDX_VAL = std::numeric_limits<ele_idx_t>::max();

void set_deduplicate(bool v);
void set_remove_self_edge(bool v);

std::shared_ptr<sparse_matrix> create_2d_matrix(data_frame::ptr df,
		const block_2d_size &block_size, const fm::scalar_type *entry_type,
		bool is_sym, const std::string &name = "");

}

#endif
