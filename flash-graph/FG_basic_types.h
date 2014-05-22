#ifndef __FG_BASIC_TYPES_H__
#define __FG_BASIC_TYPES_H__

/**
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

typedef unsigned int vsize_t;
typedef unsigned int vertex_id_t;
const vertex_id_t MAX_VERTEX_ID = UINT_MAX;
const vertex_id_t INVALID_VERTEX_ID = -1;
const size_t MAX_VERTEX_SIZE = INT_MAX;

#endif
