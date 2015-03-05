#ifndef __MEM_DATA_FRAME_H__
#define __MEM_DATA_FRAME_H__

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

#include "data_frame.h"

namespace fm
{

class vector_vector;

class mem_data_frame: public data_frame
{
	mem_data_frame(): data_frame() {
	}

	mem_data_frame(const std::vector<named_vec_t> &named_vecs): data_frame(
			named_vecs) {
	}
public:
	typedef std::shared_ptr<mem_data_frame> ptr;

	static ptr create() {
		return ptr(new mem_data_frame());
	}

	static ptr create(const std::vector<named_vec_t> &named_vecs) {
		return ptr(new mem_data_frame(named_vecs));
	}

	virtual std::shared_ptr<vector_vector> groupby(const std::string &col_name,
			gr_apply_operate<data_frame> &op) const;
	virtual bool sort(const std::string &col_name);
	virtual bool is_sorted(const std::string &col_name) const;
};

}

#endif
