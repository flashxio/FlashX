#ifndef __DATA_IO_H__
#define __DATA_IO_H__

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

/*
 * This file contains the code that imports data from other sources and
 * exports data.
 */

#include <stdlib.h>

#include <memory>
#include <vector>

namespace fm
{

class data_frame;
class scalar_type;

class line_parser
{
public:
	virtual size_t parse(const std::vector<std::string> &lines,
			data_frame &df) const = 0;
	virtual size_t get_num_cols() const = 0;
	virtual const scalar_type &get_col_type(off_t idx) const = 0;
	virtual std::string get_col_name(off_t idx) const = 0;
};

std::shared_ptr<data_frame> read_lines(const std::vector<std::string> &files,
		const line_parser &parser, bool in_mem);

}

#endif
