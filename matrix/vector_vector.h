#ifndef __VECTOR_VECTOR_H__
#define __VECTOR_VECTOR_H__

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

#include <memory>

namespace fm
{

class vector;
class mem_vector;
class scalar_type;

/*
 * This stores a vector of vectors. It's similar to the row-wise matrix,
 * but this allows each vector to have different lengths.
 */
class vector_vector
{
public:
	typedef std::shared_ptr<vector_vector> ptr;

	virtual ~vector_vector() {
	}

	virtual size_t get_num_vecs() const = 0;
	virtual size_t get_tot_num_entries() const = 0;
	virtual size_t get_length(off_t idx) const = 0;

	// We can assume each vector can be stored in memory.
	virtual bool append(const mem_vector &vec) = 0;

	/*
	 * Catenate all vectors into a single vector.
	 */
	virtual std::shared_ptr<vector> cat() const = 0;

	virtual const scalar_type &get_type() const = 0;
};

}

#endif
