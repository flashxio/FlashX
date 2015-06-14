#ifndef __VV_STORE_H__
#define __VV_STORE_H__

/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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

#include "vec_store.h"

namespace fm
{

namespace detail
{

class vv_store: public vec_store
{
	// The offsets (in #bytes) of the vectors in the data array.
	// The last offset is the end of the last vector.
	std::vector<off_t> vec_offs;
	vec_store::ptr store;
protected:
	const vec_store &get_data() const {
		return *store;
	}
	vec_store &get_data() {
		return *store;
	}
	off_t get_vec_off(off_t idx) const {
		return vec_offs[idx];
	}
public:
	typedef std::shared_ptr<vv_store> ptr;
	typedef std::shared_ptr<const vv_store> const_ptr;

	static bool is_vector_vector(const vec_store &vec) {
		return vec.get_entry_size() == 0;
	}

	static ptr cast(vec_store::ptr vec) {
		assert(is_vector_vector(*vec));
		return std::static_pointer_cast<vv_store>(vec);
	}

	vv_store(const scalar_type &type, bool in_mem);
	vv_store(const std::vector<off_t> &offs, vec_store::ptr store);

	size_t get_num_bytes(off_t idx) const {
		return vec_offs[idx + 1] - vec_offs[idx];
	}

	size_t get_length(off_t idx) const {
		return get_num_bytes(idx) / get_type().get_size();
	}

	size_t get_num_vecs() const {
		return vec_offs.size() - 1;
	}

	virtual size_t get_num_bytes() const {
		return store->get_length() * store->get_type().get_size();
	}
	vec_store::const_ptr cat() const;

	virtual vec_store::ptr deep_copy() const;

	virtual bool append(const vec_store &vec);
	virtual bool append(std::vector<vec_store::const_ptr>::const_iterator vec_it,
			std::vector<vec_store::const_ptr>::const_iterator vec_end);
};

}

}

#endif
