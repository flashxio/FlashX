#ifndef __EM_VV_STORE_H__
#define __EM_VV_STORE_H__

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

#include "vv_store.h"
#include "EM_object.h"
#include "EM_vector.h"

namespace fm
{

namespace detail
{

class EM_vv_store: public vv_store, public EM_object
{
	EM_vv_store(const scalar_type &type): vv_store(type, false) {
	}
	EM_vv_store(const std::vector<off_t> &offs,
			EM_vec_store::ptr store): vv_store(offs, store) {
	}
public:
	typedef std::shared_ptr<EM_vv_store> ptr;
	typedef std::shared_ptr<const EM_vv_store> const_ptr;

	static ptr create(const scalar_type &type) {
		return ptr(new EM_vv_store(type));
	}
	static ptr create(const std::vector<off_t> &offs, EM_vec_store::ptr store) {
		return ptr(new EM_vv_store(offs, store));
	}

	virtual safs::io_interface::ptr create_io() {
		printf("create I/O from %p in vv store\n", &get_data());
		return dynamic_cast<EM_object &>(get_data()).create_io();
	}

	virtual vec_store::ptr shallow_copy() {
		return ptr(new EM_vv_store(*this));
	}
	virtual vec_store::const_ptr shallow_copy() const {
		return const_ptr(new EM_vv_store(*this));
	}
};

}

}

#endif
