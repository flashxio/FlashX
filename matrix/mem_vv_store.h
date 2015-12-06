#ifndef __MEM_VV_STORE_H__
#define __MEM_VV_STORE_H__

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
#include "mem_vec_store.h"

namespace fm
{

class local_vv_store;

namespace detail
{

class matrix_store;

class mem_vv_store: public vv_store
{
	const smp_vec_store &get_mem_data() const {
		return static_cast<const smp_vec_store &>(get_data());
	}
	smp_vec_store &get_mem_data() {
		return static_cast<smp_vec_store &>(get_data());
	}
protected:
	mem_vv_store(const scalar_type &type): vv_store(type, true) {
	}
	mem_vv_store(const std::vector<off_t> &offs,
			mem_vec_store::ptr vec): vv_store(offs, vec) {
	}
public:
	typedef std::shared_ptr<mem_vv_store> ptr;
	typedef std::shared_ptr<const mem_vv_store> const_ptr;

	static ptr cast(vec_store::ptr store);

	static ptr create(const scalar_type &type) {
		return ptr(new mem_vv_store(type));
	}

	static ptr create(const std::vector<off_t> &offs, const scalar_type &type) {
		smp_vec_store::ptr vec = smp_vec_store::create(offs.back(), type);
		return ptr(new mem_vv_store(offs, vec));
	}

	static ptr create(const detail::raw_data_array &data,
			const std::vector<off_t> &offs, const scalar_type &type) {
		assert(offs.front() == 0 && (size_t) offs.back() == data.get_num_bytes());
		smp_vec_store::ptr vec = smp_vec_store::create(data, type);
		return ptr(new mem_vv_store(offs, vec));
	}

	virtual void set_data(const set_vv_operate &op);

	char *get_raw_arr() {
		return get_mem_data().get_raw_arr();
	}
	const char *get_raw_arr() const {
		return get_mem_data().get_raw_arr();
	}

	char*get_raw_arr(off_t idx) {
		return get_raw_arr() + get_vec_off(idx);
	}

	const char*get_raw_arr(off_t idx) const {
		return get_raw_arr() + get_vec_off(idx);
	}

	virtual vec_store::ptr shallow_copy() {
		// TODO the vector offsets may also be very large.
		return mem_vv_store::ptr(new mem_vv_store(*this));
	}
	virtual vec_store::const_ptr shallow_copy() const {
		// TODO the vector offsets may also be very large.
		return mem_vv_store::ptr(new mem_vv_store(*this));
	}
};

}

}

#endif
