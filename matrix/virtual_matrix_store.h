#ifndef __VIRTUAL_MATRIX_STORE_H__
#define __VIRTUAL_MATRIX_STORE_H__

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

#include "mem_matrix_store.h"

namespace fm
{

namespace detail
{

/*
 * This class is the base class for the classes whose data isn't stored
 * physically. The data of these classes materializes data when being
 * requested. All these classes are read-only.
 */
class virtual_matrix_store: public mem_matrix_store
{
public:
	virtual_matrix_store(size_t nrow, size_t ncol,
			const scalar_type &type): mem_matrix_store(nrow, ncol, type) {
	}

	virtual matrix_store::ptr materialize() const = 0;

	virtual void reset_data() {
		assert(0);
	}

	virtual void set_data(const set_operate &op) {
		assert(0);
	}

	virtual char *get(size_t row, size_t col) {
		assert(0);
		return NULL;
	}

	virtual std::shared_ptr<local_matrix_store> get_portion(size_t id) {
		assert(0);
		return std::shared_ptr<local_matrix_store>();
	}

	virtual bool write2file(const std::string &file_name) const {
		assert(0);
		return false;
	}
};

}

}

#endif
