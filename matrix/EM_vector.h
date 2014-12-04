#ifndef __EM_VECTOR_H__
#define __EM_VECTOR_H__

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
#include <atomic>

#include "io_interface.h"

namespace fm
{

class subvec_compute
{
public:
	typedef std::shared_ptr<subvec_compute> ptr;
	virtual ~subvec_compute() {
	}
	virtual void run(char *buf, size_t size) = 0;
};

class EM_vector;

struct fetch_vec_request
{
	size_t start;
	size_t length;
	subvec_compute::ptr compute;
};

struct set_vec_request
{
	const char *buf;
	size_t start;
	size_t length;
	subvec_compute::ptr compute;
};

/*
 * This is a per-thread data structure that helps access the vector in
 * the external memory.
 */
class EM_vector_accessor
{
	EM_vector &vec;
	safs::io_interface::ptr io;
public:
	typedef std::shared_ptr<EM_vector_accessor> ptr;

	EM_vector_accessor(EM_vector &_vec, safs::file_io_factory::shared_ptr factory);
	~EM_vector_accessor();

	/*
	 * Fetch a subvector from the original vector in [start, start + lenth).
	 * It performs computation on the fetched subvector asynchronously.
	 */
	void fetch_subvec(size_t start, size_t length, subvec_compute::ptr compute);
	void fetch_subvecs(const fetch_vec_request reqs[], size_t num);
	/*
	 * Store a subvector to the original vector in [start, start + lenth).
	 * It performs computation asynchronously when the data is stored on
	 * the original vector.
	 */
	void set_subvec(const char *buf, size_t start, size_t length,
			subvec_compute::ptr compute);
	void set_subvecs(const set_vec_request reqs[], size_t num);

	void wait4complete(int num);
	void wait4all();
};

class EM_vector
{
	std::atomic<size_t> accessor_count;
	safs::file_io_factory::shared_ptr factory;
	size_t length;
	size_t entry_size;

	EM_vector(size_t length, size_t entry_size);
public:
	typedef std::shared_ptr<EM_vector> ptr;

	static ptr create(size_t length, size_t entry_size) {
		return ptr(new EM_vector(length, entry_size));
	}

	~EM_vector();

	size_t get_size() const {
		return length;
	}

	size_t get_entry_size() const {
		return entry_size;
	}

	size_t get_byte_off(size_t entry_off) const {
		return entry_off * get_entry_size();
	}

	void resize(size_t length);

	EM_vector_accessor::ptr create_accessor() {
		accessor_count++;
		return EM_vector_accessor::ptr(new EM_vector_accessor(*this, factory));
	}

	/*
	 * Notify that a vector accessor is destroyed.
	 */
	void notify_destroy_accessor() {
		accessor_count--;
	}
};

}

#endif
