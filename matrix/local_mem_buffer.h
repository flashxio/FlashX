#ifndef __LOCAL_MEM_BUFFER_H__
#define __LOCAL_MEM_BUFFER_H__

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

#include <assert.h>

#include <atomic>
#include <memory>
#include <deque>
#include <vector>
#include <unordered_map>

#include "concurrency.h"

namespace fm
{

namespace detail
{

class local_matrix_store;

/*
 * This class keeps memory buffers in the local thread.
 *
 * If a piece of memory of the same size is being used again in the near
 * future, we can buffer the piece of memory and reuse it, so that we can
 * avoid allocating a large piece of memory from malloc.
 * It turns out to be fairly expensive to allocate a large piece of memory
 * with malloc in multiple threads. It can cause lock contention.
 *
 * If a portion of data in a matrix store will be used again (which is very
 * likely to happen in a chain of matrix computation), we should buffer it.
 *
 * The buffered data is guaranteed to be used in the same thread, so locking
 * isn't needed.
 */
class local_mem_buffer
{
	// The lock is to protect `mem_set'.
	static spin_lock mem_lock;
	// This contains all local mem buffers in the system.
	static std::vector<local_mem_buffer *> mem_set;
	static std::atomic<bool> initialized;
	static pthread_key_t mem_key;

	size_t num_allocs;
	size_t num_frees;
	std::unordered_map<size_t, std::deque<char *> > bufs;

	std::unordered_map<long, std::shared_ptr<const local_matrix_store> > portions;

	local_mem_buffer() {
		num_allocs = 0;
		num_frees = 0;
	}
	std::shared_ptr<char> _alloc(size_t num_bytes);

	void _cache_portion(long key,
			std::shared_ptr<const local_matrix_store> portion);
	std::shared_ptr<const local_matrix_store> _get_mat_portion(long key);

	void clear_local_bufs();
public:
	/*
	 * We initialize the memory buffers when the system starts to run and
	 * destroy them after system stops.
	 */
	static bool init();
	static void destroy();
	/*
	 * This function is called after each computation on data containers.
	 * We should delete the local buffer as much as possible to reduce
	 * memory consumption.
	 */
	static void clear_bufs();

	/*
	 * We cache matrix portions for EM matrix store and mapply virtual matrix
	 * store, so that we can reduce amount of data read from disks and
	 * save some computation.
	 */
	static void cache_portion(long key,
			std::shared_ptr<const local_matrix_store> portion);
	static std::shared_ptr<const local_matrix_store> get_mat_portion(long key);

	/*
	 * This function allocates memory from the memory buffer in the local thread.
	 */
	static std::shared_ptr<char> alloc(size_t num_bytes);

	~local_mem_buffer();
};

}

}

#endif
