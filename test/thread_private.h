#ifndef __THREAD_PRIVATE_H__
#define __THREAD_PRIVATE_H__

/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
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

#include "workload.h"
#include "io_interface.h"
#include "thread.h"
#include "config.h"

using namespace safs;

class cleanup_callback;
class sum_compute_allocator;
class write_compute_allocator;

/* this data structure stores the thread-private info. */
class thread_private: public thread
{
	/* the location in the thread descriptor array. */
	int idx;
	int node_id;
	workload_gen *gen;
	io_interface::ptr io;
	file_io_factory::shared_ptr factory;

	ssize_t read_bytes;
	long num_accesses;

	/* compute average number of pending IOs. */
	long tot_num_pending;
	long num_sampling;
	int max_num_pending;

	cleanup_callback *cb;
	
	struct timeval start_time, end_time;

	sum_compute_allocator *sum_alloc;
	write_compute_allocator *write_alloc;

#ifdef STATISTICS
public:
	atomic_integer num_completes;
	atomic_integer num_pending;
#endif

public:
	~thread_private() {
	}

	void init();

	int get_node_id() {
		return node_id;
	}

	thread_private(int node_id, int idx, int entry_size,
			file_io_factory::shared_ptr factory, workload_gen *gen);

	int attach2cpu();

	void run();

	ssize_t get_read_bytes();

	void print_stat();
};

#endif
