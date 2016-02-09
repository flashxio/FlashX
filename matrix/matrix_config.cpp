/*
 * Copyright 2016 Open Connectome Project (http://openconnecto.me)
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

#include <math.h>

#include <boost/format.hpp>

#include "thread.h"
#include "common.h"

#include "parameters.h"

#include "matrix_config.h"

namespace fm
{

void matrix_config::print_help()
{
	printf("Configuration parameters in matrix operations.\n");
	printf("\tthreads: the number of threads processing the matrix\n");
	printf("\tFM_prof_file: the output file containing CPU profiling\n");
	printf("\tin_mem_matrix: indicate whether to load the entire matrix to memory in advance\n");
	printf("\trow_block_size: the size of a row block (the number of rows)\n");
	printf("\trb_io_size: the size of a matrix I/O in 1D (the number of row blocks)\n");
	printf("\trb_steal_io_size: the size of a stolen matrix I/O(the number of row blocks)\n");
	printf("\tcpu_cache_size: the cpu cache size that can be used by a thread\n");
	printf("\thilbert_order: use the hilbert order\n");
	printf("\tnum_nodes: The number of NUMA nodes\n");
	printf("\tsort_buf_size: the buffer size for EM sorting\n");
	printf("\tgroupby_buf_size: the buffer size for EM groupby on vectors\n");
	printf("\tvv_groupby_buf_size: the buffer size for EM groupby on vector vectors\n");
	printf("\twrite_io_buf_size: the I/O buffer size for writing merge results\n");
	printf("\tstream_io_size: the I/O size used for streaming\n");
	printf("\tkeep_mem_buf: indicate whether to keep memory buffer for I/O in dense matrix operation\n");
}

void matrix_config::print()
{
	BOOST_LOG_TRIVIAL(info) << "Configuration parameters in matrix operations.";
	BOOST_LOG_TRIVIAL(info) << "\tSpM threads: " << num_SpM_threads;
	BOOST_LOG_TRIVIAL(info) << "\tDM threads: " << num_DM_threads;
	BOOST_LOG_TRIVIAL(info) << "\tFM_prof_file: " << prof_file;
	BOOST_LOG_TRIVIAL(info) << "\tin_mem_matrix: " << _in_mem_matrix;
	BOOST_LOG_TRIVIAL(info) << "\trow_block_size: " << row_block_size;
	BOOST_LOG_TRIVIAL(info) << "\trb_io_size: " << rb_io_size;
	BOOST_LOG_TRIVIAL(info) << "\trb_steal_io_size: " << rb_steal_io_size;
	BOOST_LOG_TRIVIAL(info) << "\tcpu_cache_size: " << cpu_cache_size;
	BOOST_LOG_TRIVIAL(info) << "\thilbert_order: " << hilbert_order;
	BOOST_LOG_TRIVIAL(info) << "\tnum_nodes: " << num_nodes;
	BOOST_LOG_TRIVIAL(info) << "\tsort_buf_size: " << sort_buf_size;
	BOOST_LOG_TRIVIAL(info) << "\tgroupby_buf_size: " << groupby_buf_size;
	BOOST_LOG_TRIVIAL(info) << "\tvv_groupby_buf_size: " << vv_groupby_buf_size;
	BOOST_LOG_TRIVIAL(info) << "\twrite_io_buf_size: " << write_io_buf_size;
	BOOST_LOG_TRIVIAL(info) << "\tstream_io_size: " << stream_io_size;
	BOOST_LOG_TRIVIAL(info) << "\tkeep_mem_buf: " << keep_mem_buf;
}

void matrix_config::init(config_map::ptr map)
{
	// `threads' sets #threads for sparse matrix and dense matrix operations.
	// But users can customize #threads for sparse matrix and dense matrix
	// separately.
	bool init_SpMT = false;
	bool init_DMT = false;
	if (map->has_option("threads")) {
		map->read_option_int("threads", num_SpM_threads);
		map->read_option_int("threads", num_DM_threads);
		init_SpMT = true;
		init_DMT = true;
	}
	if (map->has_option("SpM_threads")) {
		map->read_option_int("SpM_threads", num_SpM_threads);
		init_SpMT = true;
	}
	if (map->has_option("DM_threads")) {
		map->read_option_int("DM_threads", num_DM_threads);
		init_DMT = true;
	}
#ifdef USE_HWLOC
	if (!init_SpMT)
		num_SpM_threads = cpus.get_num_cores();
	if (!init_DMT)
		num_DM_threads = cpus.get_num_cores();
#endif

	if (map->has_option("FM_prof_file"))
		map->read_option("FM_prof_file", prof_file);
	if (map->has_option("in_mem_matrix"))
		map->read_option_bool("in_mem_matrix", _in_mem_matrix);
	if (map->has_option("row_block_size"))
		map->read_option_int("row_block_size", row_block_size);
	if (map->has_option("rb_io_size"))
		map->read_option_int("rb_io_size", rb_io_size);
	if (map->has_option("rb_steal_io_size"))
		map->read_option_int("rb_steal_io_size", rb_steal_io_size);
	if (map->has_option("cpu_cache_size"))
		map->read_option_int("cpu_cache_size", cpu_cache_size);
	if (map->has_option("hilbert_order"))
		map->read_option_bool("hilbert_order", hilbert_order);
	if (map->has_option("num_nodes"))
		map->read_option_int("num_nodes", num_nodes);
#ifdef USE_HWLOC
	else
		num_nodes = cpus.get_num_nodes();
#endif
	BOOST_LOG_TRIVIAL(info) << boost::format(
			"FlashMatrix runs on %1% threads and %2% NUMA nodes")
		% num_DM_threads % num_nodes;
	if (map->has_option("sort_buf_size")) {
		long tmp = 0;
		map->read_option_long("sort_buf_size", tmp);
		sort_buf_size = tmp;
	}
	if (map->has_option("groupby_buf_size")) {
		long tmp = 0;
		map->read_option_long("groupby_buf_size", tmp);
		groupby_buf_size = tmp;
	}
	if (map->has_option("vv_groupby_buf_size")) {
		long tmp = 0;
		map->read_option_long("vv_groupby_buf_size", tmp);
		vv_groupby_buf_size = tmp;
	}
	if (map->has_option("write_io_buf_size")) {
		long tmp = 0;
		map->read_option_long("write_io_buf_size", tmp);
		write_io_buf_size = tmp;
	}
	if (map->has_option("stream_io_size")) {
		long tmp = 0;
		map->read_option_long("stream_io_size", tmp);
		stream_io_size = tmp;
	}
	if (map->has_option("keep_mem_buf"))
		map->read_option_bool("keep_mem_buf", keep_mem_buf);
}
}
