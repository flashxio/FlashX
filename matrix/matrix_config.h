#ifndef __MATRIX_CONFIG_H__
#define __MATRIX_CONFIG_H__

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

#include "common.h"
#include "log.h"
#include "config_map.h"

#include "graph_exception.h"

namespace fm
{

/**
 * The data structure contains the configurations for matrix operations.
 */
class matrix_config
{
	// The number of threads for sparse matrix.
	int num_SpM_threads;
	// The number of threads for dense matrix.
	int num_DM_threads;
	std::string prof_file;
	bool _in_mem_matrix;
	// With 1D partition, a matrix is partitioned into row blocks.
	int row_block_size;
	// For 1D partition, each matrix I/O contains multiple row blocks.
	// The matrix I/O size in row blocks.
	int rb_io_size;
	// For 1D partition, the size of a matrix I/O stolen from another thread.
	int rb_steal_io_size;
	// The size of CPU cache that can be used by a thread. It affects
	// the super block size.
	int cpu_cache_size;
	// Indicate whether the hilbert order is enabled.
	bool hilbert_order;
	// The number of NUMA nodes.
	int num_nodes;
	// The buffer size used for external-memory sorting.
	// The number of bytes.
	size_t sort_buf_size;
	// The buffer size used for external-memory groupby on vectors.
	// The number of bytes.
	size_t groupby_buf_size;
	// The buffer size used for EM groupby on vector vectors.
	// The number of vectors.
	size_t vv_groupby_buf_size;
	// The I/O buffer size for writing merge results in sorting a vector.
	// The number of bytes.
	size_t write_io_buf_size;
	// The I/O size used for streaming.
	size_t stream_io_size;
	// Indicate whether we keep the memory buffer for I/O in dense matrix
	// operations. Allocating the memory buffer for every dense matrix operation
	// is expensive.
	bool keep_mem_buf;
	// The block size for a dense matrix.
	size_t block_size;
public:
	/**
	 * \brief The default constructor that set all configurations to
	 * their default values.
	 */
	matrix_config() {
		num_SpM_threads = 4;
		num_DM_threads = 4;
		_in_mem_matrix = false;
		row_block_size = 1024;
		rb_io_size = 1024;
		rb_steal_io_size = 1;
		cpu_cache_size = 1024 * 1024;
		hilbert_order = false;
		num_nodes = 1;
		sort_buf_size = 128 * 1024 * 1024;
		groupby_buf_size = 128 * 1024 * 1024;
		vv_groupby_buf_size = 1024 * 1024;
		write_io_buf_size = 128 * 1024 * 1024;
		stream_io_size = 128 * 1024 * 1024;
		keep_mem_buf = false;
		block_size = 32;
	}

	/**
	 * \brief Print the explanations of all configurations.
	 */
	void print_help();
	/**
	 * \brief Print the current values of all configurations.
	 */
	void print();

	/**
	 * \brief Set the configurations to the user-defined values.
	 */
	void init(config_map::ptr map);

	/**
	 * \brief Get the output file containing CPU profiling.
	 * \return the file name.
	 */
	const std::string &get_prof_file() const {
		return prof_file;
	}

	/**
	 * \brief Get the number of worker threads for sparse matrix.
	 * \return The number of worker threads for sparse matrix.
	 */
	int get_num_SpM_threads() const {
		return num_SpM_threads;
	}

	/**
	 * \brief Get the number of worker threads for dense matrix.
	 * \return The number of worker threads for dense matrix.
	 */
	int get_num_DM_threads() const {
		return num_DM_threads;
	}

	/**
	 * \brief Determine whether to use in-mem matrix data.
	 * \return true if we loads the entire matrix data in memory
	 * in advance.
	 */
	bool use_in_mem_matrix() const {
		return _in_mem_matrix;
	}

	/**
	 * \brief The size of a row block (the number of rows).
	 */
	int get_row_block_size() const {
		return row_block_size;
	}

	/**
	 * \brief The size of a matrix I/O in 1D partitioning
	 * (the number of row blocks).
	 */
	int get_rb_io_size() const {
		return rb_io_size;
	}

	/**
	 * \brief The size of a matrix I/O stolen from another thread
	 * (the number of row blocks).
	 */
	int get_rb_steal_io_size() const {
		return rb_steal_io_size;
	}

	size_t get_cpu_cache_size() const {
		return cpu_cache_size;
	}

	bool use_hilbert_order() const {
		return hilbert_order;
	}

	void set_cpu_cache_size(size_t size) {
		cpu_cache_size = size;
	}

	void set_hilbert_order(bool hilbert) {
		hilbert_order = hilbert;
	}

	void set_num_SpM_threads(int nthreads) {
		this->num_SpM_threads = nthreads;
	}

	void set_num_DM_threads(int nthreads) {
		this->num_DM_threads = nthreads;
	}

	void set_num_nodes(int num_nodes) {
		this->num_nodes = num_nodes;
	}

	int get_num_nodes() const {
		return num_nodes;
	}

	void set_sort_buf_size(size_t sort_buf_size) {
		this->sort_buf_size = sort_buf_size;
	}

	void set_groupby_buf_size(size_t groupby_buf_size) {
		this->groupby_buf_size = groupby_buf_size;
	}

	void set_vv_groupby_buf_size(size_t vv_groupby_buf_size) {
		this->vv_groupby_buf_size = vv_groupby_buf_size;
	}

	void set_write_io_buf_size(size_t write_io_buf_size) {
		this->write_io_buf_size = write_io_buf_size;
	}

	size_t get_sort_buf_size() const {
		return sort_buf_size;
	}

	size_t get_groupby_buf_size() const {
		return groupby_buf_size;
	}

	size_t get_vv_groupby_buf_size() const {
		return vv_groupby_buf_size;
	}

	size_t get_write_io_buf_size() const {
		return write_io_buf_size;
	}

	size_t get_stream_io_size() const {
		return stream_io_size;
	}

	bool is_keep_mem_buf() const {
		return keep_mem_buf;
	}

	size_t get_block_size() const {
		return block_size;
	}
};

extern matrix_config matrix_conf;

static const int PAGE_SIZE = 4096;

// TODO We need to try different range size to get better performance.
static const size_t NUMA_range_size_log = 18;

}

#endif
