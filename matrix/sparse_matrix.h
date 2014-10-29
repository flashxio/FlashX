#ifndef __SPARSE_MATRIX_H__
#define __SPARSE_MATRIX_H__

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

#include "FGlib.h"

class compute_task
{
public:
	typedef std::shared_ptr<compute_task> ptr;

	virtual ~compute_task() {
	}

	virtual void run() = 0;
	virtual io_request get_request() const = 0;
};

class task_creator
{
public:
	typedef std::shared_ptr<task_creator> ptr;

	virtual compute_task::ptr create() const = 0;
};

class row_compute_task: public compute_task
{
	size_t start_row_id;
	size_t num_rows;

	char *buf;
	size_t buf_size;

	data_loc_t loc;
	size_t io_size;
public:
	row_compute_task() {
	}

	~row_compute_task() {
		free(buf);
	}
	virtual void run();
	virtual void run_on_row(const ext_mem_undirected_vertex &v) = 0;
};

template<class T>
class row_multiply_task: public row_compute_task
{
	FG_vector<T> &input;
	FG_vector<T> &output;
public:
	row_multiply_task() {
	}

	void run_on_row(const ext_mem_undirected_vertex &v);
};

template<class T>
void row_multiply_task<T>::run_on_row(const ext_mem_undirected_vertex &v)
{
	T res = 0;
	for (size_t i = 0; i < v.get_num_edges(); i++) {
		vertex_id_t id = v.get_neighbor(i);
		res += input[id];
	}
	output[v.get_id()] = res;
}

enum part_dim_t
{
	ROW,
	COL,
	BOTH,
};

class sparse_matrix
{
	size_t nrows;
	size_t ncols;
	bool symmetric;
	file_io_factory::shared_ptr factory;
	// On which dimension(s) the matrix is partitioned.
	part_dim_t part_dim;
protected:
	sparse_matrix(file_io_factory::shared_ptr factory, size_t nrows,
			size_t ncols, bool symmetric, part_dim_t part_dim) {
		this->factory = factory;
		this->nrows = nrows;
		this->ncols = ncols;
		this->symmetric = symmetric;
		this->part_dim = part_dim;
	}

	file_io_factory::shared_ptr get_io_factory() const {
		return factory;
	}
public:
	typedef std::shared_ptr<sparse_matrix> ptr;

	virtual ~sparse_matrix() {
	}

	ptr create(FG_graph::ptr);

	virtual void compute(task_creator::ptr creator) const = 0;

	size_t get_num_rows() const {
		return nrows;
	}

	size_t get_num_cols() const {
		return ncols;
	}

	bool is_symmetric() const {
		return symmetric;
	}

	int get_file_id() const {
		return factory->get_file_id();
	}

	part_dim_t get_part_dim() const {
		return part_dim;
	}
};

class row_block
{
	off_t off;
public:
	row_block(off_t off) {
		this->off = off;
	}

	off_t get_offset() const {
		return off;
	}
};

// The number of rows in a row block.
static const int ROW_BLOCK_SIZE = 1024;

#endif
