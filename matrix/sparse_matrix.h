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
#include "matrix_io.h"

class compute_task
{
public:
	typedef std::shared_ptr<compute_task> ptr;

	virtual ~compute_task() {
	}

	virtual void run(char *buf, size_t size) = 0;
	virtual io_request get_request() const = 0;
};

class task_creator
{
public:
	typedef std::shared_ptr<task_creator> ptr;

	virtual compute_task::ptr create(const matrix_io &) const = 0;
};

class row_compute_task: public compute_task
{
	matrix_io io;
	char *buf;
	size_t buf_size;
public:
	row_compute_task(const matrix_io &_io): io(_io) {
		buf_size = ROUNDUP_PAGE(io.get_size());
		buf = (char *) valloc(buf_size);
	}

	~row_compute_task() {
		free(buf);
	}
	virtual void run(char *buf, size_t size);
	virtual void run_on_row(const ext_mem_undirected_vertex &v) = 0;
	virtual io_request get_request() const {
		return io_request(buf, io.get_loc(), buf_size, READ);
	}
};

template<class T>
class row_multiply_task: public row_compute_task
{
	const FG_vector<T> &input;
	FG_vector<T> &output;
public:
	row_multiply_task(const FG_vector<T> &_input, FG_vector<T> &_output,
			const matrix_io &_io): row_compute_task(_io), input(
				_input), output(_output) {
	}

	void run_on_row(const ext_mem_undirected_vertex &v);
};

template<class T>
void row_multiply_task<T>::run_on_row(const ext_mem_undirected_vertex &v)
{
	T res = 0;
	for (size_t i = 0; i < v.get_num_edges(); i++) {
		vertex_id_t id = v.get_neighbor(i);
		res += input.get(id);
	}
	output.set(v.get_id(), res);
}

template<class T>
class row_multiply_creator: public task_creator
{
	const FG_vector<T> &input;
	FG_vector<T> &output;

	row_multiply_creator(const FG_vector<T> &_input,
			FG_vector<T> &_output): input(_input), output(_output) {
	}
public:
	static task_creator::ptr create(const FG_vector<T> &_input,
			FG_vector<T> &_output) {
		return task_creator::ptr(new row_multiply_creator<T>(_input, _output));
	}

	virtual compute_task::ptr create(const matrix_io &io) const {
		return compute_task::ptr(new row_multiply_task<T>(input, output, io));
	}
};

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

	static ptr create(FG_graph::ptr);

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

	template<class T>
	typename FG_vector<T>::ptr multiply(typename FG_vector<T>::ptr in) const {
		if (in->get_size() != ncols) {
			BOOST_LOG_TRIVIAL(error) << boost::format(
					"the input vector has wrong length %1%. matrix ncols: %2%")
				% in->get_size() % ncols;
			return typename FG_vector<T>::ptr();
		}
		else {
			typename FG_vector<T>::ptr ret = FG_vector<T>::create(nrows);
			compute(row_multiply_creator<T>::create(*in, *ret));
			return ret;
		}
	}
};

void init_flash_matrix(config_map::ptr configs);
void destroy_flash_matrix();

#endif
