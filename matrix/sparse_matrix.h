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
#include "io_interface.h"
#include "mem_vector.h"

namespace fm
{

class compute_task
{
public:
	typedef std::shared_ptr<compute_task> ptr;

	virtual ~compute_task() {
	}

	virtual void run(char *buf, size_t size) = 0;
	virtual safs::io_request get_request() const = 0;
};

class task_creator
{
public:
	typedef std::shared_ptr<task_creator> ptr;

	virtual compute_task::ptr create(const matrix_io &) const = 0;
};

/*
 * This task performs computation on a sparse matrix in the FlashGraph format.
 */
class fg_row_compute_task: public compute_task
{
	matrix_io io;
	off_t off;
	char *buf;
	size_t buf_size;
public:
	fg_row_compute_task(const matrix_io &_io): io(_io) {
		off_t orig_off = io.get_loc().get_offset();
		off = ROUND_PAGE(orig_off);
		buf_size = ROUNDUP_PAGE(orig_off - off + io.get_size());
		buf = (char *) valloc(buf_size);
	}

	~fg_row_compute_task() {
		free(buf);
	}
	virtual void run(char *buf, size_t size);
	virtual void run_on_row(const fg::ext_mem_undirected_vertex &v) = 0;
	virtual safs::io_request get_request() const {
		return safs::io_request(buf, safs::data_loc_t(io.get_loc().get_file_id(),
					off), buf_size, READ);
	}
};

/*
 * This task performs matrix vector multiplication on a sparse matrix
 * in the FlashGraph format.
 */
template<class T>
class fg_row_multiply_task: public fg_row_compute_task
{
	const type_mem_vector<T> &input;
	type_mem_vector<T> &output;
public:
	fg_row_multiply_task(const type_mem_vector<T> &_input, type_mem_vector<T> &_output,
			const matrix_io &_io): fg_row_compute_task(_io), input(
				_input), output(_output) {
	}

	void run_on_row(const fg::ext_mem_undirected_vertex &v);
};

template<class T>
void fg_row_multiply_task<T>::run_on_row(const fg::ext_mem_undirected_vertex &v)
{
	T res = 0;
	for (size_t i = 0; i < v.get_num_edges(); i++) {
		fg::vertex_id_t id = v.get_neighbor(i);
		res += input.get(id);
	}
	output.set(v.get_id(), res);
}

template<class T>
class fg_row_multiply_creator: public task_creator
{
	const type_mem_vector<T> &input;
	type_mem_vector<T> &output;

	fg_row_multiply_creator(const type_mem_vector<T> &_input,
			type_mem_vector<T> &_output): input(_input), output(_output) {
	}
public:
	static task_creator::ptr create(const type_mem_vector<T> &_input,
			type_mem_vector<T> &_output) {
		return task_creator::ptr(new fg_row_multiply_creator<T>(_input, _output));
	}

	virtual compute_task::ptr create(const matrix_io &io) const {
		return compute_task::ptr(new fg_row_multiply_task<T>(input, output, io));
	}
};

class sparse_matrix
{
	size_t nrows;
	size_t ncols;
	bool symmetric;
	safs::file_io_factory::shared_ptr factory;
protected:
	// This constructor is used for the sparse matrix stored
	// in the FlashGraph format.
	sparse_matrix(safs::file_io_factory::shared_ptr factory, size_t num_vertices,
			bool symmetric) {
		this->factory = factory;
		this->nrows = num_vertices;
		this->ncols = num_vertices;
		this->symmetric = symmetric;
	}

	safs::file_io_factory::shared_ptr get_io_factory() const {
		return factory;
	}
public:
	typedef std::shared_ptr<sparse_matrix> ptr;

	virtual ~sparse_matrix() {
	}

	static ptr create(fg::FG_graph::ptr);

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

	virtual void transpose() = 0;

	template<class T>
	task_creator::ptr get_multiply_creator(type_mem_vector<T> &in,
			type_mem_vector<T> &out) const {
		return fg_row_multiply_creator<T>::create(in, out);
	}

	template<class T>
	typename type_mem_vector<T>::ptr multiply(typename type_mem_vector<T>::ptr in) const {
		if (in->get_length() != ncols) {
			BOOST_LOG_TRIVIAL(error) << boost::format(
					"the input vector has wrong length %1%. matrix ncols: %2%")
				% in->get_length() % ncols;
			return typename type_mem_vector<T>::ptr();
		}
		else {
			typename type_mem_vector<T>::ptr ret = type_mem_vector<T>::create(nrows);
			compute(get_multiply_creator<T>(*in, *ret));
			return ret;
		}
	}

	template<class T>
	dense_matrix::ptr multiply(dense_matrix::ptr in) const {
		if (in->get_num_rows() != ncols) {
			BOOST_LOG_TRIVIAL(error) << boost::format(
					"the input matrix has wrong #rows %1%. matrix ncols: %2%")
				% in->get_num_rows() % ncols;
			return dense_matrix::ptr();
		}
		else if (!in->is_in_mem()) {
			BOOST_LOG_TRIVIAL(error) << "SpMM doesn't support EM dense matrix";
			return dense_matrix::ptr();
		}
		else if (in->store_layout() == matrix_layout_t::L_ROW) {
			BOOST_LOG_TRIVIAL(error)
				<< "SpMM doesn't support row-wise dense matrix";
			return dense_matrix::ptr();
		}
		else {
			size_t ncol = in->get_num_cols();
			mem_col_dense_matrix::ptr col_m = mem_col_dense_matrix::cast(in);
			mem_col_dense_matrix::ptr ret = mem_col_dense_matrix::create(
					get_num_rows(), ncol, sizeof(T));
			std::vector<off_t> col_idx(1);
			for (size_t i = 0; i < ncol; i++) {
				col_idx[0] = i;
				typename type_mem_vector<T>::ptr in_col = type_mem_vector<T>::create(
						mem_dense_matrix::cast(col_m->get_cols(col_idx)));
				typename type_mem_vector<T>::ptr out_col = type_mem_vector<T>::create(
						mem_dense_matrix::cast(ret->get_cols(col_idx)));
				compute(get_multiply_creator<T>(*in_col, *out_col));
			}
			return ret;
		}
	}
};

void init_flash_matrix(config_map::ptr configs);
void destroy_flash_matrix();

}

#endif
