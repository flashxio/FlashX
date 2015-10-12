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
#include "mem_vec_store.h"
#include "NUMA_dense_matrix.h"
#include "dense_matrix.h"
#include "local_matrix_store.h"
#include "local_mem_buffer.h"

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
template<class T, class VectorType>
class fg_row_spmv_task: public fg_row_compute_task
{
	const VectorType &input;
	VectorType &output;
public:
	fg_row_spmv_task(const VectorType &_input, VectorType &_output,
			const matrix_io &_io): fg_row_compute_task(_io), input(
				_input), output(_output) {
	}

	void run_on_row(const fg::ext_mem_undirected_vertex &v);
};

template<class T, class VectorType>
void fg_row_spmv_task<T, VectorType>::run_on_row(
		const fg::ext_mem_undirected_vertex &v)
{
	T res = 0;
	for (size_t i = 0; i < v.get_num_edges(); i++) {
		fg::vertex_id_t id = v.get_neighbor(i);
		res += *(T *) input.get(id);
	}
	*(T *) output.get(v.get_id()) = res;
}

/*
 * This task performs sparse matrix dense matrix multiplication
 * in the FlashGraph format.
 * We implement this method for the sake of compatibility. It doesn't
 * run very fast.
 */
template<class T, int ROW_WIDTH>
class fg_row_spmm_task: public fg_row_compute_task
{
	const detail::mem_matrix_store &input;
	detail::mem_matrix_store &output;
public:
	fg_row_spmm_task(const detail::mem_matrix_store &_input,
			detail::mem_matrix_store &_output,
			const matrix_io &_io): fg_row_compute_task(_io),
				input(_input), output(_output) {
		assert(input.get_type() == get_scalar_type<T>());
		assert(output.get_type() == get_scalar_type<T>());
		assert(input.get_num_cols() == output.get_num_cols());
		assert(input.get_num_cols() == (size_t) ROW_WIDTH);
	}

	void run_on_row(const fg::ext_mem_undirected_vertex &v);
};

template<class T, int ROW_WIDTH>
void fg_row_spmm_task<T, ROW_WIDTH>::run_on_row(
		const fg::ext_mem_undirected_vertex &v)
{
	T res[ROW_WIDTH];
	for (size_t i = 0; i < (size_t) ROW_WIDTH; i++)
		res[i] = 0;

	for (size_t i = 0; i < v.get_num_edges(); i++) {
		fg::vertex_id_t id = v.get_neighbor(i);
		// It's fairly expensive to get a row because it requires a function
		// call on a virtual method.
		const T *row = (const T *) input.get_row(id);
		for (size_t j = 0; j < (size_t) ROW_WIDTH; j++)
			res[j] += row[j];
	}
	memcpy(output.get_row(v.get_id()), res, sizeof(T) * ROW_WIDTH);
}

template<class T>
class fg_row_spmm_task<T, 0>: public fg_row_compute_task
{
	const detail::mem_matrix_store &input;
	detail::mem_matrix_store &output;
public:
	fg_row_spmm_task(const detail::mem_matrix_store &_input,
			detail::mem_matrix_store &_output,
			const matrix_io &_io): fg_row_compute_task(_io),
				input(_input), output(_output) {
		assert(input.get_type() == get_scalar_type<T>());
		assert(output.get_type() == get_scalar_type<T>());
		assert(input.get_num_cols() == output.get_num_cols());
	}

	void run_on_row(const fg::ext_mem_undirected_vertex &v) {
		T res[input.get_num_cols()];
		for (size_t i = 0; i < input.get_num_cols(); i++)
			res[i] = 0;

		for (size_t i = 0; i < v.get_num_edges(); i++) {
			fg::vertex_id_t id = v.get_neighbor(i);
			// It's fairly expensive to get a row because it requires a function
			// call on a virtual method.
			const T *row = (const T *) input.get_row(id);
			for (size_t j = 0; j < input.get_num_cols(); j++)
				res[j] += row[j];
		}
		memcpy(output.get_row(v.get_id()), res, sizeof(T) * input.get_num_cols());
	}
};

class block_compute_task;

/*
 * The interface for processing blocks.
 */
class block_exec_order
{
public:
	typedef std::shared_ptr<block_exec_order> ptr;

	virtual bool is_valid_size(size_t height, size_t width) const = 0;

	// The vector of blocks should be from multiple block rows.
	virtual bool exec(block_compute_task &task,
			std::vector<const sparse_block_2d *> &blocks) const = 0;
};

/*
 * A compute task runs on a 2D-partitioned block of a sparse matrix.
 */
class block_compute_task: public compute_task
{
	block_exec_order::ptr exec_order;
	// The block rows in the buffer.
	std::vector<char *> block_rows;
	matrix_io io;
	off_t off;
	detail::local_mem_buffer::irreg_buf_t buf;
	size_t real_io_size;
protected:
	block_2d_size block_size;
public:
	block_compute_task(const matrix_io &_io, const sparse_matrix &mat,
			block_exec_order::ptr order);
	~block_compute_task();

	const matrix_io &get_io() const {
		return io;
	}

	virtual void run(char *buf, size_t size);

	virtual safs::io_request get_request() const {
		return safs::io_request(buf.second.get(),
				safs::data_loc_t(io.get_loc().get_file_id(), off),
				real_io_size, READ);
	}

	virtual void run_on_block(const sparse_block_2d &block) = 0;
	virtual void notify_complete() = 0;
};

class block_spmm_task: public block_compute_task
{
	const detail::mem_matrix_store &input;
	detail::mem_matrix_store &output;

	/*
	 * A task is responsible for processing the entire block rows.
	 * The result is stored in out_part.
	 */
	detail::local_row_matrix_store::ptr out_part;
public:
	block_spmm_task(const detail::mem_matrix_store &_input,
			detail::mem_matrix_store &_output, const matrix_io &io,
			const sparse_matrix &mat, block_exec_order::ptr order);

	const detail::mem_matrix_store &get_in_matrix() const {
		return input;
	}
	const detail::mem_matrix_store &get_out_matrix() const {
		return output;
	}

	/*
	 * Get the rows in the input matrix required by SpMM.
	 */
	const char *get_in_rows(size_t in_row, size_t num_rows);

	/*
	 * Get the rows in the output matrix required by SpMM.
	 */
	char *get_out_rows(size_t out_row, size_t num_rows);

	/*
	 * The entire block rows have been processed, the data in out_part is
	 * the final result, we can now save it in the output matrix.
	 */
	void notify_complete();
};

template<class T, int ROW_WIDTH>
class block_spmm_task_impl: public block_spmm_task
{
	rp_edge_iterator run_on_row_part(rp_edge_iterator it,
			const T *in_rows, T *out_rows) {
		size_t row_idx = it.get_rel_row_idx();
		T *dest_row = out_rows + ROW_WIDTH * row_idx;
		while (it.has_next()) {
			size_t col_idx = it.next();
			const T *src_row = in_rows + ROW_WIDTH * col_idx;
			for (size_t j = 0; j < ROW_WIDTH; j++)
				dest_row[j] += src_row[j];
		}
		return it;
	}
	void run_on_coo(const local_coo_t *coos, size_t num,
			const T *in_rows, T *out_rows) {
		for (size_t i = 0; i < num; i++) {
			local_coo_t coo = coos[i];
			const T *src_row = in_rows + ROW_WIDTH * coo.second /* col_idx */;
			T *dest_row = out_rows + ROW_WIDTH * coo.first /* row_idx */ ;
			for (size_t j = 0; j < ROW_WIDTH; j++)
				dest_row[j] += src_row[j];
		}
	}
public:
	block_spmm_task_impl(const detail::mem_matrix_store &input,
			detail::mem_matrix_store &output, const matrix_io &io,
			const sparse_matrix &mat,
			block_exec_order::ptr order): block_spmm_task(input,
				output, io, mat, order) {
		assert(ROW_WIDTH == output.get_num_cols());
	}

	void run_on_block(const sparse_block_2d &block) {
		if (block.is_empty())
			return;

		size_t in_row_start = block.get_block_col_idx() * block_size.get_num_cols();
		size_t num_in_rows = std::min(block_size.get_num_cols(),
				get_in_matrix().get_num_rows() - in_row_start);
		const char *in_rows = get_in_rows(in_row_start, num_in_rows);

		size_t out_row_start = block.get_block_row_idx() * block_size.get_num_rows();
		size_t num_out_rows = std::min(block_size.get_num_rows(),
				get_out_matrix().get_num_rows() - out_row_start);
		char *out_rows = get_out_rows(out_row_start, num_out_rows);

		if (block.has_rparts()) {
			rp_edge_iterator it = block.get_first_edge_iterator();
			while (!block.is_rparts_end(it)) {
				it = run_on_row_part(it, (const T *) in_rows, (T *) out_rows);
				it = block.get_next_edge_iterator(it);
			}
		}
		run_on_coo(block.get_coo_start(), block.get_num_coo_vals(),
				(const T *) in_rows, (T *) out_rows);
	}
};

template<class T>
class block_spmm_task_impl<T, 0>: public block_spmm_task
{
	rp_edge_iterator run_on_row_part(rp_edge_iterator it,
			const T *in_rows, T *out_rows) {
		size_t row_idx = it.get_rel_row_idx();
		size_t row_width = get_out_matrix().get_num_cols();
		T *dest_row = out_rows + row_width * row_idx;
		while (it.has_next()) {
			size_t col_idx = it.next();
			const T *src_row = in_rows + row_width * col_idx;
			for (size_t j = 0; j < row_width; j++)
				dest_row[j] += src_row[j];
		}
		return it;
	}
	void run_on_coo(const local_coo_t *coos, size_t num,
			const T *in_rows, T *out_rows) {
		size_t row_width = get_out_matrix().get_num_cols();
		for (size_t i = 0; i < num; i++) {
			local_coo_t coo = coos[i];
			const T *src_row = in_rows + row_width * coo.second /* col_idx */;
			T *dest_row = out_rows + row_width * coo.first /* row_idx */ ;
			for (size_t j = 0; j < row_width; j++)
				dest_row[j] += src_row[j];
		}
	}
public:
	block_spmm_task_impl(const detail::mem_matrix_store &input,
			detail::mem_matrix_store &output, const matrix_io &io,
			const sparse_matrix &mat,
			block_exec_order::ptr order): block_spmm_task(input,
				output, io, mat, order) {
	}

	void run_on_block(const sparse_block_2d &block) {
		if (block.is_empty())
			return;

		size_t in_row_start = block.get_block_col_idx() * block_size.get_num_cols();
		size_t num_in_rows = std::min(block_size.get_num_cols(),
				get_in_matrix().get_num_rows() - in_row_start);
		const char *in_rows = get_in_rows(in_row_start, num_in_rows);

		size_t out_row_start = block.get_block_row_idx() * block_size.get_num_rows();
		size_t num_out_rows = std::min(block_size.get_num_rows(),
				get_out_matrix().get_num_rows() - out_row_start);
		char *out_rows = get_out_rows(out_row_start, num_out_rows);

		if (block.has_rparts()) {
			rp_edge_iterator it = block.get_first_edge_iterator();
			while (!block.is_rparts_end(it)) {
				it = run_on_row_part(it, (const T *) in_rows, (T *) out_rows);
				it = block.get_next_edge_iterator(it);
			}
		}
		run_on_coo(block.get_coo_start(), block.get_num_coo_vals(),
				(const T *) in_rows, (T *) out_rows);
	}
};

/*
 * This task performs matrix vector multiplication on a sparse matrix in
 * a native format with 2D partitioning.
 */
template<class T>
class block_spmv_task: public block_compute_task
{
	const detail::mem_vec_store &input;
	detail::mem_vec_store &output;

	rp_edge_iterator run_on_row_part(rp_edge_iterator it,
			const T *in_arr, T *out_arr) {
		T sum = 0;
		while (it.has_next()) {
			size_t rel_col_idx = it.next();
			sum += in_arr[rel_col_idx];
		}
		out_arr[it.get_rel_row_idx()] += sum;
		return it;
	}
	void run_on_coo(const local_coo_t *coos, size_t num,
			const T *in_arr, T *out_arr) {
		for (size_t i = 0; i < num; i++) {
			local_coo_t coo = coos[i];
			out_arr[coo.first /* row_idx */] += in_arr[coo.second /* col_idx */];
		}
	}
public:
	block_spmv_task(const detail::mem_vec_store &_input,
			detail::mem_vec_store &_output, const matrix_io &_io,
			const sparse_matrix &mat, block_exec_order::ptr order): block_compute_task(
				_io, mat, order), input(_input), output(_output) {
	}

	void run_on_block(const sparse_block_2d &block) {
		if (block.is_empty())
			return;

		size_t start_col_idx
			= block.get_block_col_idx() * block_size.get_num_cols();
		size_t num_cols = std::min(block_size.get_num_cols(),
				input.get_length() - start_col_idx);
		size_t start_row_idx
			= block.get_block_row_idx() * block_size.get_num_rows();
		size_t num_rows = std::min(block_size.get_num_rows(),
				output.get_length() - start_row_idx);
		const char *in_buf = input.get_sub_arr(start_col_idx,
				start_col_idx + num_cols);
		char *out_buf = output.get_sub_arr(start_row_idx,
				start_row_idx + num_rows);
		assert(in_buf);
		assert(out_buf);
		if (block.has_rparts()) {
			rp_edge_iterator it = block.get_first_edge_iterator();
			while (!block.is_rparts_end(it)) {
				it = run_on_row_part(it, (const T *) in_buf, (T *) out_buf);
				it = block.get_next_edge_iterator(it);
			}
		}
		run_on_coo(block.get_coo_start(), block.get_num_coo_vals(),
				(const T *) in_buf, (T *) out_buf);
	}

	void notify_complete() {
	}
};

template<class T, class VectorType>
class fg_row_spmv_creator: public task_creator
{
	const VectorType &input;
	VectorType &output;

	fg_row_spmv_creator(const VectorType &_input,
			VectorType &_output): input(_input), output(_output) {
	}
public:
	static task_creator::ptr create(const VectorType &_input,
			VectorType &_output) {
		if (_input.get_type() != get_scalar_type<T>()
				|| _output.get_type() != get_scalar_type<T>()) {
			BOOST_LOG_TRIVIAL(error) << "wrong vector type in spmv creator";
			return task_creator::ptr();
		}
		return task_creator::ptr(new fg_row_spmv_creator<T, VectorType>(
					_input, _output));
	}

	virtual compute_task::ptr create(const matrix_io &io) const {
		return compute_task::ptr(new fg_row_spmv_task<T, VectorType>(input,
					output, io));
	}
};

template<class T>
class b2d_spmv_creator: public task_creator
{
	const detail::mem_vec_store &input;
	detail::mem_vec_store &output;
	const sparse_matrix &mat;
	block_exec_order::ptr order;

	b2d_spmv_creator(const detail::mem_vec_store &_input,
			detail::mem_vec_store &_output, const sparse_matrix &_mat);
public:
	static task_creator::ptr create(const detail::mem_vec_store &_input,
			detail::mem_vec_store &_output, const sparse_matrix &mat) {
		if (_input.get_type() != get_scalar_type<T>()
				|| _output.get_type() != get_scalar_type<T>()) {
			BOOST_LOG_TRIVIAL(error) << "wrong vector type in spmv creator";
			return task_creator::ptr();
		}
		return task_creator::ptr(new b2d_spmv_creator<T>(_input, _output, mat));
	}

	virtual compute_task::ptr create(const matrix_io &io) const {
		return compute_task::ptr(new block_spmv_task<T>(input, output, io, mat,
					order));
	}
};

template<class T>
class spmm_creator: public task_creator
{
	const detail::mem_matrix_store &input;
	detail::mem_matrix_store &output;
	const sparse_matrix &mat;
	block_exec_order::ptr order;

	spmm_creator(const detail::mem_matrix_store &_input,
			detail::mem_matrix_store &_output, const sparse_matrix &_mat);
public:
	static task_creator::ptr create(const detail::mem_matrix_store &_input,
			detail::mem_matrix_store &_output, const sparse_matrix &mat) {
		if (_input.get_type() != get_scalar_type<T>()
				|| _output.get_type() != get_scalar_type<T>()) {
			BOOST_LOG_TRIVIAL(error) << "wrong matrix type in spmm creator";
			return task_creator::ptr();
		}
		return task_creator::ptr(new spmm_creator<T>(_input, _output, mat));
	}

	virtual compute_task::ptr create(const matrix_io &io) const {
		if (order) {
			switch (output.get_num_cols()) {
				case 1:
					return compute_task::ptr(new block_spmm_task_impl<T, 1>(
								input, output, io, mat, order));
				case 2:
					return compute_task::ptr(new block_spmm_task_impl<T, 2>(
								input, output, io, mat, order));
				case 4:
					return compute_task::ptr(new block_spmm_task_impl<T, 4>(
								input, output, io, mat, order));
				case 8:
					return compute_task::ptr(new block_spmm_task_impl<T, 8>(
								input, output, io, mat, order));
				case 16:
					return compute_task::ptr(new block_spmm_task_impl<T, 16>(
								input, output, io, mat, order));
				case 32:
					return compute_task::ptr(new block_spmm_task_impl<T, 32>(
								input, output, io, mat, order));
				case 64:
					return compute_task::ptr(new block_spmm_task_impl<T, 64>(
								input, output, io, mat, order));
				case 128:
					return compute_task::ptr(new block_spmm_task_impl<T, 128>(
								input, output, io, mat, order));
				default:
					return compute_task::ptr(new block_spmm_task_impl<T, 0>(
								input, output, io, mat, order));
			}
		}
		else {
			switch (output.get_num_cols()) {
				case 1:
					return compute_task::ptr(new fg_row_spmm_task<T, 1>(input,
								output, io));
				case 2:
					return compute_task::ptr(new fg_row_spmm_task<T, 2>(input,
								output, io));
				case 4:
					return compute_task::ptr(new fg_row_spmm_task<T, 4>(input,
								output, io));
				case 8:
					return compute_task::ptr(new fg_row_spmm_task<T, 8>(input,
								output, io));
				case 16:
					return compute_task::ptr(new fg_row_spmm_task<T, 16>(input,
								output, io));
				case 32:
					return compute_task::ptr(new fg_row_spmm_task<T, 32>(input,
								output, io));
				case 64:
					return compute_task::ptr(new fg_row_spmm_task<T, 64>(input,
								output, io));
				case 128:
					return compute_task::ptr(new fg_row_spmm_task<T, 128>(input,
								output, io));
				default:
					return compute_task::ptr(new fg_row_spmm_task<T, 0>(input,
								output, io));
			}
		}
	}
};

/*
 * The non-zero entries in a sparse matrix are organized in blocks.
 * When multiplying the sparse matrix with other dense matrices, we traverse
 * the blocks in a certain order to increase CPU cache hits. To do so, we
 * need to group multiple adjacent blocks in a super block, which has
 * the same number of rows and columns. The size of the super block directly
 * affects CPU cache hits and is affected by the size of the dense matrix.
 * This function estimates the super block size in terms of the number of
 * block rows. The goal is to have the rows of the dense matrix involved
 * in the computation of the super block always in the CPU cache.
 * `entry_size' is the size of a row in the dense matrix.
 */
static inline size_t cal_super_block_size(const block_2d_size &block_size,
		size_t entry_size)
{
	// 1MB gives the best performance.
	// Maybe the reason is that I run eight threads in a processor, which has
	// a L3 cache of 8MB. Therefore, each thread gets 1MB in L3.
	size_t size = matrix_conf.get_cpu_cache_size() / entry_size
		/ block_size.get_num_rows() / 2;
	size_t max_size
		= detail::mem_matrix_store::CHUNK_SIZE / block_size.get_num_rows();
	return std::min(std::max(size, 1UL), max_size);
}

/*
 * This is a base class for a sparse matrix. It provides a set of functions
 * to perform computation on the sparse matrix. Currently, it has matrix
 * vector multiplication and matrix matrix multiplication.
 * We assume the sparse matrix is stored in external memory. If the matrix
 * is in memory, we can use in_mem_io to access the sparse matrix in memory
 * while reusing all the existing code for computation.
 */
class sparse_matrix
{
	// Whether the matrix is represented by the FlashGraph format.
	bool is_fg;
	size_t nrows;
	size_t ncols;
	bool symmetric;

	template<class T, class VectorType>
	task_creator::ptr get_fg_multiply_creator(const VectorType &in,
			VectorType &out) const {
		return fg_row_spmv_creator<T, VectorType>::create(in, out);
	}

	template<class T>
	task_creator::ptr get_multiply_creator(const detail::mem_vec_store &in,
			detail::mem_vec_store &out) const {
		return b2d_spmv_creator<T>::create(in, out, *this);
	}

	template<class T>
	task_creator::ptr get_multiply_creator(const detail::mem_matrix_store &in,
			detail::mem_matrix_store &out) const {
		return spmm_creator<T>::create(in, out, *this);
	}
protected:
	// This constructor is used for the sparse matrix stored
	// in the FlashGraph format.
	sparse_matrix(size_t num_vertices, bool symmetric) {
		this->nrows = num_vertices;
		this->ncols = num_vertices;
		this->symmetric = symmetric;
		this->is_fg = true;
	}

	sparse_matrix(size_t nrows, size_t ncols, bool symmetric) {
		this->symmetric = symmetric;
		this->is_fg = false;
		this->nrows = nrows;
		this->ncols = ncols;
	}
public:
	typedef std::shared_ptr<sparse_matrix> ptr;

	virtual ~sparse_matrix() {
	}

	/*
	 * This creates a sparse matrix sotred in the FlashGraph format.
	 */
	static ptr create(fg::FG_graph::ptr);
	/*
	 * This create a symmetric sparse matrix partitioned in 2D dimensions.
	 * The sparse matrix is stored in memory.
	 */
	static ptr create(SpM_2d_index::ptr, SpM_2d_storage::ptr);
	/*
	 * This creates an asymmetric sparse matrix partitioned in 2D dimensions.
	 * The sparse matrix is stored in memory.
	 */
	static ptr create(SpM_2d_index::ptr index, SpM_2d_storage::ptr mat,
			SpM_2d_index::ptr t_index, SpM_2d_storage::ptr t_mat);

	/*
	 * These two functions creates a sparse matrix partitioned in 2D dimensions
	 * from the data stored in SAFS.
	 */
	static ptr create(SpM_2d_index::ptr index,
			safs::file_io_factory::shared_ptr mat_io_fac);
	static ptr create(SpM_2d_index::ptr index,
			safs::file_io_factory::shared_ptr mat_io_fac,
			SpM_2d_index::ptr t_index,
			safs::file_io_factory::shared_ptr t_mat_io_fac);

	bool is_fg_matrix() const {
		return is_fg;
	}

	/*
	 * When customizing computatin on a sparse matrix (regardless of
	 * the storage format and the computation), we need to define two things:
	 * how data is accessed and what computation runs on the data fetched
	 * from the external memory. matrix_io_generator defines data access
	 * and compute_task defines the computation.
	 * `num_block_rows' will affect the behavior of matrix_io_generator.
	 * More explanation of `num_block_rows' is seen in the comments of
	 * `init_io_gens'.
	 */
	void compute(task_creator::ptr creator, size_t num_block_rows) const;
	/*
	 * This method defines how data in the matrix is accessed.
	 * Since we need to perform computation in parallel, we need to have
	 * a matrix I/O generator for each thread.
	 * `num_block_rows' defines how many block rows are accessed in an I/O.
	 */
	virtual void init_io_gens(size_t num_block_rows,
			std::vector<matrix_io_generator::ptr> &io_gens) const = 0;

	virtual safs::file_io_factory::shared_ptr get_io_factory() const = 0;
	/*
	 * Get the block size.
	 * For a row-major or column-major matrix, a block has only one row
	 * or one column.
	 */
	virtual const block_2d_size &get_block_size() const = 0;
	/*
	 * Get the offsets of the block rows.
	 */
	virtual void get_block_row_offs(const std::vector<off_t> &block_row_idxs,
			std::vector<off_t> &offs) const = 0;
	/*
	 * The subclass defines the order of processing a set of blocks when
	 * multiplying this sparse matrix with a dense matrix.
	 * All blocks are organized in a rectangular area.
	 */
	virtual block_exec_order::ptr get_multiply_order(
			size_t num_block_rows, size_t num_block_cols) const = 0;

	size_t get_num_rows() const {
		return nrows;
	}

	size_t get_num_cols() const {
		return ncols;
	}

	bool is_symmetric() const {
		return symmetric;
	}

	virtual void transpose() {
		size_t tmp = nrows;
		nrows = ncols;
		ncols = tmp;
	}

	/*
	 * This version of SpMV allows users to provide the output vector.
	 * It requires users to initialize the output vector.
	 */
	template<class T>
	bool multiply(const detail::mem_vec_store &in,
			detail::mem_vec_store &out) const {
		if (in.get_length() != ncols) {
			BOOST_LOG_TRIVIAL(error) << boost::format(
					"the input vector has wrong length %1%. matrix ncols: %2%")
				% in.get_length() % ncols;
			return false;
		}
		if (is_fg && in.get_num_nodes() >= 0)
			compute(get_fg_multiply_creator<T, detail::NUMA_vec_store>(
						dynamic_cast<const detail::NUMA_vec_store &>(in),
						dynamic_cast<detail::NUMA_vec_store &>(out)),
					cal_super_block_size(get_block_size(), sizeof(T)));
		else if (is_fg && in.get_num_nodes() < 0)
			compute(get_fg_multiply_creator<T, detail::smp_vec_store>(
						dynamic_cast<const detail::smp_vec_store &>(in),
						dynamic_cast<detail::smp_vec_store &>(out)),
					cal_super_block_size(get_block_size(), sizeof(T)));
		else
			compute(get_multiply_creator<T>(in, out),
					cal_super_block_size(get_block_size(), sizeof(T)));
		return true;
	}

	/*
	 * This version of SpMM allows users to provide the output matrix.
	 * It requires users to initialize the output matrix.
	 */
	template<class T>
	bool multiply(const detail::matrix_store &in, detail::matrix_store &out) const {
		if (in.get_num_rows() != ncols
				|| in.get_num_cols() != out.get_num_cols()
				|| out.get_num_rows() != this->get_num_rows()) {
			BOOST_LOG_TRIVIAL(error) <<
					"the input and output matrix have incompatible dimensions";
			return false;
		}
		else if (!in.is_in_mem() || !out.is_in_mem()) {
			BOOST_LOG_TRIVIAL(error) << "SpMM doesn't support EM dense matrix";
			return false;
		}
		// We allow the output matrix in column major.
		// TODO we should also allow the input matrix to be column-major.
		else if (in.store_layout() != matrix_layout_t::L_ROW) {
			BOOST_LOG_TRIVIAL(error) << "the dense matrix isn't row major.";
			return false;
		}

		compute(get_multiply_creator<T>(
					dynamic_cast<const detail::mem_matrix_store &>(in),
					dynamic_cast<detail::mem_matrix_store &>(out)),
				cal_super_block_size(get_block_size(),
					sizeof(T) * in.get_num_cols()));
		return true;
	}
};

template<class T>
b2d_spmv_creator<T>::b2d_spmv_creator(const detail::mem_vec_store &_input,
		detail::mem_vec_store &_output, const sparse_matrix &_mat): input(
			_input), output(_output), mat(_mat)
{
	// We only handle the case the element size is 2^n.
	assert(1 << ((size_t) log2(sizeof(T))) == sizeof(T));
	size_t sb_size = cal_super_block_size(mat.get_block_size(), sizeof(T));
	order = mat.get_multiply_order(sb_size, sb_size);
}

template<class T>
spmm_creator<T>::spmm_creator(const detail::mem_matrix_store &_input,
		detail::mem_matrix_store &_output,
		const sparse_matrix &_mat): input(_input), output(_output), mat(_mat)
{
	// This initialization only for 2D partitioned sparse matrix.
	if (!mat.is_fg_matrix()) {
		// We only handle the case the element size is 2^n.
		assert(1 << ((size_t) log2(sizeof(T))) == sizeof(T));
		// Hilbert curve requires that there are 2^n block rows and block columns.
		// If the input matrix doesn't have 2^n columns, we have to find a number
		// of 2^n that is close to the number of columns in the input matrix
		// to calculate the super block size, so that the super block has 2^n
		// block rows and columns.
		size_t num_cols = input.get_num_cols();
		num_cols = 1 << ((size_t) log2(num_cols));
		size_t sb_size = cal_super_block_size(mat.get_block_size(),
				sizeof(T) * num_cols);
		assert(input.get_num_cols() == output.get_num_cols());
		order = mat.get_multiply_order(sb_size, sb_size);
	}
}

void init_flash_matrix(config_map::ptr configs);
void destroy_flash_matrix();

}

#endif
