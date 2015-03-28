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
#include "NUMA_vector.h"
#include "NUMA_dense_matrix.h"

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
class fg_row_spmv_task: public fg_row_compute_task
{
	const NUMA_vector &input;
	NUMA_vector &output;
public:
	fg_row_spmv_task(const NUMA_vector &_input, NUMA_vector &_output,
			const matrix_io &_io): fg_row_compute_task(_io), input(
				_input), output(_output) {
	}

	void run_on_row(const fg::ext_mem_undirected_vertex &v);
};

template<class T>
void fg_row_spmv_task<T>::run_on_row(const fg::ext_mem_undirected_vertex &v)
{
	T res = 0;
	for (size_t i = 0; i < v.get_num_edges(); i++) {
		fg::vertex_id_t id = v.get_neighbor(i);
		res += input.get<T>(id);
	}
	output.set<T>(v.get_id(), res);
}

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
	char *buf;
	size_t buf_size;
protected:
	block_2d_size block_size;
public:
	block_compute_task(const matrix_io &_io, const sparse_matrix &mat,
			block_exec_order::ptr order);

	~block_compute_task() {
		free(buf);
	}

	virtual void run(char *buf, size_t size);

	virtual safs::io_request get_request() const {
		return safs::io_request(buf, safs::data_loc_t(io.get_loc().get_file_id(),
					off), buf_size, READ);
	}

	virtual void run_on_block(const sparse_block_2d &block) = 0;
};

template<class T>
class block_spmm_task: public block_compute_task
{
	const NUMA_row_tall_dense_matrix &input;
	NUMA_row_tall_dense_matrix &output;

	rp_edge_iterator run_on_row_part(rp_edge_iterator it,
			const T *in_rows, T *out_rows) {
		size_t row_idx = it.get_rel_row_idx();
		size_t row_width = output.get_num_cols();
		T *dest_row = out_rows + row_width * row_idx;
		while (it.has_next()) {
			size_t col_idx = it.next();
			const T *src_row = in_rows + row_width * col_idx;
			for (size_t j = 0; j < row_width; j++)
				dest_row[j] += src_row[j];
		}
		return it;
	}
public:
	block_spmm_task(const NUMA_row_tall_dense_matrix &_input,
			NUMA_row_tall_dense_matrix &_output, const matrix_io &_io,
			const sparse_matrix &mat, block_exec_order::ptr order): block_compute_task(
				_io, mat, order), input(_input), output(_output) {
	}

	void run_on_block(const sparse_block_2d &block) {
		size_t start_col_idx
			= block.get_block_col_idx() * block_size.get_num_cols();
		size_t start_row_idx
			= block.get_block_row_idx() * block_size.get_num_rows();
		const char *in_rows = input.get_rows(start_col_idx,
				start_col_idx + block_size.get_num_cols());
		char *out_rows = output.get_rows(start_row_idx,
				start_row_idx + block_size.get_num_rows());
		rp_edge_iterator it = block.get_first_edge_iterator();
		while (!block.is_block_end(it)) {
			it = run_on_row_part(it, (const T *) in_rows, (T *) out_rows);
			it = block.get_next_edge_iterator(it);
		}
	}
};

/*
 * This task performs matrix vector multiplication on a sparse matrix in
 * a native format with 2D partitioning.
 */
template<class T>
class block_spmv_task: public block_compute_task
{
	const NUMA_vector &input;
	NUMA_vector &output;

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
public:
	block_spmv_task(const NUMA_vector &_input,
			NUMA_vector &_output, const matrix_io &_io,
			const sparse_matrix &mat, block_exec_order::ptr order): block_compute_task(
				_io, mat, order), input(_input), output(_output) {
	}

	void run_on_block(const sparse_block_2d &block) {
		size_t start_col_idx
			= block.get_block_col_idx() * block_size.get_num_cols();
		size_t start_row_idx
			= block.get_block_row_idx() * block_size.get_num_rows();
		rp_edge_iterator it = block.get_first_edge_iterator();
		const char *in_buf = input.get_sub_arr(start_col_idx,
				start_col_idx + block_size.get_num_cols());
		char *out_buf = output.get_sub_arr(start_row_idx,
				start_row_idx + block_size.get_num_rows());
		assert(in_buf);
		assert(out_buf);
		while (!block.is_block_end(it)) {
			it = run_on_row_part(it, (const T *) in_buf, (T *) out_buf);
			it = block.get_next_edge_iterator(it);
		}
	}
};

template<class T>
class fg_row_spmv_creator: public task_creator
{
	const NUMA_vector &input;
	NUMA_vector &output;

	fg_row_spmv_creator(const NUMA_vector &_input,
			NUMA_vector &_output): input(_input), output(_output) {
	}
public:
	static task_creator::ptr create(const NUMA_vector &_input,
			NUMA_vector &_output) {
		if (_input.get_type() != get_scalar_type<T>()
				|| _output.get_type() != get_scalar_type<T>()) {
			BOOST_LOG_TRIVIAL(error) << "wrong vector type in spmv creator";
			return task_creator::ptr();
		}
		return task_creator::ptr(new fg_row_spmv_creator<T>(_input, _output));
	}

	virtual compute_task::ptr create(const matrix_io &io) const {
		return compute_task::ptr(new fg_row_spmv_task<T>(input, output, io));
	}
};

template<class T>
class b2d_spmv_creator: public task_creator
{
	const NUMA_vector &input;
	NUMA_vector &output;
	const sparse_matrix &mat;
	block_exec_order::ptr order;

	b2d_spmv_creator(const NUMA_vector &_input,
			NUMA_vector &_output, const sparse_matrix &_mat);
public:
	static task_creator::ptr create(const NUMA_vector &_input,
			NUMA_vector &_output, const sparse_matrix &mat) {
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
class b2d_spmm_creator: public task_creator
{
	const NUMA_row_tall_dense_matrix &input;
	NUMA_row_tall_dense_matrix &output;
	const sparse_matrix &mat;
	block_exec_order::ptr order;

	b2d_spmm_creator(const NUMA_row_tall_dense_matrix &_input,
			NUMA_row_tall_dense_matrix &_output, const sparse_matrix &_mat);
public:
	static task_creator::ptr create(const NUMA_row_tall_dense_matrix &_input,
			NUMA_row_tall_dense_matrix &_output, const sparse_matrix &mat) {
		if (_input.get_type() != get_scalar_type<T>()
				|| _output.get_type() != get_scalar_type<T>()) {
			BOOST_LOG_TRIVIAL(error) << "wrong matrix type in spmm creator";
			return task_creator::ptr();
		}
		return task_creator::ptr(new b2d_spmm_creator<T>(_input, _output, mat));
	}

	virtual compute_task::ptr create(const matrix_io &io) const {
		return compute_task::ptr(new block_spmm_task<T>(input, output,
					io, mat, order));
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
		/ block_size.get_num_rows();
	return std::max(size, 1UL);
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

	template<class T>
	task_creator::ptr get_multiply_creator(const NUMA_vector &in,
			NUMA_vector &out) const {
		if (is_fg)
			return fg_row_spmv_creator<T>::create(in, out);
		else
			return b2d_spmv_creator<T>::create(in, out, *this);
	}

	template<class T>
	task_creator::ptr get_multiply_creator(const NUMA_row_tall_dense_matrix &in,
			NUMA_row_tall_dense_matrix &out) const {
		assert(!is_fg);
		return b2d_spmm_creator<T>::create(in, out, *this);
	}

	template<class T>
	void multiply_matrix(const NUMA_row_tall_dense_matrix &row_m,
			NUMA_row_tall_dense_matrix &ret) const {
		compute(get_multiply_creator<T>(row_m, ret),
				cal_super_block_size(get_block_size(),
					sizeof(T) * row_m.get_num_cols()));
	}

#if 0
	template<class T>
	void multiply_matrix(const mem_col_dense_matrix &col_m,
			mem_col_dense_matrix &ret) const {
		size_t ncol = col_m.get_num_cols();
		std::vector<off_t> col_idx(1);
		for (size_t i = 0; i < ncol; i++) {
			col_idx[0] = i;
			mem_vector::ptr in_col = mem_vector::create(
					mem_dense_matrix::cast(col_m.get_cols(col_idx)));
			mem_vector::ptr out_col = mem_vector::create(
					mem_dense_matrix::cast(ret.get_cols(col_idx)));
			compute(get_multiply_creator<T>(*in_col, *out_col),
				cal_super_block_size(get_block_size(),
					sizeof(T) * col_m.get_num_cols()));
		}
	}
#endif
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
	bool multiply(const NUMA_vector &in, NUMA_vector &out) const {
		if (in.get_length() != ncols) {
			BOOST_LOG_TRIVIAL(error) << boost::format(
					"the input vector has wrong length %1%. matrix ncols: %2%")
				% in.get_length() % ncols;
			return false;
		}
		compute(get_multiply_creator<T>(in, out),
				cal_super_block_size(get_block_size(), sizeof(T)));
		return true;
	}

	/*
	 * This version of SpMV allocates the output vector.
	 */
	template<class T>
	NUMA_vector::ptr multiply(NUMA_vector::ptr in) const {
		if (in->get_length() != ncols) {
			BOOST_LOG_TRIVIAL(error) << boost::format(
					"the input vector has wrong length %1%. matrix ncols: %2%")
				% in->get_length() % ncols;
			return NUMA_vector::ptr();
		}
		else {
			NUMA_vector::ptr ret = NUMA_vector::create(nrows,
					in->get_num_nodes(), get_scalar_type<T>());
			ret->reset_data();
			multiply<T>(*in, *ret);
			return ret;
		}
	}

	/*
	 * This version of SpMM allows users to provide the output matrix.
	 * It requires users to initialize the output matrix.
	 */
	template<class T>
	bool multiply(const dense_matrix &in, dense_matrix &out) const {
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

		if (in.store_layout() == matrix_layout_t::L_ROW) {
			if (out.store_layout() != matrix_layout_t::L_ROW) {
				BOOST_LOG_TRIVIAL(error) << "wrong matrix layout for output matrix";
				return false;
			}
			// TODO how are we going to tell the difference between mem_dense_matrix
			// and NUMA dense matrix.
			multiply_matrix<T>((const NUMA_row_tall_dense_matrix &) in,
					(NUMA_row_tall_dense_matrix &) out);
			return true;
		}
		else {
			if (out.store_layout() != matrix_layout_t::L_COL) {
				BOOST_LOG_TRIVIAL(error) << "wrong matrix layout for output matrix";
				return false;
			}
#if 0
			multiply_matrix<T>((const mem_col_dense_matrix &) in,
					(mem_col_dense_matrix &) out);
#endif
			return true;
		}
	}

	/*
	 * This version of SpMM allocates the output matrix.
	 */
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

		if (in->store_layout() == matrix_layout_t::L_ROW) {
			// TODO how are we going to tell the difference between mem_dense_matrix
			// and NUMA dense matrix.
			const NUMA_row_tall_dense_matrix &in_mat
				= (const NUMA_row_tall_dense_matrix &) *in;
			NUMA_row_tall_dense_matrix::ptr ret = NUMA_row_tall_dense_matrix::create(
					get_num_rows(), in->get_num_cols(), in_mat.get_num_nodes(),
					get_scalar_type<T>());
			ret->reset_data();
			multiply_matrix<T>(in_mat, *ret);
			return ret;
		}
		else {
#if 0
			mem_col_dense_matrix::ptr ret = mem_col_dense_matrix::create(
					get_num_rows(), in->get_num_cols(), get_scalar_type<T>());
			ret->reset_data();
			multiply_matrix<T>((const mem_col_dense_matrix &) *in, *ret);
			return ret;
#endif
		}
	}
};

template<class T>
b2d_spmv_creator<T>::b2d_spmv_creator(const NUMA_vector &_input,
		NUMA_vector &_output, const sparse_matrix &_mat): input(
			_input), output(_output), mat(_mat)
{
	// We only handle the case the element size is 2^n.
	assert(1 << ((size_t) log2(sizeof(T))) == sizeof(T));
	size_t sb_size = cal_super_block_size(mat.get_block_size(), sizeof(T));
	order = mat.get_multiply_order(sb_size, sb_size);
}

template<class T>
b2d_spmm_creator<T>::b2d_spmm_creator(const NUMA_row_tall_dense_matrix &_input,
		NUMA_row_tall_dense_matrix &_output, const sparse_matrix &_mat): input(
			_input), output(_output), mat(_mat)
{
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

void init_flash_matrix(config_map::ptr configs);
void destroy_flash_matrix();

}

#endif
