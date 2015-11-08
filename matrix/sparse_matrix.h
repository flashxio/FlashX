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
 * This task performs sparse matrix dense matrix multiplication
 * in the FlashGraph format.
 * We implement this method for the sake of compatibility. It doesn't
 * run very fast.
 */
template<class DenseType, class SparseType, int ROW_WIDTH>
class fg_row_spmm_task: public fg_row_compute_task
{
	const detail::mem_matrix_store &input;
	detail::mem_matrix_store &output;
public:
	fg_row_spmm_task(const detail::mem_matrix_store &_input,
			detail::mem_matrix_store &_output,
			const matrix_io &_io): fg_row_compute_task(_io),
				input(_input), output(_output) {
		assert(input.get_type() == get_scalar_type<DenseType>());
		assert(output.get_type() == get_scalar_type<DenseType>());
		assert(input.get_num_cols() == output.get_num_cols());
		assert(input.get_num_cols() == (size_t) ROW_WIDTH);
	}

	void run_on_row(const fg::ext_mem_undirected_vertex &v);
};

template<class DenseType, class SparseType, int ROW_WIDTH>
void fg_row_spmm_task<DenseType, SparseType, ROW_WIDTH>::run_on_row(
		const fg::ext_mem_undirected_vertex &v)
{
	DenseType res[ROW_WIDTH];
	for (size_t i = 0; i < (size_t) ROW_WIDTH; i++)
		res[i] = 0;

	bool has_val = v.has_edge_data();
	for (size_t i = 0; i < v.get_num_edges(); i++) {
		fg::vertex_id_t id = v.get_neighbor(i);
		SparseType data = 1;
		if (has_val)
			data = v.get_edge_data<SparseType>(i);
		// TODO It's fairly expensive to get a row because it requires a function
		// call on a virtual method.
		const DenseType *row = (const DenseType *) input.get_row(id);
		for (size_t j = 0; j < (size_t) ROW_WIDTH; j++)
			res[j] += row[j] * data;
	}
	memcpy(output.get_row(v.get_id()), res, sizeof(DenseType) * ROW_WIDTH);
}

template<class DenseType, class SparseType>
class fg_row_spmm_task<DenseType, SparseType, 0>: public fg_row_compute_task
{
	const detail::mem_matrix_store &input;
	detail::mem_matrix_store &output;
public:
	fg_row_spmm_task(const detail::mem_matrix_store &_input,
			detail::mem_matrix_store &_output,
			const matrix_io &_io): fg_row_compute_task(_io),
				input(_input), output(_output) {
		assert(input.get_type() == get_scalar_type<DenseType>());
		assert(output.get_type() == get_scalar_type<DenseType>());
		assert(input.get_num_cols() == output.get_num_cols());
	}

	void run_on_row(const fg::ext_mem_undirected_vertex &v) {
		DenseType res[input.get_num_cols()];
		for (size_t i = 0; i < input.get_num_cols(); i++)
			res[i] = 0;

		bool has_val = v.has_edge_data();
		for (size_t i = 0; i < v.get_num_edges(); i++) {
			fg::vertex_id_t id = v.get_neighbor(i);
			SparseType data = 1;
			if (has_val)
				data = v.get_edge_data<SparseType>(i);
			// It's fairly expensive to get a row because it requires a function
			// call on a virtual method.
			const DenseType *row = (const DenseType *) input.get_row(id);
			for (size_t j = 0; j < input.get_num_cols(); j++)
				res[j] += row[j] * data;
		}
		memcpy(output.get_row(v.get_id()), res,
				sizeof(DenseType) * input.get_num_cols());
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
	size_t entry_size;
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
	// The size of the non-zero entries.
	size_t entry_size;

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

	size_t get_entry_size() const {
		return entry_size;
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

template<class DenseType, class SparseType, int ROW_WIDTH>
class row_part_func
{
	size_t row_width;
public:
	row_part_func(size_t row_width) {
		this->row_width = row_width;
		assert(ROW_WIDTH == row_width);
	}

	rp_edge_iterator operator()(rp_edge_iterator it,
			const char *_in_rows, char *_out_rows) {
		const DenseType *in_rows = (const DenseType *) _in_rows;
		DenseType *out_rows = (DenseType *) _out_rows;
		size_t row_idx = it.get_rel_row_idx();
		DenseType *dest_row = out_rows + ROW_WIDTH * row_idx;
		bool has_val = it.get_entry_size() > 0;
		while (it.has_next()) {
			SparseType data = 1;
			if (has_val)
				data = it.get_curr_data<SparseType>();
			size_t col_idx = it.next();
			const DenseType *src_row = in_rows + ROW_WIDTH * col_idx;
			for (size_t j = 0; j < ROW_WIDTH; j++)
				dest_row[j] += src_row[j] * data;
		}
		return it;
	}
};

template<class DenseType, class SparseType>
class row_part_func<DenseType, SparseType, 0>
{
	size_t row_width;
public:
	row_part_func(size_t row_width) {
		this->row_width = row_width;
	}

	rp_edge_iterator operator()(rp_edge_iterator it,
			const char *_in_rows, char *_out_rows) {
		const DenseType *in_rows = (const DenseType *) _in_rows;
		DenseType *out_rows = (DenseType *) _out_rows;
		size_t row_idx = it.get_rel_row_idx();
		DenseType *dest_row = out_rows + row_width * row_idx;
		bool has_val = it.get_entry_size() > 0;
		while (it.has_next()) {
			SparseType data = 1;
			if (has_val)
				data = it.get_curr_data<SparseType>();
			size_t col_idx = it.next();
			const DenseType *src_row = in_rows + row_width * col_idx;
			for (size_t j = 0; j < row_width; j++)
				dest_row[j] += src_row[j] * data;
		}
		return it;
	}
};

template<class DenseType, class SparseType, int ROW_WIDTH>
class coo_func
{
	size_t row_width;
public:
	coo_func(size_t row_width) {
		this->row_width = row_width;
		assert(ROW_WIDTH == row_width);
	}

	void operator()(const local_coo_t *coos, const char *_coo_vals,
			size_t num, const char *_in_rows, char *_out_rows) {
		const SparseType *coo_vals = (const SparseType *) _coo_vals;
		const DenseType *in_rows = (const DenseType *) _in_rows;
		DenseType *out_rows = (DenseType *) _out_rows;
		for (size_t i = 0; i < num; i++) {
			local_coo_t coo = coos[i];
			SparseType data = 1;
			if (coo_vals)
				data = coo_vals[i];
			const DenseType *src_row = in_rows + ROW_WIDTH * coo.get_col_idx();
			DenseType *dest_row = out_rows + ROW_WIDTH * coo.get_row_idx();
			for (size_t j = 0; j < ROW_WIDTH; j++)
				dest_row[j] += src_row[j] * data;
		}
	}
};

template<class DenseType, class SparseType>
class coo_func<DenseType, SparseType, 0>
{
	size_t row_width;
public:
	coo_func(size_t row_width) {
		this->row_width = row_width;
	}

	void operator()(const local_coo_t *coos, const char *_coo_vals,
			size_t num, const char *_in_rows, char *_out_rows) {
		const SparseType *coo_vals = (const SparseType *) _coo_vals;
		const DenseType *in_rows = (const DenseType *) _in_rows;
		DenseType *out_rows = (DenseType *) _out_rows;
		for (size_t i = 0; i < num; i++) {
			local_coo_t coo = coos[i];
			SparseType data = 1;
			if (coo_vals)
				data = coo_vals[i];
			const DenseType *src_row = in_rows + row_width * coo.get_col_idx();
			DenseType *dest_row = out_rows + row_width * coo.get_row_idx();
			for (size_t j = 0; j < row_width; j++)
				dest_row[j] += src_row[j] * data;
		}
	}
};

template<class RpFuncType, class COOFuncType>
class block_spmm_task_impl: public block_spmm_task
{
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
		size_t row_width = get_out_matrix().get_num_cols();

		if (block.has_rparts()) {
			RpFuncType rp_func(row_width);
			rp_edge_iterator it = block.get_first_edge_iterator(get_entry_size());
			while (!block.is_rparts_end(it)) {
				it = rp_func(it, in_rows, out_rows);
				it = block.get_next_edge_iterator(it, get_entry_size());
			}
		}
		const char *coo_vals = NULL;
		COOFuncType coo_func(row_width);
		if (get_entry_size() > 0)
			coo_vals = block.get_coo_val_start(get_entry_size());
		coo_func(block.get_coo_start(), coo_vals, block.get_num_coo_vals(),
				in_rows, out_rows);
	}
};

template<class DenseType, class SparseType>
class spmm_creator: public task_creator
{
	const detail::mem_matrix_store &input;
	detail::mem_matrix_store &output;
	const sparse_matrix &mat;
	block_exec_order::ptr order;

	spmm_creator(const detail::mem_matrix_store &_input,
			detail::mem_matrix_store &_output, const sparse_matrix &_mat);

	template<int ROW_WIDTH>
	compute_task::ptr create_block_compute_task(const matrix_io &io) const {
		typedef row_part_func<DenseType, SparseType, ROW_WIDTH> width_rp_func;
		typedef coo_func<DenseType, SparseType, ROW_WIDTH> width_coo_func;
		return compute_task::ptr(
				new block_spmm_task_impl<width_rp_func,  width_coo_func>(
					input, output, io, mat, order));
	}
public:
	static task_creator::ptr create(const detail::mem_matrix_store &_input,
			detail::mem_matrix_store &_output, const sparse_matrix &mat) {
		if (_input.get_type() != get_scalar_type<DenseType>()
				|| _output.get_type() != get_scalar_type<DenseType>()) {
			BOOST_LOG_TRIVIAL(error) << "wrong matrix type in spmm creator";
			return task_creator::ptr();
		}
		return task_creator::ptr(new spmm_creator<DenseType, SparseType>(
					_input, _output, mat));
	}

	virtual compute_task::ptr create(const matrix_io &io) const {
		if (order) {
			switch (output.get_num_cols()) {
				case 1: return create_block_compute_task<1>(io);
				case 2: return create_block_compute_task<2>(io);
				case 4: return create_block_compute_task<4>(io);
				case 8: return create_block_compute_task<8>(io);
				case 16: return create_block_compute_task<16>(io);
				default: return create_block_compute_task<0>(io);
			}
		}
		else {
			switch (output.get_num_cols()) {
				case 1: return compute_task::ptr(new fg_row_spmm_task<DenseType,
								SparseType, 1>(input, output, io));
				case 2: return compute_task::ptr(new fg_row_spmm_task<DenseType,
								SparseType, 2>(input, output, io));
				case 4: return compute_task::ptr(new fg_row_spmm_task<DenseType,
								SparseType, 4>(input, output, io));
				case 8: return compute_task::ptr(new fg_row_spmm_task<DenseType,
								SparseType, 8>(input, output, io));
				case 16: return compute_task::ptr(new fg_row_spmm_task<DenseType,
								 SparseType, 16>(input, output, io));
				default: return compute_task::ptr(new fg_row_spmm_task<DenseType,
								 SparseType, 0>(input, output, io));
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
	// The type of a non-zero entry.
	const scalar_type *entry_type;
	bool symmetric;

	template<class DenseType, class SparseType>
	task_creator::ptr get_multiply_creator(const detail::mem_matrix_store &in,
			detail::mem_matrix_store &out) const {
		return spmm_creator<DenseType, SparseType>::create(in, out, *this);
	}
protected:
	// This constructor is used for the sparse matrix stored
	// in the FlashGraph format.
	sparse_matrix(size_t num_vertices, const scalar_type *entry_type,
			bool symmetric) {
		this->nrows = num_vertices;
		this->ncols = num_vertices;
		this->entry_type = entry_type;
		this->symmetric = symmetric;
		this->is_fg = true;
	}

	sparse_matrix(size_t nrows, size_t ncols, const scalar_type *entry_type,
			bool symmetric) {
		this->symmetric = symmetric;
		this->is_fg = false;
		this->nrows = nrows;
		this->ncols = ncols;
		this->entry_type = entry_type;
	}
public:
	typedef std::shared_ptr<sparse_matrix> ptr;

	virtual ~sparse_matrix() {
	}

	/*
	 * This creates a sparse matrix sotred in the FlashGraph format.
	 */
	static ptr create(fg::FG_graph::ptr, const scalar_type *entry_type);
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

	/*
	 * The size of a non-zero entry.
	 */
	size_t get_entry_size() const {
		// Binary matrix doesn't need to store the non-zero entry.
		if (entry_type == NULL || *entry_type == get_scalar_type<bool>())
			return 0;
		else
			return entry_type->get_size();
	}
	template<class T>
	bool is_type() const {
		if (entry_type == NULL)
			return false;
		else
			return *entry_type == get_scalar_type<T>();
	}

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
	 * This version of SpMM allows users to provide the output matrix.
	 * It requires users to initialize the output matrix.
	 */
	template<class DenseType, class SparseType>
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

		compute(get_multiply_creator<DenseType, SparseType>(
					dynamic_cast<const detail::mem_matrix_store &>(in),
					dynamic_cast<detail::mem_matrix_store &>(out)),
				cal_super_block_size(get_block_size(),
					sizeof(DenseType) * in.get_num_cols()));
		return true;
	}

	template<class DenseType, class SparseType>
	bool multiply(const detail::vec_store &in, detail::vec_store &out) const {
		detail::matrix_store::const_ptr in_mat = in.conv2mat(in.get_length(),
				1, true);
		detail::matrix_store::const_ptr out_mat = out.conv2mat(out.get_length(),
				1, true);
		return multiply<DenseType, SparseType>(*in_mat,
				const_cast<detail::matrix_store &>(*out_mat));
	}
};

template<class DenseType, class SparseType>
spmm_creator<DenseType, SparseType>::spmm_creator(
		const detail::mem_matrix_store &_input, detail::mem_matrix_store &_output,
		const sparse_matrix &_mat): input(_input), output(_output), mat(_mat)
{
	// This initialization only for 2D partitioned sparse matrix.
	if (!mat.is_fg_matrix()) {
		// We only handle the case the element size is 2^n.
		assert(1 << ((size_t) log2(sizeof(DenseType))) == sizeof(DenseType));
		// Hilbert curve requires that there are 2^n block rows and block columns.
		// If the input matrix doesn't have 2^n columns, we have to find a number
		// of 2^n that is close to the number of columns in the input matrix
		// to calculate the super block size, so that the super block has 2^n
		// block rows and columns.
		size_t num_cols = input.get_num_cols();
		num_cols = 1 << ((size_t) log2(num_cols));
		size_t sb_size = cal_super_block_size(mat.get_block_size(),
				sizeof(DenseType) * num_cols);
		assert(input.get_num_cols() == output.get_num_cols());
		order = mat.get_multiply_order(sb_size, sb_size);
	}
}

void init_flash_matrix(config_map::ptr configs);
void destroy_flash_matrix();

}

#endif
