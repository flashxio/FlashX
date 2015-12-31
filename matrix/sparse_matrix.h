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
#include "EM_object.h"

namespace fm
{

namespace detail
{

class compute_task: public detail::portion_compute
{
public:
	typedef std::shared_ptr<compute_task> ptr;

	virtual safs::io_request get_request() const = 0;
};

class task_creator
{
public:
	typedef std::shared_ptr<task_creator> ptr;

	virtual compute_task::ptr create(const matrix_io &) const = 0;
	virtual bool set_data(detail::matrix_store::const_ptr in,
			detail::matrix_store::ptr out) = 0;
	virtual std::vector<detail::EM_object *> get_EM_objs() = 0;
	virtual bool is_complete() const = 0;
};

class row_portions
{
	size_t portion_size_log;
	size_t portion_mask;
	size_t num_cols;
	size_t entry_size;
	size_t tot_num_rows;
	std::vector<local_matrix_store::const_ptr> portions;
	std::vector<const char *> raw_portions;

	row_portions();
public:
	typedef std::shared_ptr<row_portions> ptr;

	static ptr create(matrix_store::const_ptr mat);

	const scalar_type &get_type() const {
		return portions[0]->get_type();
	}

	/*
	 * This method is called very frequently. Let's avoid the overhead
	 * of function calls by inlining the function.
	 */
	const char *get_row(size_t row_idx) const {
		size_t local_row_idx = row_idx & portion_mask;
		size_t portion_idx = row_idx >> portion_size_log;
		return raw_portions[portion_idx]
			+ local_row_idx * num_cols * entry_size;
	}
	const char *get_rows(size_t start_row, size_t end_row) const;
	size_t get_num_cols() const {
		return num_cols;
	}
	size_t get_tot_num_rows() const {
		return tot_num_rows;
	}
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
protected:
	// rows in the input matrix.
	row_portions::ptr in_row_portions;
	// rows in the output matrix.
	local_matrix_store::ptr out_rows;
	char *raw_out_rows;
public:
	fg_row_compute_task(matrix_store &output, const matrix_io &_io,
			row_portions::ptr in_row_portions): io(_io) {
		off_t orig_off = io.get_loc().get_offset();
		off = ROUND_PAGE(orig_off);
		buf_size = ROUNDUP_PAGE(orig_off - off + io.get_size());
		buf = (char *) valloc(buf_size);

		this->in_row_portions = in_row_portions;
		this->out_rows = output.get_portion(_io.get_top_left().get_row_idx(),
				0, _io.get_num_rows(), output.get_num_cols());
		assert(out_rows);
		this->raw_out_rows = out_rows->get_raw_arr();
		assert(raw_out_rows);
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
	matrix_store &output;
public:
	fg_row_spmm_task(row_portions::ptr in_row_portions,
			matrix_store &_output, const matrix_io &_io): fg_row_compute_task(
				_output, _io, in_row_portions), output(_output) {
		assert(in_row_portions->get_type() == get_scalar_type<DenseType>());
		assert(output.get_type() == get_scalar_type<DenseType>());
		assert(in_row_portions->get_num_cols() == output.get_num_cols());
		assert(in_row_portions->get_num_cols() == (size_t) ROW_WIDTH);
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
		const DenseType *row = (const DenseType *) in_row_portions->get_row(id);
		for (size_t j = 0; j < (size_t) ROW_WIDTH; j++)
			res[j] += row[j] * data;
	}
	size_t rel_row_idx = v.get_id() - out_rows->get_global_start_row();
	assert(rel_row_idx < out_rows->get_num_rows());
	char *row = raw_out_rows + rel_row_idx * ROW_WIDTH * sizeof(DenseType);
	memcpy(row, res, sizeof(DenseType) * ROW_WIDTH);
}

template<class DenseType, class SparseType>
class fg_row_spmm_task<DenseType, SparseType, 0>: public fg_row_compute_task
{
	matrix_store &output;
public:
	fg_row_spmm_task(row_portions::ptr in_row_portions,
			matrix_store &_output, const matrix_io &_io): fg_row_compute_task(
				_output, _io, in_row_portions), output(_output) {
		assert(in_row_portions->get_type() == get_scalar_type<DenseType>());
		assert(output.get_type() == get_scalar_type<DenseType>());
		assert(in_row_portions->get_num_cols() == output.get_num_cols());
	}

	void run_on_row(const fg::ext_mem_undirected_vertex &v) {
		DenseType res[in_row_portions->get_num_cols()];
		for (size_t i = 0; i < in_row_portions->get_num_cols(); i++)
			res[i] = 0;

		bool has_val = v.has_edge_data();
		for (size_t i = 0; i < v.get_num_edges(); i++) {
			fg::vertex_id_t id = v.get_neighbor(i);
			SparseType data = 1;
			if (has_val)
				data = v.get_edge_data<SparseType>(i);
			const DenseType *row = (const DenseType *) in_row_portions->get_row(id);
			for (size_t j = 0; j < in_row_portions->get_num_cols(); j++)
				res[j] += row[j] * data;
		}
		size_t rel_row_idx = v.get_id() - out_rows->get_global_start_row();
		assert(rel_row_idx < out_rows->get_num_rows());
		char *row = raw_out_rows
			+ rel_row_idx * output.get_num_cols() * sizeof(DenseType);
		memcpy(row, res, sizeof(DenseType) * output.get_num_cols());
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
 * A compute task runs on multiple 2D-partitioned block rows of a sparse matrix.
 */
class block_compute_task: public compute_task
{
	block_exec_order::ptr exec_order;
	// The block rows in the buffer.
	std::vector<char *> block_rows;
	matrix_io io;
	off_t off;
	local_mem_buffer::irreg_buf_t buf;
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

/*
 * This class helps to stream data to the output dense matrix from SpMM
 * on disks, so the dense matrix is a very tall and skinny matrix.
 * This assumes that the data comes from multiple threads without a specific
 * order. But once we order the incoming data, we can stream data to disks
 * sequentially.
 * This data structure is shared by multiple threads.
 */
class EM_matrix_stream
{
	matrix_store::ptr mat;

	class filled_local_store
	{
		local_matrix_store::ptr data;
		std::atomic<size_t> num_filled_rows;
	public:
		typedef std::shared_ptr<filled_local_store> ptr;

		filled_local_store(local_matrix_store::ptr data) {
			this->data = data;
			num_filled_rows = 0;
		}
		// If the write completely fill the local buffer, it returns true.
		bool write(local_matrix_store::const_ptr portion,
				off_t global_start_row, off_t global_start_col);

		local_matrix_store::const_ptr get_whole_portion() const {
			return data;
		}

		off_t get_global_start_row() const {
			return data->get_global_start_row();
		}
	};

	pthread_spinlock_t lock;
	// This keeps the buffers with partial data in EM matrix portions.
	// If an EM matrix portion is complete, the portion is flushed to disks
	// and it is deleted from the hashtable.
	std::unordered_map<off_t, filled_local_store::ptr> portion_bufs;

	EM_matrix_stream(matrix_store::ptr mat) {
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		this->mat = mat;
	}
public:
	typedef std::shared_ptr<EM_matrix_stream> ptr;

	static ptr create(matrix_store::ptr mat) {
		return ptr(new EM_matrix_stream(mat));
	}

	void write_async(local_matrix_store::const_ptr portion,
			off_t start_row, off_t start_col);
	bool is_complete() const;
};

class block_spmm_task: public block_compute_task
{
	// The size of the non-zero entries.
	size_t entry_size;

	matrix_store &output;
	EM_matrix_stream::ptr output_stream;

	/*
	 * A task is responsible for processing the entire block rows.
	 * The result is stored in out_part.
	 */
	local_row_matrix_store::ptr out_part;
protected:
	row_portions::ptr in_row_portions;
public:
	block_spmm_task(row_portions::ptr in_row_portions,
			matrix_store &_output, EM_matrix_stream::ptr stream,
			const matrix_io &io, const sparse_matrix &mat,
			block_exec_order::ptr order);

	const matrix_store &get_out_matrix() const {
		return output;
	}

	size_t get_entry_size() const {
		return entry_size;
	}

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
	block_spmm_task_impl(row_portions::ptr in_row_portions,
			matrix_store &output, EM_matrix_stream::ptr stream
			, const matrix_io &io, const sparse_matrix &mat,
			block_exec_order::ptr order): block_spmm_task(in_row_portions,
				output, stream, io, mat, order) {
	}

	void run_on_block(const sparse_block_2d &block) {
		if (block.is_empty())
			return;

		size_t in_row_start = block.get_block_col_idx() * block_size.get_num_cols();
		size_t num_in_rows = std::min(block_size.get_num_cols(),
				in_row_portions->get_tot_num_rows() - in_row_start);
		const char *in_rows = in_row_portions->get_rows(in_row_start,
				in_row_start + num_in_rows);

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
	row_portions::ptr in_row_portions;
	matrix_store::ptr output;
	EM_matrix_stream::ptr output_stream;
	const sparse_matrix &mat;
	block_exec_order::ptr order;

	spmm_creator(const sparse_matrix &_mat, size_t num_in_cols);

	template<int ROW_WIDTH>
	compute_task::ptr create_block_compute_task(const matrix_io &io) const {
		typedef row_part_func<DenseType, SparseType, ROW_WIDTH> width_rp_func;
		typedef coo_func<DenseType, SparseType, ROW_WIDTH> width_coo_func;
		return compute_task::ptr(
				new block_spmm_task_impl<width_rp_func,  width_coo_func>(
					in_row_portions, *output, output_stream, io, mat, order));
	}
public:
	static task_creator::ptr create(const sparse_matrix &mat,
			size_t num_in_cols) {
		return task_creator::ptr(new spmm_creator<DenseType, SparseType>(
					mat, num_in_cols));
	}

	virtual bool is_complete() const {
		if (output_stream == NULL)
			return true;
		else
			return output_stream->is_complete();
	}

	virtual bool set_data(matrix_store::const_ptr in,
			matrix_store::ptr out);

	virtual std::vector<EM_object *> get_EM_objs() {
		std::vector<EM_object *> ret;
		if (!output->is_in_mem()) {
			const EM_object *obj = dynamic_cast<const EM_object *>(output.get());
			ret.push_back(const_cast<EM_object *>(obj));
		}
		return ret;
	}

	virtual compute_task::ptr create(const matrix_io &io) const {
		if (order) {
			switch (output->get_num_cols()) {
				case 1: return create_block_compute_task<1>(io);
				case 2: return create_block_compute_task<2>(io);
				case 4: return create_block_compute_task<4>(io);
				case 8: return create_block_compute_task<8>(io);
				case 16: return create_block_compute_task<16>(io);
				default: return create_block_compute_task<0>(io);
			}
		}
		else {
			assert(in_row_portions);
			switch (output->get_num_cols()) {
				case 1: return compute_task::ptr(new fg_row_spmm_task<DenseType,
								SparseType, 1>(in_row_portions, *output, io));
				case 2: return compute_task::ptr(new fg_row_spmm_task<DenseType,
								SparseType, 2>(in_row_portions, *output, io));
				case 4: return compute_task::ptr(new fg_row_spmm_task<DenseType,
								SparseType, 4>(in_row_portions, *output, io));
				case 8: return compute_task::ptr(new fg_row_spmm_task<DenseType,
								SparseType, 8>(in_row_portions, *output, io));
				case 16: return compute_task::ptr(new fg_row_spmm_task<DenseType,
								 SparseType, 16>(in_row_portions, *output, io));
				default: return compute_task::ptr(new fg_row_spmm_task<DenseType,
								 SparseType, 0>(in_row_portions, *output, io));
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
		= mem_matrix_store::CHUNK_SIZE / block_size.get_num_rows();
	return std::min(std::max(size, 1UL), max_size);
}

}

/*
 * This is a base class for a sparse matrix. It provides a set of functions
 * to perform computation on the sparse matrix. Currently, it has matrix
 * vector multiplication and matrix matrix multiplication.
 * We assume the sparse matrix is stored in external memory. If the matrix
 * is in memory, we can use in_mem_io to access the sparse matrix in memory
 * while reusing all the existing code for computation.
 */
class sparse_matrix: public detail::EM_object
{
	// Whether the matrix is represented by the FlashGraph format.
	bool is_fg;
	size_t nrows;
	size_t ncols;
	// The type of a non-zero entry.
	const scalar_type *entry_type;
	bool symmetric;
	detail::EM_object::io_set::ptr ios;

	template<class DenseType, class SparseType>
	detail::task_creator::ptr get_multiply_creator(size_t num_in_cols) const {
		return detail::spmm_creator<DenseType, SparseType>::create(*this,
				num_in_cols);
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

	void reset_ios() {
		ios = NULL;
	}

	void _transpose() {
		size_t tmp = nrows;
		nrows = ncols;
		ncols = tmp;
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

	std::vector<safs::io_interface::ptr> create_ios() const;

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
	void compute(detail::task_creator::ptr creator, size_t num_block_rows,
			size_t min_num_brows) const;
	/*
	 * This method defines how data in the matrix is accessed.
	 * `num_block_rows' defines how many block rows are accessed in an I/O.
	 */
	virtual matrix_io_generator::ptr create_io_gen(size_t num_block_rows,
			size_t min_num_brows) const = 0;

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
	virtual detail::block_exec_order::ptr get_multiply_order(
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

	virtual sparse_matrix::ptr transpose() const = 0;

	bool multiply(detail::matrix_store::const_ptr in,
			detail::matrix_store::ptr out, detail::task_creator::ptr create) const;

	/*
	 * This version of SpMM allows users to provide the output matrix.
	 * It requires users to initialize the output matrix.
	 */
	template<class DenseType, class SparseType>
	bool multiply(detail::matrix_store::const_ptr in,
			detail::matrix_store::ptr out) const {
		return multiply(in, out, get_multiply_creator<DenseType, SparseType>(
					in->get_num_cols()));
	}

	template<class DenseType, class SparseType>
	bool multiply(detail::vec_store::const_ptr in,
			detail::vec_store::ptr out) const {
		return multiply<DenseType, SparseType>(
				in->conv2mat(in->get_length(), 1, true),
				out->conv2mat(out->get_length(), 1, true));
	}
};

namespace detail
{

template<class DenseType, class SparseType>
spmm_creator<DenseType, SparseType>::spmm_creator(const sparse_matrix &_mat,
		size_t num_in_cols): mat(_mat)
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
		size_t num_cols = num_in_cols;
		num_cols = 1 << ((size_t) log2(num_cols));
		size_t sb_size = cal_super_block_size(mat.get_block_size(),
				sizeof(DenseType) * num_cols);
		order = mat.get_multiply_order(sb_size, sb_size);
	}
}

template<class DenseType, class SparseType>
bool spmm_creator<DenseType, SparseType>::set_data(
		matrix_store::const_ptr in, matrix_store::ptr out)
{
	if (in->get_type() != get_scalar_type<DenseType>()
			|| out->get_type() != get_scalar_type<DenseType>()) {
		BOOST_LOG_TRIVIAL(error) << "wrong matrix type in spmm creator";
		return false;
	}
	this->output = out;
	this->in_row_portions = row_portions::create(in);
	if (in_row_portions == NULL)
		return false;
	if (!output->is_in_mem())
		output_stream = EM_matrix_stream::create(output);
	return true;
}

}

void init_flash_matrix(config_map::ptr configs);
void destroy_flash_matrix();

}

#endif
