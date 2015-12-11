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

#include <pthread.h>

#include "sparse_matrix.h"
#include "matrix_io.h"
#include "matrix_config.h"

namespace fm
{

/*
 * This represents a collection of row blocks that we usually access with
 * a single I/O request. However, due to load balancing, we sometimes
 * need to access the collection of row blocks with many smaller I/O
 * requests. In this case, each row block is accessed by one I/O request.
 */
class large_row_io
{
	// A reference to the row block vector of the graph.
	const std::vector<row_block> *blocks;
	// The offset of the first row block in the row block vector.
	off_t first_row_block;
	size_t num_row_blocks;
	size_t tot_num_rows;
public:
	large_row_io(const std::vector<row_block> &_blocks, off_t first_row_block,
			size_t num_row_blocks, size_t tot_num_rows) {
		this->blocks = &_blocks;
		this->first_row_block = first_row_block;
		this->num_row_blocks = num_row_blocks;
		this->tot_num_rows = tot_num_rows;
	}

	/*
	 * Create an I/O request that access all row blocks in the collection.
	 */
	matrix_io get_io(size_t tot_num_cols, int file_id) {
		off_t first_row_id = first_row_block * matrix_conf.get_row_block_size();
		matrix_loc top_left(first_row_id, 0);
		size_t num_rows = tot_num_rows;
		off_t first_row_offset = blocks->at(first_row_block).get_offset();
		size_t size = blocks->at(first_row_block + num_row_blocks).get_offset()
			- first_row_offset;
		matrix_io ret(top_left, num_rows, tot_num_cols,
					safs::data_loc_t(file_id, first_row_offset), size);
		first_row_block += num_row_blocks;
		num_row_blocks = 0;
		tot_num_rows = 0;
		return ret;
	}

	/*
	 * Create an I/O request that accesses a single row block.
	 */
	matrix_io get_sub_io(size_t tot_num_cols, int file_id) {
		off_t first_row_id = first_row_block * matrix_conf.get_row_block_size();
		matrix_loc top_left(first_row_id, 0);
		size_t num_curr_row_blocks = std::min(num_row_blocks,
				(size_t) matrix_conf.get_rb_steal_io_size());
		size_t num_rows = std::min(tot_num_rows,
				(size_t) num_curr_row_blocks * matrix_conf.get_row_block_size());
		off_t first_row_offset = blocks->at(first_row_block).get_offset();
		size_t size = blocks->at(first_row_block + num_curr_row_blocks).get_offset()
			- first_row_offset;
		matrix_io ret(top_left, num_rows,
				tot_num_cols, safs::data_loc_t(file_id, first_row_offset), size);
		first_row_block += num_curr_row_blocks;
		num_row_blocks -= num_curr_row_blocks;
		tot_num_rows -= num_rows;
		return ret;
	}

	bool has_data() const {
		return num_row_blocks > 0;
	}
};

/*
 * The I/O generator that access a matrix on disks by rows.
 * An I/O generator are assigned a number of row blocks that it can accesses.
 * Each thread has an I/O generator and gets I/O requests from its own I/O
 * generator. When load balancing kicks in, a thread will try to steal I/O requests
 * from other threads' I/O generators.
 */
class row_io_generator: public matrix_io_generator
{
	std::vector<large_row_io> ios;
	// The current offset in the row_block vector.
	volatile off_t curr_io_off;
	int file_id;
	size_t tot_num_cols;

	pthread_spinlock_t lock;
public:
	row_io_generator(const std::vector<row_block> &_blocks, size_t tot_num_rows,
			size_t tot_num_cols, int file_id, const row_block_mapper &mapper);

	/*
	 * This method is called by the worker thread that owns the I/O generator
	 * to get I/O requests.
	 */
	virtual matrix_io get_next_io() {
		// TODO it should return a smaller I/O to improve load balancing at
		// the end of SpMM.
		matrix_io ret;
		pthread_spin_lock(&lock);
		// It's possible that all IOs have been stolen.
		// We have to check it.
		if ((size_t) curr_io_off < ios.size()) {
			assert(ios[curr_io_off].has_data());
			ret = ios[curr_io_off++].get_io(tot_num_cols, file_id);
			assert(ret.is_valid());
		}
		pthread_spin_unlock(&lock);
		return ret;
	}

	virtual bool has_next_io() {
		// The last entry in the row_block vector is the size of the matrix
		// file.
		return (size_t) curr_io_off < ios.size();
	}
};

row_io_generator::row_io_generator(const std::vector<row_block> &blocks,
		size_t tot_num_rows, size_t tot_num_cols, int file_id,
		const row_block_mapper &mapper)
{
	pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	// blocks[blocks.size() - 1] is an empty block. It only indicates
	// the end of the matrix file.
	// blocks[blocks.size() - 2] is the last row block. It's possible that
	// it's smaller than the full-size row block.
	for (size_t i = 0; i < mapper.get_num_ranges(); i++) {
		size_t num_row_blocks = mapper.get_range(i).num;
		off_t rb_off = mapper.get_range(i).idx;
		size_t num_rows = std::min(num_row_blocks * matrix_conf.get_row_block_size(),
				tot_num_rows - rb_off * matrix_conf.get_row_block_size());
		ios.emplace_back(blocks, rb_off, num_row_blocks, num_rows);
	}
	curr_io_off = 0;
	this->tot_num_cols = tot_num_cols;
	this->file_id = file_id;
}

/*
 * The I/O generator for accessing 2D-partitioned blocks.
 * A generated I/O request may access multiple 2D-partitioned blocks.
 * Each I/O generator still accesses all blocks in a block row.
 */
class b2d_io_generator: public matrix_io_generator
{
	pthread_spinlock_t lock;
	size_t brow_idx;
	size_t min_num_brows;
	size_t io_num_brows;
	// This is the starting point where we give a single block row for
	// each assignment for load balancing. This points to a block row.
	size_t balance_point;

	SpM_2d_index::ptr idx;
	block_2d_size block_size;
	int file_id;
	size_t tot_num_cols;
public:
	b2d_io_generator(SpM_2d_index::ptr idx, int file_id, size_t io_num_brows,
			size_t min_num_brows);

	virtual matrix_io get_next_io();

	virtual bool has_next_io() {
		return brow_idx < idx->get_num_block_rows();
	}
};

b2d_io_generator::b2d_io_generator(SpM_2d_index::ptr idx, int file_id,
		size_t io_num_brows, size_t min_num_brows): block_size(
			idx->get_header().get_2d_block_size())
{
	this->idx = idx;
	this->io_num_brows = io_num_brows;
	this->min_num_brows = min_num_brows;
	if (io_num_brows > min_num_brows) {
		size_t num_threads = detail::mem_thread_pool::get_global_num_threads();
		// We reserve some block rows in the matrix that are assigned one blow
		// row at a time.
		if (idx->get_num_block_rows() > io_num_brows * num_threads)
			balance_point = idx->get_num_block_rows() - num_threads;
		else
			balance_point = 0;
	}
	else
		balance_point = idx->get_num_block_rows();
	this->tot_num_cols = idx->get_header().get_num_cols();
	this->file_id = file_id;
	this->brow_idx = 0;
	pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
}

matrix_io b2d_io_generator::get_next_io()
{
	matrix_io ret;
	pthread_spin_lock(&lock);
	// It's possible that all IOs have been stolen.
	// We have to check it.
	if (has_next_io()) {
		// If we have passed the balancing point, we now assign one block row
		// to a thread each time.
		if (brow_idx >= balance_point)
			io_num_brows = min_num_brows;

		matrix_loc mat_loc(brow_idx * block_size.get_num_rows(), 0);
		safs::data_loc_t data_loc(file_id, idx->get_block_row_off(brow_idx));
		size_t num_brows = std::min(io_num_brows,
				idx->get_num_block_rows() - brow_idx);
		size_t brow_size = idx->get_block_row_off(brow_idx
				+ num_brows) - idx->get_block_row_off(brow_idx);
		ret = matrix_io(mat_loc, block_size.get_num_rows() * num_brows,
				tot_num_cols, data_loc, brow_size);
		brow_idx += num_brows;
	}
	pthread_spin_unlock(&lock);
	return ret;
}

matrix_io_generator::ptr matrix_io_generator::create(
		const std::vector<row_block> &_blocks, size_t tot_num_rows,
		size_t tot_num_cols, int file_id, const row_block_mapper &mapper)
{
	return matrix_io_generator::ptr(new row_io_generator(_blocks,
				tot_num_rows, tot_num_cols, file_id, mapper));
}

matrix_io_generator::ptr matrix_io_generator::create(SpM_2d_index::ptr idx,
		int file_id, size_t io_num_brows, size_t min_num_brows)
{
	return matrix_io_generator::ptr(new b2d_io_generator(idx, file_id,
				io_num_brows, min_num_brows));
}

void row_block_mapper::init(size_t num_rbs, size_t range_size)
{
	// We now need to jump to the next row block that belongs to
	// the current I/O generator.
	for (size_t i = 0; i < num_rbs; i += range_size) {
		struct rb_range range;
		range.idx = i;
		range.num = std::min(range_size, num_rbs - i);
		ranges.push_back(range);
	}
}

row_block_mapper::row_block_mapper(const SpM_2d_index &index, size_t range_size)
{
	init(index.get_num_block_rows(), range_size);
}

}
