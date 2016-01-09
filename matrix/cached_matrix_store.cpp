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

#include <boost/format.hpp>

#include "cached_matrix_store.h"
#include "combined_matrix_store.h"
#include "EM_dense_matrix.h"
#include "local_matrix_store.h"

namespace fm
{

namespace detail
{

static std::vector<std::weak_ptr<cached_matrix_store> > cached_mats;

cached_matrix_store::ptr cached_matrix_store::create(size_t num_rows,
		size_t num_cols, int num_nodes, const scalar_type &type,
		size_t num_cached_vecs, matrix_layout_t cached_layout)
{
	cached_matrix_store::ptr ret(new cached_matrix_store(num_rows,
				num_cols, num_nodes, type, num_cached_vecs, cached_layout));
	cached_mats.push_back(ret);
	return ret;
}

void cached_matrix_store::drop_all_cache()
{
	for (size_t i = 0; i < cached_mats.size(); i++) {
		cached_matrix_store::ptr store = cached_mats[i].lock();
		if (store)
			store->drop_cache();
	}
	cached_mats.clear();
}

static inline matrix_layout_t decide_layout(size_t num_rows, size_t num_cols)
{
	if (num_rows < num_cols)
		return matrix_layout_t::L_ROW;
	else
		return matrix_layout_t::L_COL;
}

// This is created for transpose.
cached_matrix_store::cached_matrix_store(size_t num_rows, size_t num_cols,
		const scalar_type &type): matrix_store(num_rows, num_cols, false, type)
{
}

cached_matrix_store::cached_matrix_store(size_t num_rows, size_t num_cols,
		int num_nodes, const scalar_type &type, size_t num_cached_vecs,
		matrix_layout_t cached_layout): matrix_store(num_rows, num_cols, false, type)
{
	assert(num_cached_vecs > 0);
	em_buf = EM_matrix_store::create(num_rows, num_cols,
			decide_layout(num_rows, num_cols), type);
	matrix_store::const_ptr uncached;
	if (em_buf->is_wide()) {
		num_cached_vecs = std::min(num_cached_vecs, num_rows);
		cached_buf = mem_matrix_store::create(num_cached_vecs, num_cols,
				cached_layout, type, num_nodes);
		std::vector<off_t> row_offs(num_rows - num_cached_vecs);
		for (size_t i = 0; i < row_offs.size(); i++)
			row_offs[i] = num_cached_vecs + i;
		if (!row_offs.empty()) {
			uncached = em_buf->get_rows(row_offs);
			assert(uncached);
		}
	}
	else {
		num_cached_vecs = std::min(num_cached_vecs, num_cols);
		cached_buf = mem_matrix_store::create(num_rows, num_cached_vecs,
				cached_layout, type, num_nodes);
		std::vector<off_t> col_offs(num_cols - num_cached_vecs);
		for (size_t i = 0; i < col_offs.size(); i++)
			col_offs[i] = num_cached_vecs + i;
		if (!col_offs.empty()) {
			uncached = em_buf->get_cols(col_offs);
			assert(uncached);
		}
	}
	// If the entire matrix is cached.
	if (uncached == NULL)
		mixed = cached_buf;
	else {
		std::vector<matrix_store::const_ptr> mats(2);
		mats[0] = cached_buf;
		mats[1] = uncached;
		mixed = combined_matrix_store::create(mats, cached_buf->store_layout());
	}
	cached = cached_buf;
	em_store = em_buf;
}

std::string cached_matrix_store::get_name() const
{
	return boost::str(boost::format("cached%1%_mat(%2%)")
			% get_num_cached_vecs() % em_store->get_name());
}

matrix_store::const_ptr cached_matrix_store::transpose() const
{
	cached_matrix_store *store = new cached_matrix_store(get_num_cols(),
			get_num_rows(), get_type());
	store->mixed = mixed->transpose();
	if (cached)
		store->cached = cached->transpose();
	store->em_store = EM_matrix_store::cast(em_store->transpose());
	// We don't need to set em_buf and cached_buf because these are only
	// required for write data to the matrix store.
	cached_matrix_store::ptr ret(store);
	if (cached)
		cached_mats.push_back(ret);
	return ret;
}

void cached_matrix_store::write_portion_async(
		local_matrix_store::const_ptr portion, off_t start_row, off_t start_col)
{
	assert(em_buf);
	assert(cached_buf);

	if (is_wide()) {
		local_matrix_store::const_ptr sub_portion = portion->get_portion(0, 0,
				cached_buf->get_num_rows(), portion->get_num_cols());
		local_matrix_store &mutable_sub
			= const_cast<local_matrix_store &>(*sub_portion);
		size_t mem_portion_size = cached_buf->get_portion_size().second;
		for (size_t lstart_col = 0; lstart_col < portion->get_num_cols();
				lstart_col += mem_portion_size) {
			size_t lnum_cols = std::min(mem_portion_size,
					portion->get_num_cols() - lstart_col);
			mutable_sub.resize(0, lstart_col, sub_portion->get_num_rows(),
					lnum_cols);
			cached_buf->write_portion_async(sub_portion, start_row,
					start_col + lstart_col);
		}
	}
	else {
		local_matrix_store::const_ptr sub_portion = portion->get_portion(0, 0,
				portion->get_num_rows(), cached_buf->get_num_cols());
		local_matrix_store &mutable_sub
			= const_cast<local_matrix_store &>(*sub_portion);
		size_t mem_portion_size = cached_buf->get_portion_size().first;
		for (size_t lstart_row = 0; lstart_row < portion->get_num_rows();
				lstart_row += mem_portion_size) {
			size_t lnum_rows = std::min(mem_portion_size,
					portion->get_num_rows() - lstart_row);
			mutable_sub.resize(lstart_row, 0, lnum_rows,
					sub_portion->get_num_cols());
			cached_buf->write_portion_async(sub_portion, start_row + lstart_row,
					start_col);
		}
	}
	em_buf->write_portion_async(portion, start_row, start_col);
}

matrix_store::const_ptr cached_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	if (!is_wide()) {
		BOOST_LOG_TRIVIAL(error) << "can't get rows from a tall matrix\n";
		return matrix_store::const_ptr();
	}

	// If all rows are cached.
	size_t num_cached = get_num_cached_vecs();
	bool all_cached = true;
	for (size_t i = 0; i < idxs.size(); i++)
		if ((size_t) idxs[i] >= num_cached)
			all_cached = false;
	// If we want to access all rows of the matrix, we can return
	// the entire matrix directly.
	if (all_cached && idxs.size() == cached->get_num_rows()
			&& std::is_sorted(idxs.begin(), idxs.end()))
		return cached;
	else if (all_cached)
		return cached->get_rows(idxs);

	return mixed->get_rows(idxs);
}

void cached_matrix_store::set_cache_portion(bool cache_portion)
{
	const_cast<matrix_store &>(*mixed).set_cache_portion(cache_portion);
}

}

}
