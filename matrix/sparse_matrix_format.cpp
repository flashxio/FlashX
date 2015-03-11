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

#include <boost/format.hpp>

#include "log.h"
#include "native_file.h"
#include "in_mem_io.h"

#include "sparse_matrix_format.h"
#include "matrix_config.h"
#include "factor.h"
#include "generic_type.h"

namespace fm
{

void sparse_block_2d::verify(const block_2d_size &block_size) const
{
	row_part_iterator it = get_iterator();
	size_t rel_row_id = 0;
	size_t num_rows = 0;
	while (it.has_next()) {
		const sparse_row_part &part = it.next();
		assert(part.get_num_non_zeros() <= block_size.get_num_cols());
		if (part.get_rel_row_idx() > 0)
			assert(rel_row_id < part.get_rel_row_idx());
		rel_row_id = part.get_rel_row_idx();
		num_rows++;
	}
	assert(num_rows <= block_size.get_num_rows());
}

void sparse_block_2d::append(const sparse_row_part &part)
{
	char *end = row_parts + rparts_size;
	memcpy(end, &part, part.get_size());
	assert(((size_t) rparts_size) + part.get_size()
			<= std::numeric_limits<uint32_t>::max());
	rparts_size += part.get_size();
}

SpM_2d_index::ptr SpM_2d_index::create(const matrix_header &header,
		const std::vector<off_t> &offs)
{
	block_2d_size block_size = header.get_2d_block_size();
	if (offs.size() != block_size.cal_num_block_rows(header.get_num_rows()) + 1) {
		BOOST_LOG_TRIVIAL(error) << "There are an incorrect number of offsets";
		return SpM_2d_index::ptr();
	}
	void *buf = NULL;
	int ret = posix_memalign(&buf, PAGE_SIZE, get_size(offs.size()));
	if (ret) {
		BOOST_LOG_TRIVIAL(error) << "Can't allocate memory";
		return SpM_2d_index::ptr();
	}
	SpM_2d_index *idx = new (buf) SpM_2d_index(header);
	for (size_t i = 0; i < offs.size(); i++) 
		idx->offs[i] = offs[i];
	return SpM_2d_index::ptr(idx, deleter());
}

size_t SpM_2d_index::get_num_block_rows() const
{
	block_2d_size block_size = header.get_2d_block_size();
	return block_size.cal_num_block_rows(header.get_num_rows());
}

void SpM_2d_index::dump(const std::string &file) const
{
	FILE *f = fopen(file.c_str(), "w");
	if (f == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format("can't open %1%: %2%")
			% file % strerror(errno);
		return;
	}

	size_t ret = fwrite(this, get_size(get_num_entries()), 1, f);
	if (ret == 0) {
		BOOST_LOG_TRIVIAL(error) << boost::format("can't write to %1%: %2%")
			% file % strerror(errno);
		fclose(f);
		return;
	}
	fclose(f);
}

off_t SpM_2d_index::get_block_row_off(size_t idx) const
{
	if (idx >= get_num_entries()) {
		BOOST_LOG_TRIVIAL(error) << "get_block_row_off: out of range.";
		return -1;
	}
	return offs[idx];
}

SpM_2d_index::ptr SpM_2d_index::load(const std::string &idx_file)
{
	size_t size = safs::native_file(idx_file).get_size();
	char *data = (char *) malloc(size);
	FILE *f = fopen(idx_file.c_str(), "r");
	if (f == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format("can't open %1%: %2%")
			% idx_file % strerror(errno);
		return SpM_2d_index::ptr();
	}
	size_t ret = fread(data, size, 1, f);
	if (ret == 0) {
		BOOST_LOG_TRIVIAL(error) << boost::format("can't read %1%: %2%")
			% idx_file % strerror(errno);
		fclose(f);
		return SpM_2d_index::ptr();
	}

	fclose(f);
	return ptr((SpM_2d_index *) data, deleter());
}

void SpM_2d_storage::verify() const
{
	block_2d_size block_size = index->get_header().get_2d_block_size();
	for (size_t i = 0; i < get_num_block_rows(); i++) {
		block_row_iterator brow_it = get_block_row_it(i);
		while (brow_it.has_next()) {
			const sparse_block_2d &block = brow_it.next();
			block.verify(block_size);
		}
	}
}

SpM_2d_storage::ptr SpM_2d_storage::load(const std::string &mat_file,
			SpM_2d_index::ptr index)
{
	size_t size = safs::native_file(mat_file).get_size();
	char *data = new char[size];
	FILE *f = fopen(mat_file.c_str(), "r");
	if (f == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format("can't open %1%: %2%")
			% mat_file % strerror(errno);
		return SpM_2d_storage::ptr();
	}
	size_t ret = fread(data, size, 1, f);
	if (ret == 0) {
		BOOST_LOG_TRIVIAL(error) << boost::format("can't read %1%: %2%")
			% mat_file % strerror(errno);
		fclose(f);
		return SpM_2d_storage::ptr();
	}
	fclose(f);
	return ptr(new SpM_2d_storage(std::unique_ptr<char[]>(data), index,
				mat_file));
}

SpM_2d_storage::ptr SpM_2d_storage::create(const matrix_header &header,
		const vector_vector &vv, SpM_2d_index::ptr index)
{
	if (!vv.is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The vector of vectors has to be in memory";
		return SpM_2d_storage::ptr();
	}
	mem_vector::ptr vec = mem_vector::cast(vv.cat());
	if (vec->get_type().get_type() == get_type<char>()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The vector of vectors contains a wrong type of data";
		return SpM_2d_storage::ptr();
	}

	size_t size = sizeof(header) + vec->get_length();
	char *data = new char[size];
	*(matrix_header *) data = header;
	memcpy(data + sizeof(header), vec->get_raw_arr(), vec->get_length());
	return ptr(new SpM_2d_storage(std::unique_ptr<char[]>(data), index,
				"anonymous"));
}

safs::file_io_factory::shared_ptr SpM_2d_storage::create_io_factory() const
{
	return safs::file_io_factory::shared_ptr(new safs::in_mem_io_factory(
				data.get(), mat_file_id, mat_name));
}

////////// The following code constructs a 2D-partitioned sparse matrix/////////

class set_2d_label_operate: public type_set_operate<factor_value_t>
{
	block_2d_size block_size;
public:
	set_2d_label_operate(const block_2d_size &_size): block_size(_size) {
	}

	virtual void set(factor_value_t *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		assert(col_idx == 0);
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = (row_idx + i) / block_size.get_num_rows();
	}
};

class part_2d_apply_operate: public gr_apply_operate<sub_vector_vector>
{
	size_t row_len;
	block_2d_size block_size;
public:
	part_2d_apply_operate(const block_2d_size &_size,
			size_t row_len): block_size(_size) {
		this->row_len = row_len;
	}

	void run(const void *key, const sub_vector_vector &val,
			mem_vector &out) const;

	const scalar_type &get_key_type() const {
		return get_scalar_type<factor_value_t>();
	}

	const scalar_type &get_output_type() const {
		return get_scalar_type<char>();
	}

	size_t get_num_out_eles() const {
		return 1;
	}
};

void part_2d_apply_operate::run(const void *key, const sub_vector_vector &val,
		mem_vector &out) const
{
	size_t block_height = block_size.get_num_rows();
	size_t block_width = block_size.get_num_cols();
	size_t num_blocks = ceil(((double) row_len) / block_width);
	factor_value_t block_row_id = *(const factor_value_t *) key;
	printf("block row id: %d, #blocks: %ld\n", block_row_id, num_blocks);
	size_t tot_num_non_zeros = 0;
	size_t max_row_parts = 0;
	for (size_t i = 0; i < val.get_num_vecs(); i++) {
		const fg::ext_mem_undirected_vertex *v
			= (const fg::ext_mem_undirected_vertex *) val.get_raw_arr(i);
		assert(v->get_id() / block_height == (size_t) block_row_id);
		tot_num_non_zeros += v->get_num_edges();
		// I definitely over estimate the number of row parts.
		// If a row doesn't have many non-zero entries, I assume that
		// the non-zero entries distribute evenly across all row parts.
		max_row_parts += std::min(num_blocks, v->get_num_edges());
	}

	std::vector<size_t> neigh_idxs(val.get_num_vecs());
	// The maximal size of a block.
	size_t max_block_size
		// Even if a block is empty, its header still exists. The size is
		// accurate.
		= sizeof(sparse_block_2d) * num_blocks
		// The size for row part headers is highly over estimated.
		+ sizeof(sparse_row_part) * max_row_parts
		// The size is accurate.
		+ sparse_row_part::get_col_entry_size() * tot_num_non_zeros;
	out.resize(max_block_size);
	size_t curr_size = 0;
	// The maximal size of a row part.
	size_t max_row_size = sparse_row_part::get_size(block_width);
	std::unique_ptr<char[]> buf = std::unique_ptr<char[]>(new char[max_row_size]);
	size_t num_non_zeros = 0;
	for (size_t col_idx = 0; col_idx < row_len; col_idx += block_width) {
		sparse_block_2d *block
			= new (out.get_raw_arr() + curr_size) sparse_block_2d(
					block_row_id, col_idx / block_width);
		for (size_t row_idx = 0; row_idx < val.get_num_vecs(); row_idx++) {
			const fg::ext_mem_undirected_vertex *v
				= (const fg::ext_mem_undirected_vertex *) val.get_raw_arr(row_idx);
			// If the vertex has no more edges left.
			if (neigh_idxs[row_idx] >= v->get_num_edges())
				continue;
			assert(v->get_neighbor(neigh_idxs[row_idx]) >= col_idx);
			// If the vertex has no edges that fall in the range.
			if (v->get_neighbor(neigh_idxs[row_idx]) >= col_idx + block_width)
				continue;

			sparse_row_part *part = new (buf.get()) sparse_row_part(row_idx);
			size_t idx = neigh_idxs[row_idx];
			for (; idx < v->get_num_edges()
					&& v->get_neighbor(idx) < col_idx + block_width; idx++)
				part->add(block_size, v->get_neighbor(idx));
			assert(part->get_size() <= max_row_size);
			neigh_idxs[row_idx] = idx;
			num_non_zeros += part->get_num_non_zeros();
			assert(block->get_size() + part->get_size()
					<= max_block_size - curr_size);
			block->append(*part);
		}
		// Only the non-empty blocks exist in a block row.
		if (!block->is_empty()) {
			curr_size += block->get_size();
			block->verify(block_size);
		}
	}
	out.resize(curr_size);
}

std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr> create_2d_matrix(
		vector_vector::ptr adjs, const block_2d_size &block_size)
{
	size_t num_rows = adjs->get_num_vecs();
	factor f(ceil(((double) num_rows) / block_size.get_num_rows()));
	factor_vector::ptr labels = factor_vector::create(f, num_rows);
	labels->set_data(set_2d_label_operate(block_size));
	vector_vector::ptr res = adjs->groupby(*labels,
			part_2d_apply_operate(block_size, num_rows));

	matrix_header mheader(matrix_type::SPARSE, 0, num_rows, num_rows,
			matrix_layout_t::L_ROW_2D, prim_type::P_BOOL, block_size);

	// Construct the index file of the adjacency matrix.
	std::vector<off_t> offsets(res->get_num_vecs() + 1);
	off_t off = sizeof(mheader);
	for (size_t i = 0; i < res->get_num_vecs(); i++) {
		offsets[i] = off;
		off += res->get_length(i);
	}
	offsets[res->get_num_vecs()] = off;
	SpM_2d_index::ptr idx = SpM_2d_index::create(mheader, offsets);

	return std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr>(
			idx, SpM_2d_storage::create(mheader, *res, idx));
}

void export_2d_matrix(vector_vector::ptr adjs, const block_2d_size &block_size,
		const std::string &mat_file, const std::string &mat_idx_file)
{
	size_t num_rows = adjs->get_num_vecs();
	factor f(ceil(((double) num_rows) / block_size.get_num_rows()));
	factor_vector::ptr labels = factor_vector::create(f, num_rows);
	labels->set_data(set_2d_label_operate(block_size));
	vector_vector::ptr res = adjs->groupby(*labels,
			part_2d_apply_operate(block_size, num_rows));

	matrix_header mheader(matrix_type::SPARSE, 0, num_rows, num_rows,
			matrix_layout_t::L_ROW_2D, prim_type::P_BOOL, block_size);
	FILE *f_2d = fopen(mat_file.c_str(), "w");
	if (f_2d == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format("open %1%: %2%")
			% mat_file % strerror(errno);
		return;
	}
	fwrite(&mheader, sizeof(mheader), 1, f_2d);
	bool ret = mem_vector::cast(res->cat())->export2(f_2d);
	assert(ret);
	fclose(f_2d);

	// Construct the index file of the adjacency matrix.
	std::vector<off_t> offsets(res->get_num_vecs() + 1);
	off_t off = sizeof(mheader);
	for (size_t i = 0; i < res->get_num_vecs(); i++) {
		offsets[i] = off;
		off += res->get_length(i);
	}
	offsets[res->get_num_vecs()] = off;
	SpM_2d_index::ptr mindex = SpM_2d_index::create(mheader, offsets);
	mindex->dump(mat_idx_file);
}

}
