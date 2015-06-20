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
#include "safs_file.h"
#include "in_mem_io.h"

#include "sparse_matrix_format.h"
#include "matrix_config.h"
#include "local_vec_store.h"

namespace fm
{

void sparse_block_2d::verify(const block_2d_size &block_size) const
{
	size_t rel_row_id = 0;
	size_t num_rows = 0;
	rp_edge_iterator it = get_first_edge_iterator();
	while (!is_block_end(it)) {
		size_t num_nz = 0;
		while (it.has_next()) {
			num_nz++;
			it.next();
		}
		assert(num_nz <= block_size.get_num_cols());
		if (it.get_rel_row_idx() > 0)
			assert(rel_row_id < it.get_rel_row_idx());
		rel_row_id = it.get_rel_row_idx();
		num_rows++;
		it = get_next_edge_iterator(it);
	}
	assert(num_rows <= block_size.get_num_rows());
}

void sparse_block_2d::append(const sparse_row_part &part, size_t part_size)
{
	char *end = row_parts + rparts_size;

	memcpy(end, &part, part_size);
	assert(((size_t) rparts_size) + part_size
			<= std::numeric_limits<uint32_t>::max());
	rparts_size += part_size;
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

void SpM_2d_index::safs_dump(const std::string &file) const
{
	size_t size = get_size(get_num_entries());
	safs::safs_file f(safs::get_sys_RAID_conf(), file);
	bool ret = f.create_file(size);
	assert(ret);

	safs::file_io_factory::shared_ptr io_fac = safs::create_io_factory(
			file, safs::REMOTE_ACCESS);
	if (io_fac == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"can't create io factory for %1%") % file;
		return;
	}

	safs::io_interface::ptr io = create_io(io_fac, thread::get_curr_thread());
	if (io == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"can't create io instance for %1%") % file;
		return;
	}

	local_buf_vec_store buf(0, ROUNDUP(size, PAGE_SIZE),
			get_scalar_type<char>(), -1);
	memcpy(buf.get_raw_arr(), this, size);
	safs::data_loc_t loc(io_fac->get_file_id(), 0);
	safs::io_request req(buf.get_raw_arr(), loc, buf.get_length(), WRITE);
	io->access(&req, 1);
	io->wait4complete(1);
}

off_t SpM_2d_index::get_block_row_off(size_t idx) const
{
	if (idx >= get_num_entries()) {
		BOOST_LOG_TRIVIAL(error) << "get_block_row_off: out of range.";
		return -1;
	}
	return offs[idx];
}

SpM_2d_index::ptr SpM_2d_index::safs_load(const std::string &idx_file)
{
	safs::file_io_factory::shared_ptr io_fac = safs::create_io_factory(
			idx_file, safs::GLOBAL_CACHE_ACCESS);
	if (io_fac == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"can't create io factory for %1%") % idx_file;
		return SpM_2d_index::ptr();
	}

	safs::io_interface::ptr io = create_io(io_fac, thread::get_curr_thread());
	if (io == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"can't create io instance for %1%") % idx_file;
		return SpM_2d_index::ptr();
	}

	size_t size = safs::safs_file(safs::get_sys_RAID_conf(),
			idx_file).get_size();
	char *data = (char *) malloc(size);
	io->access(data, 0, size, READ);
	SpM_2d_index *index = (SpM_2d_index *) data;
	index->header.verify();
	return ptr(index, deleter());
}

SpM_2d_index::ptr SpM_2d_index::load(const std::string &idx_file)
{
	FILE *f = fopen(idx_file.c_str(), "r");
	if (f == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format("can't open %1%: %2%")
			% idx_file % strerror(errno);
		return SpM_2d_index::ptr();
	}

	size_t size = safs::native_file(idx_file).get_size();
	char *data = (char *) malloc(size);
	size_t ret = fread(data, size, 1, f);
	if (ret == 0) {
		BOOST_LOG_TRIVIAL(error) << boost::format("can't read %1%: %2%")
			% idx_file % strerror(errno);
		fclose(f);
		free(data);
		return SpM_2d_index::ptr();
	}

	fclose(f);
	SpM_2d_index *index = (SpM_2d_index *) data;
	index->header.verify();
	return ptr(index, deleter());
}

void SpM_2d_storage::verify() const
{
	matrix_header *header = (matrix_header *) data.get();
	header->verify();
	block_2d_size block_size = index->get_header().get_2d_block_size();
#pragma omp parallel for
	for (size_t i = 0; i < get_num_block_rows(); i++) {
		block_row_iterator brow_it = get_block_row_it(i);
		long prev_block_col_idx = -1;
		size_t num_blocks = 0;
		while (brow_it.has_next()) {
			const sparse_block_2d &block = brow_it.next();
			assert(block.get_block_row_idx() == i);
			assert((long) block.get_block_col_idx() > prev_block_col_idx);
			prev_block_col_idx = block.get_block_col_idx();
			block.verify(block_size);
			num_blocks++;
		}
	}
}

SpM_2d_storage::ptr SpM_2d_storage::safs_load(const std::string &mat_file,
		SpM_2d_index::ptr index)
{
	size_t size = safs::safs_file(safs::get_sys_RAID_conf(), mat_file).get_size();
	char *data = NULL;
	int mret = posix_memalign((void **) &data, PAGE_SIZE, size);
	BOOST_VERIFY(mret == 0);

	safs::file_io_factory::shared_ptr io_fac = safs::create_io_factory(
			mat_file, safs::GLOBAL_CACHE_ACCESS);
	if (io_fac == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"can't create io factory for %1%") % mat_file;
		return SpM_2d_storage::ptr();
	}
	safs::io_interface::ptr io = create_io(io_fac, thread::get_curr_thread());
	if (io == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"can't create io instance for %1%") % mat_file;
		return SpM_2d_storage::ptr();
	}

	io->access(data, 0, size, READ);
	matrix_header *header = (matrix_header *) data;
	header->verify();
	return ptr(new SpM_2d_storage(std::shared_ptr<char>(data, deleter()),
				index, mat_file));
}

SpM_2d_storage::ptr SpM_2d_storage::load(const std::string &mat_file,
			SpM_2d_index::ptr index)
{
	size_t size = safs::native_file(mat_file).get_size();
	char *data = NULL;
	int mret = posix_memalign((void **) &data, PAGE_SIZE, size);
	BOOST_VERIFY(mret == 0);
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
	matrix_header *header = (matrix_header *) data;
	header->verify();
	return ptr(new SpM_2d_storage(std::shared_ptr<char>(data, deleter()),
				index, mat_file));
}

SpM_2d_storage::ptr SpM_2d_storage::create(const matrix_header &header,
		const vector_vector &vv, SpM_2d_index::ptr index)
{
	if (!vv.is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The vector of vectors has to be in memory";
		return SpM_2d_storage::ptr();
	}
	vector::ptr vec = vv.cat();
	if (vec->get_type().get_type() != get_type<char>()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The vector of vectors contains a wrong type of data";
		return SpM_2d_storage::ptr();
	}

	size_t size = sizeof(header) + vec->get_length();
	// The sparse matrix multiplication accesses data in pages. We have to
	// make sure the array that stores the sparse matrix is aligned to
	// page size.
	char *data = NULL;
	int ret = posix_memalign((void **) &data, PAGE_SIZE,
			ROUNDUP(size, PAGE_SIZE));
	BOOST_VERIFY(ret == 0);
	*(matrix_header *) data = header;
	memcpy(data + sizeof(header),
			dynamic_cast<const detail::mem_vec_store &>(vec->get_data()).get_raw_arr(),
			vec->get_length());
	return ptr(new SpM_2d_storage(std::shared_ptr<char>(data, deleter()),
				index, "anonymous"));
}

safs::file_io_factory::shared_ptr SpM_2d_storage::create_io_factory() const
{
	return safs::file_io_factory::shared_ptr(new safs::in_mem_io_factory(
				data, mat_file_id, mat_name));
}

}
