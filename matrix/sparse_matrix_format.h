#ifndef __SPARSE_MATRIX_FORMAT_H__
#define __SPARSE_MATRIX_FORMAT_H__

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

#include "vertex.h"
#include "matrix_header.h"
#include "vector_vector.h"

namespace safs
{
	class file_io_factory;
}

namespace fm
{

/*
 * This iterates the edges in a row part.
 */
class rp_edge_iterator
{
	uint16_t rel_row_idx;
	// This always points to the beginning of the row part.
	uint16_t *rel_col_idx_start;
	// This points to the current location of the iterator on the row part.
	uint16_t *rel_col_idx_p;
	// This points to the non-zero value in the current location.
	const char *data_p;
	size_t entry_size;
public:
	rp_edge_iterator() {
		rel_row_idx = 0;
		rel_col_idx_start = NULL;
		rel_col_idx_p = NULL;
		data_p = NULL;
		entry_size = 0;
	}

	rp_edge_iterator(uint16_t rel_row_idx, uint16_t *rel_col_idx_start) {
		this->rel_row_idx = rel_row_idx;
		this->rel_col_idx_start = rel_col_idx_start;
		this->rel_col_idx_p = rel_col_idx_start;
		this->data_p = NULL;
		this->entry_size = 0;
	}

	rp_edge_iterator(uint16_t rel_row_idx, uint16_t *rel_col_idx_start,
			const char *data_p, size_t entry_size) {
		this->rel_row_idx = rel_row_idx;
		this->rel_col_idx_start = rel_col_idx_start;
		this->rel_col_idx_p = rel_col_idx_start;
		this->data_p = data_p;
		this->entry_size = entry_size;
	}

	bool is_valid() const {
		return rel_col_idx_start != NULL;
	}

	size_t get_rel_row_idx() const {
		return rel_row_idx;
	}

	off_t get_offset() const {
		return rel_col_idx_p - rel_col_idx_start;
	}

	bool has_next() const {
		// The highest bit of the relative col idx has to be 0.
		return *rel_col_idx_p <= (size_t) std::numeric_limits<int16_t>::max();
	};

	// This returns the relative column index.
	size_t get_curr() const {
		return *rel_col_idx_p;
	}

	template<class T>
	T get_curr_data() const {
		assert(data_p && entry_size == sizeof(T));
		return *(const T *) data_p;
	}

	// This returns the relative column index.
	size_t next() {
		size_t ret = *rel_col_idx_p;
		rel_col_idx_p++;
		data_p += entry_size;
		return ret;
	}

	void append(const block_2d_size &block_size, size_t col_idx) {
		*rel_col_idx_p = col_idx & block_size.get_ncol_mask();
		assert(*rel_col_idx_p <= (size_t) std::numeric_limits<int16_t>::max());
		rel_col_idx_p++;
	}

	const char *get_curr_addr() const {
		return (const char *) rel_col_idx_p;
	}

	const char *get_curr_data() const {
		return data_p;
	}
};

typedef std::pair<size_t, size_t> coo_nz_t;

/*
 * This stores the header of a row part. It doesn't contain any attributes.
 */
class sparse_row_part
{
	static const int num_bits = sizeof(uint16_t) * 8;
	uint16_t rel_row_idx;
	uint16_t rel_col_idxs[0];

	// The row part can only be allocated in the heap. The copy constructor
	// doesn't make any sense for this class.
	sparse_row_part(const sparse_row_part &part) = delete;
public:
	static size_t get_size(size_t num_non_zeros) {
		return sizeof(sparse_row_part) + sizeof(uint16_t) * num_non_zeros;
	}

	static size_t get_num_entries(size_t size) {
		assert((size - sizeof(sparse_row_part)) % sizeof(uint16_t) == 0);
		return (size - sizeof(sparse_row_part)) / sizeof(uint16_t);
	}

	static size_t get_row_id_size() {
		return sizeof(rel_row_idx);
	}

	static size_t get_col_entry_size() {
		return sizeof(rel_col_idxs[0]);
	}

	sparse_row_part(uint16_t rel_row_idx) {
		this->rel_row_idx = rel_row_idx;
		// The highest bit indicates that a new row part starts.
		this->rel_row_idx |= 1 << (num_bits - 1);
	}

	size_t get_rel_row_idx() const {
		static const int mask = (1 << (num_bits - 1)) - 1;
		return rel_row_idx & mask;
	}

	size_t get_rel_col_idx(size_t idx) const {
		return rel_col_idxs[idx];
	}

	rp_edge_iterator get_edge_iterator() {
		return rp_edge_iterator(get_rel_row_idx(), rel_col_idxs);
	}

	rp_edge_iterator get_edge_iterator(const char *data, size_t entry_size) {
		return rp_edge_iterator(get_rel_row_idx(), rel_col_idxs, data,
				entry_size);
	}
};

/*
 * This is used to store the rows with a single non-zero entry.
 */
class local_coo_t
{
	static const int num_bits = sizeof(uint16_t) * 8;
	// A block has at most 2^15 rows and cols.
	// The most significant bit in the row index is always 1, so this can
	// be interpreted as a row part as well.
	uint16_t row_idx;
	uint16_t col_idx;
public:
	local_coo_t(uint16_t row_idx, uint16_t col_idx) {
		this->row_idx = row_idx | (1 << (num_bits - 1));
		this->col_idx = col_idx;
	}

	uint16_t get_row_idx() const {
		static const int mask = (1 << (num_bits - 1)) - 1;
		return row_idx & mask;
	}

	uint16_t get_col_idx() const {
		return col_idx;
	}
};

/*
 * This stores the block of a sparse matrix partitioned in 2D.
 * The block is stored in row major.
 */
class sparse_block_2d
{
	uint32_t block_row_idx;
	uint32_t block_col_idx;
	// The total number of non-zero entries.
	uint32_t nnz;
	// The number of rows with non-zero entries.
	uint16_t nrow;
	// The number of rows with a single non-zero entry.
	uint16_t num_coo_vals;
	// This is where the row parts are serialized.
	char row_parts[0];

	// The 2D-partitioned block has to be allowed in the heap. The copy
	// constructor doesn't make sense for it.
	sparse_block_2d(const sparse_block_2d &block) = delete;

	/*
	 * Row header size including the COO region.
	 */
	size_t get_rindex_size() const {
		// The space used by row ids.
		return nrow * sparse_row_part::get_row_id_size()
			// The space used by col index entries in the row part.
			+ nnz * sparse_row_part::get_col_entry_size()
			// The empty row part that indicates the end of row-part region.
			+ sparse_row_part::get_row_id_size();
		// TODO should I align to the size of non-zero entry type?
	}

	/*
	 * The row header size excluding the COO region.
	 */
	size_t get_rheader_size() const {
		// The space used by row ids. The empty row part is also included.
		return (nrow - num_coo_vals) * sparse_row_part::get_row_id_size()
			// The space used by col index entries in the row part.
			+ (nnz - num_coo_vals) * sparse_row_part::get_col_entry_size()
			// The empty row part that indicates the end of row-part region.
			+ sparse_row_part::get_row_id_size();
	}

	sparse_row_part *get_rpart_end() {
		return (sparse_row_part *) (row_parts + get_rheader_size()
				// Let's exclude the last empty row part.
				- sparse_row_part::get_row_id_size());
	}
	local_coo_t *get_coo_start() {
		return (local_coo_t *) (row_parts + get_rheader_size());
	}

	char *get_nz_data() {
		return (char *) (row_parts + get_rindex_size());
	}

	const char *get_nz_data() const {
		return (char *) (row_parts + get_rindex_size());
	}
public:
	sparse_block_2d(uint32_t block_row_idx, uint32_t block_col_idx) {
		this->block_row_idx = block_row_idx;
		this->block_col_idx = block_col_idx;
		nnz = 0;
		nrow = 0;
		num_coo_vals = 0;
	}

	size_t get_block_row_idx() const {
		return block_row_idx;
	}

	size_t get_block_col_idx() const {
		return block_col_idx;
	}

	bool is_empty() const {
		return nnz == 0;
	}

	size_t get_size(size_t entry_size) const {
		if (entry_size == 0)
			return sizeof(*this) + get_rindex_size();
		else
			return sizeof(*this) + get_rindex_size() + entry_size * nnz;
	}

	/*
	 * This doesn't count the COO entries even though they also follow
	 * the row part format.
	 */
	bool has_rparts() const {
		return nnz - num_coo_vals > 0;
	}

	rp_edge_iterator get_first_edge_iterator() const {
		assert(has_rparts());
		// Discard the const qualifier
		// TODO I should make a const edge iterator
		sparse_row_part *rp = (sparse_row_part *) row_parts;
		return rp->get_edge_iterator();
	}

	rp_edge_iterator get_first_edge_iterator(size_t entry_size) const {
		assert(has_rparts());
		// Discard the const qualifier
		// TODO I should make a const edge iterator
		sparse_row_part *rp = (sparse_row_part *) row_parts;
		return rp->get_edge_iterator(get_nz_data(), entry_size);
	}

	rp_edge_iterator get_next_edge_iterator(const rp_edge_iterator &it) const {
		assert(!it.has_next());
		// TODO I should make a const edge iterator
		sparse_row_part *rp = (sparse_row_part *) it.get_curr_addr();
		return rp->get_edge_iterator();
	}

	rp_edge_iterator get_next_edge_iterator(const rp_edge_iterator &it,
			size_t entry_size) const {
		assert(!it.has_next());
		// TODO I should make a const edge iterator
		sparse_row_part *rp = (sparse_row_part *) it.get_curr_addr();
		return rp->get_edge_iterator(it.get_curr_data(), entry_size);
	}

	size_t get_num_coo_vals() const {
		return num_coo_vals;
	}
	const local_coo_t *get_coo_start() const {
		return (local_coo_t *) (row_parts + get_rheader_size());
	}

	void append(const sparse_row_part &part, size_t part_size);

	void add_coo(const std::vector<coo_nz_t> &nnz,
			const block_2d_size &block_size);

	/*
	 * Finalize the construction of the block.
	 * We need to add an empty row in the end if there aren't COO entries
	 * in the block, so the edge iterator can notice the end of the last row.
	 * We also need to add the non-zero values.
	 */
	void finalize(const char *data, size_t num_bytes);

	/*
	 * Get all non-zero entries in the block.
	 * This is used for testing only.
	 */
	std::vector<coo_nz_t> get_non_zeros(const block_2d_size &block_size) const;

	size_t get_nnz() const {
		return nnz;
	}

	void verify(const block_2d_size &block_size) const;
};

/*
 * This allows users to iterate the blocks in a block row.
 */
class block_row_iterator
{
	const sparse_block_2d *block;
	const sparse_block_2d *end;
public:
	/*
	 * `first' is the first block in the block row.
	 * `end' is the end of the last block in the block row.
	 */
	block_row_iterator(const sparse_block_2d *first, const sparse_block_2d *end) {
		block = first;
		this->end = end;
	}

	bool has_next() const {
		return block < end;
	}

	const sparse_block_2d &get_curr() const {
		return *block;
	}

	const sparse_block_2d &next(size_t entry_size) {
		const sparse_block_2d *orig = block;
		block = (const sparse_block_2d *) (((const char *) block)
			+ block->get_size(entry_size));
		return *orig;
	}
};

/*
 * This indexes a sparse matrix for easy access to the 2D-partitioned blocks.
 */
class SpM_2d_index
{
	struct deleter {
		void operator()(SpM_2d_index *p) const{
			free(p);
		}
	};

	matrix_header header;
	off_t offs[0];

	SpM_2d_index(const matrix_header &_header): header(_header) {
	}

	/*
	 * The number of offset entries in the index.
	 */
	size_t get_num_entries() const {
		return get_num_block_rows() + 1;
	}
public:
	typedef std::shared_ptr<SpM_2d_index> ptr;

	/*
	 * Calculate the storage size of the index.
	 */
	static size_t get_size(size_t num_entries) {
		return sizeof(SpM_2d_index)
			+ num_entries * sizeof(SpM_2d_index::offs[0]);
	}

	static ptr create(const matrix_header &header,
			const std::vector<off_t> &offs);
	static ptr load(const std::string &idx_file);
	static ptr safs_load(const std::string &idx_file);

	void verify() const;
	size_t get_num_block_rows() const;
	void dump(const std::string &file) const;
	void safs_dump(const std::string &file) const;
	off_t get_block_row_off(size_t idx) const;
	const matrix_header &get_header() const {
		return header;
	}
};

class vector_vector;

class SpM_2d_storage
{
	struct deleter {
		void operator()(char *p) const {
			free(p);
		}
	};

	std::string mat_name;
	int mat_file_id;

	std::shared_ptr<char> data;
	SpM_2d_index::ptr index;

	SpM_2d_storage(std::shared_ptr<char> data,
			SpM_2d_index::ptr index, const std::string mat_name) {
		this->data = data;
		this->index = index;
		this->mat_name = mat_name;
		this->mat_file_id = -1;
	}
public:
	typedef std::shared_ptr<SpM_2d_storage> ptr;

	static ptr safs_load(const std::string &mat_file,
			SpM_2d_index::ptr index);
	static ptr load(const std::string &mat_file,
			SpM_2d_index::ptr index);
	static ptr create(const matrix_header &header, const vector_vector &vv,
			SpM_2d_index::ptr index);

	static void verify(SpM_2d_index::ptr index, const std::string &mat_file);

	block_row_iterator get_block_row_it(size_t idx) const {
		char *start = data.get() + index->get_block_row_off(idx);
		char *end = data.get() + index->get_block_row_off(idx + 1);
		return block_row_iterator((const sparse_block_2d *) start,
				(const sparse_block_2d *) end);
	}

	size_t get_num_block_rows() const {
		return index->get_num_block_rows();
	}

	void verify() const;

	std::shared_ptr<safs::file_io_factory> create_io_factory() const;
};

}

#endif
