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

#include "fm_utils.h"

#include "factor.h"
#include "generic_type.h"
#include "local_vec_store.h"
#include "data_frame.h"
#include "vector_vector.h"
#include "local_vv_store.h"
#include "mem_vv_store.h"
#include "EM_vv_store.h"
#include "EM_vector.h"
#include "local_vec_store.h"
#include "sparse_matrix.h"

namespace fm
{

static bool deduplicate = false;
// remove self edges. It's enabled by default.
static bool remove_selfe = true;

void set_deduplicate(bool v)
{
	deduplicate = v;
}

void set_remove_self_edge(bool v)
{
	remove_selfe = v;
}

/*
 * This is a simplifed ext_mem_undirected_vertex.
 */
class serialized_row
{
	ele_idx_t id;
	uint32_t attr_size;
	ele_idx_t num_edges;
	ele_idx_t nz_idxs[0];
public:
	static size_t cal_row_size(size_t nnz, size_t attr_size) {
		return sizeof(serialized_row) + sizeof(nz_idxs[0]) * nnz + attr_size * nnz;
	}

	serialized_row(ele_idx_t id, uint32_t data_size, size_t num_edges) {
		this->id = id;
		this->attr_size = data_size;
		this->num_edges = num_edges;
	}

	size_t get_size() const {
		return sizeof(serialized_row) + sizeof(nz_idxs[0]) * num_edges
			+ attr_size * num_edges;
	}

	size_t get_row_id() const {
		return id;
	}

	size_t get_nnz() const {
		return num_edges;
	}

	size_t get_nz_idx(size_t i) const {
		return nz_idxs[i];
	}

	char *get_data(size_t i) {
		return (char *) (nz_idxs + num_edges) + i * attr_size;
	}

	const char *get_data(size_t i) const {
		return (char *) (nz_idxs + num_edges) + i * attr_size;
	}

	void set_val(size_t i, size_t nz_idx) {
		assert(i < num_edges);
		nz_idxs[i] = nz_idx;
	}

	void set_val(size_t i, size_t nz_idx, const char *val) {
		assert(i < num_edges);
		nz_idxs[i] = nz_idx;
		memcpy(get_data(i), val, attr_size);
	}
};

/*
 * This applies to a data frame and generate a sparse matrix with 2D
 * partitioning.
 */
class SpM_apply_operate: public gr_apply_operate<sub_data_frame>
{
	matrix_header mheader;
	block_2d_size block_size;
	size_t num_cols;
	size_t attr_size;
public:
	SpM_apply_operate(const matrix_header &header) {
		this->mheader = header;
		this->block_size = header.get_2d_block_size();
		this->num_cols = header.get_num_cols();
		this->attr_size = header.get_entry_size();
	}

	size_t get_attr_size() const {
		return attr_size;
	}

	virtual bool ignore_key(const void *key) const {
		ele_idx_t vid = *(const ele_idx_t *) key;
		return vid == INVALID_IDX_VAL;
	}

	void run(const void *key, const sub_data_frame &val,
			local_vec_store &out) const;
	virtual void run_row(const void *key, const sub_data_frame &val,
			local_vec_store &out) const = 0;
	void run_blocks(const void *key, const local_vv_store &val,
			local_vec_store &out) const;

	const scalar_type &get_key_type() const {
		return get_scalar_type<ele_idx_t>();
	}

	const scalar_type &get_output_type() const {
		return get_scalar_type<char>();
	}

	size_t get_num_out_eles() const {
		return 1;
	}
};

/*
 * This works for a binary sparse matrix.
 */
class binary_apply_operate: public SpM_apply_operate
{
public:
	binary_apply_operate(const matrix_header &header): SpM_apply_operate(header) {
	}
	virtual void run_row(const void *key, const sub_data_frame &val,
			local_vec_store &out) const;
};

/*
 * This works for a sparse matrix with arbitrary non-zero values.
 */
template<class T>
class attr_apply_operate: public SpM_apply_operate
{
	typedef std::pair<ele_idx_t, T> edge_type;

	struct edge_compare {
		bool operator()(const edge_type &e1, const edge_type &e2) {
			return e1.first < e2.first;
		}
	};

	struct edge_predicate {
		bool operator()(const edge_type &e1, const edge_type &e2) {
			return e1.first == e2.first;
		}
	};
public:
	attr_apply_operate(const matrix_header &header): SpM_apply_operate(header) {
		assert(sizeof(T) == header.get_entry_size());
	}
	virtual void run_row(const void *key, const sub_data_frame &val,
			local_vec_store &out) const;
};

struct block_info
{
	// The number of non-zero entries.
	uint32_t nnz;
	// The number of rows with non-zero entries.
	uint16_t nrow;
	// The number of rows with a single non-zero entry.
	uint16_t num_coos;

	block_info() {
		nnz = 0;
		nrow = 0;
		num_coos = 0;
	}
};

struct block_pointers
{
	sparse_block_2d *block;
	// The pointer to the coo region in the block.
	local_coo_t *coos;
	// The pointer to the values of SCSR in the block.
	char *nz_vals;
	// The pointer to the values of coos in the block.
	char *coo_vals;

	block_pointers() {
		block = NULL;
		coos = NULL;
		nz_vals = NULL;
		coo_vals = NULL;
	}
};

static void add_nz(block_pointers &ps, const serialized_row &v,
		const std::vector<ele_idx_t> &neighs,
		block_2d_size block_size, size_t nz_size,
		char *buf)
{
	// There is only one edge from the vertex in the block
	// We store it in the COO region.
	if (neighs.size() == 1) {
		// Add index
		*ps.coos = local_coo_t(v.get_row_id() & block_size.get_nrow_mask(),
				neighs[0] & block_size.get_ncol_mask());
		// Add the non-zero value.
		if (nz_size > 0)
			memcpy(ps.coo_vals, v.get_data(0), nz_size);

		// move to the next one.
		ps.coos++;
		if (nz_size > 0)
			ps.coo_vals += nz_size;
	}
	// There are more than one edge from the vertex.
	// We store it in SCSR.
	else if (neighs.size() > 1) {
		// Add the index
		size_t row_idx = v.get_row_id() & block_size.get_nrow_mask();
		sparse_row_part *part = new (buf) sparse_row_part(row_idx);
		rp_edge_iterator edge_it = part->get_edge_iterator();
		for (size_t k = 0; k < neighs.size(); k++)
			edge_it.append(block_size, neighs[k]);
		ps.block->append(*part, sparse_row_part::get_size(neighs.size()));
		// Add the non-zero values.
		if (nz_size > 0)
			memcpy(ps.nz_vals, v.get_data(0), nz_size * neighs.size());

		// Move to the next one.
		ps.nz_vals += nz_size * neighs.size();
	}
}

size_t collect_block_info(size_t block_row_id, const local_vv_store &val,
		block_2d_size block_size, size_t nz_size,
		std::vector<block_info> &block_infos)
{
	// We calculate the number of bytes in each block accurately.
	for (size_t i = 0; i < val.get_num_vecs(); i++) {
		const serialized_row *v = (const serialized_row *) val.get_raw_arr(i);
		assert(val.get_length(i) == v->get_size());
		assert((v->get_row_id() >> block_size.get_nrow_log()) == block_row_id);

		// The id of the current block.
		size_t curr_bid = 0;
		// The number of edges in the block.
		size_t num_edges_block = 0;
		for (size_t j = 0; j < v->get_nnz(); j++) {
			ele_idx_t id = v->get_nz_idx(j);
			size_t block_id = id >> block_size.get_ncol_log();
			if (curr_bid == block_id)
				num_edges_block++;
			// this can happen for the first block.
			else if (num_edges_block == 0) {
				curr_bid = block_id;
				num_edges_block = 1;
			}
			else {
				block_infos[curr_bid].nrow++;
				block_infos[curr_bid].nnz += num_edges_block;
				// There is only one edge from the vertex in the block
				// We store it in the COO region.
				if (num_edges_block == 1)
					block_infos[curr_bid].num_coos++;
				curr_bid = block_id;
				num_edges_block = 1;
			}
		}

		// For the last block
		if (num_edges_block > 0) {
			block_infos[curr_bid].nrow++;
			block_infos[curr_bid].nnz += num_edges_block;
			// There is only one edge from the vertex in the block
			// We store it in the COO region.
			if (num_edges_block == 1)
				block_infos[curr_bid].num_coos++;
		}
	}

	// Get the total number of bytes in the output vector.
	size_t tot_num_bytes = 0;
	for (size_t i = 0; i < block_infos.size(); i++) {
		// We don't store empty blocks.
		if (block_infos[i].nnz == 0)
			continue;

		sparse_block_2d block(0, 0, block_infos[i].nnz, block_infos[i].nrow,
				block_infos[i].num_coos);
		tot_num_bytes += block.get_size(nz_size);
	}

	return tot_num_bytes;
}

void SpM_apply_operate::run_blocks(const void *key, const local_vv_store &val,
		local_vec_store &out) const
{
	factor_value_t block_row_id = *(const factor_value_t *) key;

	size_t num_blocks = div_ceil(num_cols, block_size.get_num_cols());
	std::vector<block_info> block_infos(num_blocks);
	size_t tot_num_bytes = collect_block_info(block_row_id, val, block_size,
			attr_size, block_infos);

	// If the block row doesn't have any non-zero entries, let's insert an
	// empty block row, so the matrix index can work correctly.
	if (tot_num_bytes == 0) {
		if (block_row_id == 0) {
			out.resize(sizeof(sparse_block_2d) + sizeof(matrix_header));
			memcpy(out.get_raw_arr(), &mheader, sizeof(mheader));
		}
		else
			out.resize(sizeof(sparse_block_2d));
		sparse_block_2d *block = new (out.get_raw_arr()) sparse_block_2d(
					block_row_id, 0);
		assert(block->is_empty());
		return;
	}

	char *arr_start;
	if (block_row_id == 0) {
		out.resize(tot_num_bytes + sizeof(matrix_header));
		arr_start = out.get_raw_arr() + sizeof(matrix_header);
		memcpy(out.get_raw_arr(), &mheader, sizeof(mheader));
	}
	else {
		out.resize(tot_num_bytes);
		arr_start = out.get_raw_arr();
	}

	// Get the location where we can fill data to.
	std::vector<block_pointers> blocks(num_blocks);
	size_t curr_size = 0;
	for (size_t i = 0; i < block_infos.size(); i++) {
		// If the block doesn't have non-zero entries, we will skip it.
		if (block_infos[i].nnz == 0)
			continue;

		sparse_block_2d *block
			= new (arr_start + curr_size) sparse_block_2d(
					block_row_id, i, block_infos[i].nnz, block_infos[i].nrow,
					block_infos[i].num_coos);
		curr_size += block->get_size(attr_size);
		blocks[i].coos = block->get_coo_start();
		blocks[i].coo_vals = block->get_coo_val_start(attr_size);
		blocks[i].nz_vals = block->get_nz_data();
		// For now we only include single-entry rows.
		// We'll add other non-zero entries later.
		blocks[i].block = new (block) sparse_block_2d(
				block_row_id, i, block_infos[i].num_coos,
				block_infos[i].num_coos, block_infos[i].num_coos);
	}

	size_t max_row_size = sparse_row_part::get_size(block_size.get_num_cols());
	std::unique_ptr<char[]> buf = std::unique_ptr<char[]>(new char[max_row_size]);
	// Serialize data.
	std::vector<ele_idx_t> neighs;
	for (size_t i = 0; i < val.get_num_vecs(); i++) {
		const serialized_row *v = (const serialized_row *) val.get_raw_arr(i);

		size_t curr_bid = 0;
		neighs.clear();
		for (size_t j = 0; j < v->get_nnz(); j++) {
			ele_idx_t id = v->get_nz_idx(j);
			size_t block_id = id >> block_size.get_ncol_log();
			if (curr_bid == block_id)
				neighs.push_back(id);
			else {
				if (neighs.size() > block_size.get_num_cols()) {
					BOOST_LOG_TRIVIAL(error)
						<< "ERROR! There are more neighbors than the block size";
					return;
				}
				add_nz(blocks[curr_bid], *v, neighs, block_size, attr_size,
						buf.get());
				curr_bid = block_id;
				neighs.clear();
				neighs.push_back(id);
			}
		}

		if (neighs.size() > block_size.get_num_cols()) {
			BOOST_LOG_TRIVIAL(error)
				<< "ERROR! There are more neighbors than the block size";
			return;
		}
		add_nz(blocks[curr_bid], *v, neighs, block_size, attr_size, buf.get());
	}

	// Finalize each block.
	for (size_t i = 0; i < blocks.size(); i++) {
		if (blocks[i].block) {
			// We have fill non-zero values to the output vector,
			// so we can pass NULL here.
			blocks[i].block->finalize(NULL, 0);
//			blocks[i].block->verify(block_size);
		}
	}
}

void binary_apply_operate::run_row(const void *key, const sub_data_frame &val,
		local_vec_store &out) const
{
	ele_idx_t vid = *(const ele_idx_t *) key;
	if (vid == INVALID_IDX_VAL) {
		out.resize(0);
		return;
	}

	assert(val.size() == 2);
	assert(out.is_type<char>());
	// The data frame is sorted based on the first vector and now we need
	// to access the entries in the second vector.
	const local_vec_store &vec = *val[1];
	assert(vec.get_type() == get_scalar_type<ele_idx_t>());

	// I added an invalid edge for each vertex.
	// The invalid edge is the maximal integer.
	size_t num_edges = vec.get_length() - 1;
	// TODO we actually don't need to alloate memory multiple times.
	std::unique_ptr<ele_idx_t[]> edge_buf
		= std::unique_ptr<ele_idx_t[]>(new ele_idx_t[num_edges]);
	size_t edge_idx = 0;
	for (size_t i = 0; i < vec.get_length(); i++) {
		if (vec.get<ele_idx_t>(i) == INVALID_IDX_VAL
				// skip self-edges.
				|| (remove_selfe && vec.get<ele_idx_t>(i) == vid))
			continue;
		edge_buf[edge_idx++] = vec.get<ele_idx_t>(i);
	}
	assert(edge_idx <= num_edges);
	// TODO we need to shuffle non-zero values accordingly.
	// If there are self-edges, edge_idx has the actual number of edges.
	num_edges = edge_idx;
	std::sort(edge_buf.get(), edge_buf.get() + num_edges);
	if (deduplicate) {
		ele_idx_t *end = std::unique(edge_buf.get(),
				edge_buf.get() + num_edges);
		num_edges = end - edge_buf.get();
	}
	size_t size = serialized_row::cal_row_size(num_edges, get_attr_size());
	out.resize(size);

	new (out.get_raw_arr()) serialized_row(vid, get_attr_size(), num_edges);
	serialized_row *row = (serialized_row *) out.get_raw_arr();
	for (size_t i = 0; i < num_edges; i++)
		row->set_val(i, edge_buf[i]);
}

template<class T>
void attr_apply_operate<T>::run_row(const void *key, const sub_data_frame &val,
		local_vec_store &out) const
{
	ele_idx_t vid = *(const ele_idx_t *) key;
	if (vid == INVALID_IDX_VAL) {
		out.resize(0);
		return;
	}

	assert(val.size() == 3);
	assert(out.is_type<char>());
	// The data frame is sorted based on the first vector and now we need
	// to access the entries in the second vector.
	const local_vec_store &vec = *val[1];
	const local_vec_store &attr_vec = *val[2];
	assert(vec.get_type() == get_scalar_type<ele_idx_t>());

	// I added an invalid edge for each vertex.
	// The invalid edge is the maximal integer.
	ele_idx_t num_edges = vec.get_length() - 1;
	// TODO we actually don't need to alloate memory multiple times.
	std::unique_ptr<edge_type[]> edge_buf
		= std::unique_ptr<edge_type[]>(new edge_type[num_edges]);
	size_t edge_idx = 0;
	for (size_t i = 0; i < vec.get_length(); i++) {
		if (vec.get<ele_idx_t>(i) == INVALID_IDX_VAL
				// skip self-edges.
				|| (remove_selfe && vec.get<ele_idx_t>(i) == vid))
			continue;
		edge_buf[edge_idx].first = vec.get<ele_idx_t>(i);
		edge_buf[edge_idx].second = attr_vec.get<T>(i);
		edge_idx++;
	}
	assert(edge_idx <= num_edges);
	// If there are self-edges, edge_idx has the actual number of edges.
	num_edges = edge_idx;
	std::sort(edge_buf.get(), edge_buf.get() + num_edges, edge_compare());
	if (deduplicate) {
		edge_type *end = std::unique(edge_buf.get(),
				edge_buf.get() + num_edges, edge_predicate());
		num_edges = end - edge_buf.get();
	}
	assert(get_attr_size() == attr_vec.get_entry_size());
	size_t size = serialized_row::cal_row_size(num_edges, get_attr_size());
	out.resize(size);

	new (out.get_raw_arr()) serialized_row(vid, get_attr_size(), num_edges);
	serialized_row *row = (serialized_row *) out.get_raw_arr();
	for (size_t i = 0; i < num_edges; i++)
		row->set_val(i, edge_buf[i].first, (const char *) &edge_buf[i].second);
}

static void expose_portion(const sub_data_frame &sub_df, off_t loc, size_t length)
{
	for (size_t i = 0; i < sub_df.size(); i++)
		const_cast<local_vec_store *>(sub_df[i].get())->expose_sub_vec(loc, length);
}

/*
 * The data frame contains the non-zero entries in a sparse matrix.
 * We group non-zero entries from multiple contiguous rows.
 * The first column of `df' stores the row Ids of non-zero entries.
 *
 * When we convert an edge list into a sparse matrix with 2D partitioning,
 * it has two stages: we first construct rows in the local memory buffer
 * (the rows have the same format as serialized FlashGraph vertex format)
 * and then convert them into blocks.
 */
void SpM_apply_operate::run(const void *key, const sub_data_frame &df,
		local_vec_store &out) const
{
	detail::mem_vv_store::ptr vv_buf = detail::mem_vv_store::create(
			get_scalar_type<char>());
	local_buf_vec_store::ptr row_buf(new local_buf_vec_store(0,
				1, get_scalar_type<char>(), -1));
	sub_data_frame tmp_df(df.begin() + 1, df.end());
	auto sorted_col = tmp_df[0];
	agg_operate::const_ptr find_next
		= sorted_col->get_type().get_agg_ops().get_find_next();
	size_t loc = 0;
	size_t col_len = sorted_col->get_length();
	// We can't search a vv store.
	assert(!detail::vv_store::is_vector_vector(*sorted_col));
	const char *start = sorted_col->get_raw_arr();
	size_t entry_size = sorted_col->get_type().get_size();
	size_t orig_local_start = tmp_df[0]->get_local_start();
	size_t orig_len = tmp_df[0]->get_length();
	while (loc < col_len) {
		size_t curr_length = col_len - loc;
		const char *curr_ptr = start + entry_size * loc;
		size_t rel_end;
		find_next->runAgg(curr_length, curr_ptr, &rel_end);
		// This expose a portion of the data frame.
		expose_portion(tmp_df, orig_local_start + loc, rel_end);
		// The first argument is the key and the second one is the value
		// (a data frame)
		run_row(curr_ptr, tmp_df, *row_buf);
		if (row_buf->get_length() > 0)
			vv_buf->append(*row_buf);
		loc += rel_end;
	}
	expose_portion(tmp_df, orig_local_start, orig_len);
	if (vv_buf->get_num_vecs() > 0) {
		local_vv_store::const_ptr lvv
			= std::dynamic_pointer_cast<const local_vv_store>(
					vv_buf->get_portion(0, vv_buf->get_num_vecs()));
		run_blocks(key, *lvv, out);
	}
	else
		out.resize(0);
}

namespace
{

struct unit4
{
	char data[4];
};

struct unit8
{
	char data[8];
};

struct empty_deleter {
	void operator()(detail::vec_store *addr) {
	}
};

}

std::pair<fm::SpM_2d_index::ptr, vector_vector::ptr> create_2d_matrix(
		data_frame::const_ptr df, const block_2d_size &block_size, size_t num_rows,
		size_t num_cols, const fm::scalar_type *entry_type)
{
	if (df->get_num_vecs() == 1) {
		BOOST_LOG_TRIVIAL(error)
			<< "The data frame needs to have at least two columns";
		return std::pair<fm::SpM_2d_index::ptr, vector_vector::ptr>();
	}
	struct timeval start, end;
	gettimeofday(&start, NULL);
	auto sorted_df = df->sort(df->get_vec_name(0));
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to sort the edge list\n",
			time_diff(start, end));

	size_t attr_size = 0;
	detail::vec_store::const_ptr row_idxs = sorted_df->get_vec(0);
	detail::vec_store::const_ptr col_idxs = sorted_df->get_vec(1);
	// Construct groupby keys
	// TODO in the future, we should evaluate this lazily.
	vector::ptr key = vector::create(row_idxs);
	dense_matrix::ptr key_mat = key->conv2mat(key->get_length(), 1, false);
	scalar_variable::ptr var(new scalar_variable_impl<ele_idx_t>(
				block_size.get_num_rows()));
	key_mat = key_mat->apply_scalar(var,
			bulk_operate::conv2ptr(get_scalar_type<ele_idx_t>().get_basic_ops().get_divide()));
	key_mat = key_mat->cast_ele_type(get_scalar_type<ele_idx_t>());
	key = key_mat->conv2vec();

	data_frame::ptr groupby_df = data_frame::create();
	// TODO this is a very ugly solution.
	groupby_df->add_vec("key",
			std::shared_ptr<detail::vec_store>(const_cast<detail::vec_store *>(key->get_raw_store().get()),
				empty_deleter()));
	groupby_df->add_vec("x",
			std::shared_ptr<detail::vec_store>(const_cast<detail::vec_store *>(row_idxs.get()),
				empty_deleter()));
	groupby_df->add_vec("y",
			std::shared_ptr<detail::vec_store>(const_cast<detail::vec_store *>(col_idxs.get()),
				empty_deleter()));
	if (sorted_df->get_num_vecs() > 2) {
		attr_size = sorted_df->get_vec(2)->get_type().get_size();
		groupby_df->add_vec("attr",
				std::shared_ptr<detail::vec_store>(const_cast<detail::vec_store *>(sorted_df->get_vec(2).get()),
					empty_deleter()));
	}

	prim_type type = prim_type::P_BOOL;
	if (attr_size)
		type = groupby_df->get_vec("attr")->get_type().get_type();
	matrix_header mheader(matrix_type::SPARSE, attr_size, num_rows,
			num_cols, matrix_layout_t::L_ROW_2D, type, block_size);
	mheader.verify();

	gettimeofday(&start, NULL);
	vector_vector::ptr res;
	if (attr_size == 0) {
		std::unique_ptr<SpM_apply_operate> op(
				new binary_apply_operate(mheader));
		res = groupby_df->groupby("key", *op);
	}
	// Instead of giving the real data type, we give a type that indicates
	// the size of the edge data size. Actually, we don't interpret data type
	// here. Only the data size matters.
	else if (attr_size == 4) {
		std::unique_ptr<attr_apply_operate<unit4> > op(
				new attr_apply_operate<unit4>(mheader));
		res = groupby_df->groupby("key", *op);
	}
	else if (attr_size == 8) {
		std::unique_ptr<attr_apply_operate<unit8> > op(
				new attr_apply_operate<unit8>(mheader));
		res = groupby_df->groupby("key", *op);
	}
	else {
		BOOST_LOG_TRIVIAL(error)
			<< "The edge attribute has an unsupported type";
		return std::pair<fm::SpM_2d_index::ptr, vector_vector::ptr>();
	}
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to groupby the edge list.\n",
			time_diff(start, end));

	// Construct the index file of the adjacency matrix.
	std::vector<off_t> offsets(res->get_num_vecs() + 1);
	off_t off = 0;
	for (size_t i = 0; i < res->get_num_vecs(); i++) {
		offsets[i] = off;
		off += res->get_length(i);
	}
	// The first block row actually contains the matrix header.
	// We need to skip the matrix header.
	offsets[0] = sizeof(mheader);
	offsets[res->get_num_vecs()] = off;
	SpM_2d_index::ptr idx = SpM_2d_index::create(mheader, offsets);

	return std::pair<SpM_2d_index::ptr, vector_vector::ptr>(idx, res);
}

static inline bool set_persistent(detail::vec_store::const_ptr vec,
		const std::string name)
{
	detail::EM_vec_store::const_ptr em_vec
		= std::dynamic_pointer_cast<const detail::EM_vec_store>(vec);
	detail::EM_vec_store &em_vec1 = const_cast<detail::EM_vec_store &>(*em_vec);
	return em_vec1.set_persistent(name);
}

std::shared_ptr<sparse_matrix> create_2d_matrix(data_frame::ptr df,
		const block_2d_size &block_size, const fm::scalar_type *entry_type,
		bool is_sym, const std::string &name)
{
	ele_idx_t max_vid;
	{
		vector::ptr vec = vector::create(df->get_vec(0));
		assert(vec->get_type() == get_scalar_type<ele_idx_t>());
		max_vid = vec->max<ele_idx_t>();
		vec = vector::create(df->get_vec(1));
		assert(vec->get_type() == get_scalar_type<ele_idx_t>());
		max_vid = std::max(max_vid, vec->max<ele_idx_t>());
	}
	printf("max id: %d\n", max_vid);

	detail::vec_store::ptr seq_vec = detail::create_seq_vec_store<ele_idx_t>(
			0, max_vid, 1);
	detail::vec_store::ptr rep_vec = detail::create_rep_vec_store<ele_idx_t>(
			max_vid + 1, INVALID_IDX_VAL);
	detail::vec_store::ptr attr_extra;
	detail::vec_store::ptr attr_vec;
	if (df->get_num_vecs() > 2)
		attr_vec = df->get_vec(2);
	if (attr_vec && attr_vec->get_type() == get_scalar_type<int>())
		attr_extra = detail::create_rep_vec_store<int>(max_vid + 1, 0);
	else if (attr_vec && attr_vec->get_type() == get_scalar_type<long>())
		attr_extra = detail::create_rep_vec_store<long>(max_vid + 1, 0);
	else if (attr_vec && attr_vec->get_type() == get_scalar_type<float>())
		attr_extra = detail::create_rep_vec_store<float>(max_vid + 1, 0);
	else if (attr_vec && attr_vec->get_type() == get_scalar_type<double>())
		attr_extra = detail::create_rep_vec_store<double>(max_vid + 1, 0);
	else if (attr_vec) {
		BOOST_LOG_TRIVIAL(error) << "unknown attribute type";
		return sparse_matrix::ptr();
	}
	assert(seq_vec->get_length() == rep_vec->get_length());

	// I artificially add an invalid out-edge for each vertex, so it's
	// guaranteed that each vertex exists in the adjacency lists.
	data_frame::ptr new_df = data_frame::create();
	new_df->add_vec(df->get_vec_name(0), seq_vec);
	new_df->add_vec(df->get_vec_name(1), rep_vec);
	if (df->get_num_vecs() > 2) {
		assert(attr_extra);
		new_df->add_vec(df->get_vec_name(2), attr_extra);
	}
	df->append(new_df);

	// I artificially add an invalid in-edge for each vertex.
	new_df = data_frame::create();
	new_df->add_vec(df->get_vec_name(1), seq_vec);
	new_df->add_vec(df->get_vec_name(0), rep_vec);
	if (df->get_num_vecs() > 2) {
		assert(attr_extra);
		new_df->add_vec(df->get_vec_name(2), attr_extra);
	}
	df->append(new_df);

	std::string spm_name = name;
	if (spm_name == "")
		spm_name = std::string("spm") + gen_rand_name(10) + ".mat";

	auto out_mat = create_2d_matrix(df, block_size, max_vid + 1, max_vid + 1,
			entry_type);
	if (!is_sym) {
		data_frame::ptr reversed_df = data_frame::create();
		reversed_df->add_vec(df->get_vec_name(1), df->get_vec(1));
		reversed_df->add_vec(df->get_vec_name(0), df->get_vec(0));
		for (size_t i = 2; i < df->get_num_vecs(); i++)
			reversed_df->add_vec(df->get_vec_name(i), df->get_vec(i));
		auto in_mat = create_2d_matrix(reversed_df, block_size, max_vid + 1, max_vid + 1,
				entry_type);

		if (out_mat.second->get_raw_store()->is_in_mem()) {
			assert(in_mat.second->get_raw_store()->is_in_mem());
			SpM_2d_storage::ptr out_store = SpM_2d_storage::create(
					*out_mat.second, out_mat.first);
			SpM_2d_storage::ptr in_store = SpM_2d_storage::create(
					*in_mat.second, in_mat.first);
			return sparse_matrix::create(out_mat.first, out_store,
					in_mat.first, in_store);
		}
		else {
			vector::ptr out_vec = out_mat.second->cat();
			vector::ptr in_vec = in_mat.second->cat();
			std::string out_spm_name = std::string("out_") + spm_name;
			std::string in_spm_name = std::string("in_") + spm_name;
			bool ret1 = set_persistent(out_vec->get_raw_store(), out_spm_name);
			bool ret2 = set_persistent(in_vec->get_raw_store(), in_spm_name);
			if (!ret1 || !ret2)
				return sparse_matrix::ptr();
			safs::file_io_factory::shared_ptr out_factory = safs::create_io_factory(
					out_spm_name, safs::REMOTE_ACCESS);
			safs::file_io_factory::shared_ptr in_factory = safs::create_io_factory(
					in_spm_name, safs::REMOTE_ACCESS);
			return sparse_matrix::create(out_mat.first, out_factory,
					in_mat.first, in_factory);
		}
	}
	else {
		if (out_mat.second->get_raw_store()->is_in_mem()) {
			SpM_2d_storage::ptr store = SpM_2d_storage::create(
					*out_mat.second, out_mat.first);
			return sparse_matrix::create(out_mat.first, store);
		}
		else {
			vector::ptr vec = out_mat.second->cat();
			bool ret = set_persistent(vec->get_raw_store(), spm_name);
			if (!ret)
				return sparse_matrix::ptr();
			safs::file_io_factory::shared_ptr factory = safs::create_io_factory(
					spm_name, safs::REMOTE_ACCESS);
			return sparse_matrix::create(out_mat.first, factory);
		}
	}
}

}
