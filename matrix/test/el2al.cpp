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

#include <string.h>
#include <assert.h>

#include <vector>
#include <string>

#include "log.h"
#include "FG_basic_types.h"
#include "vertex.h"
#include "graph_file_header.h"
#include "vertex_index.h"

#include "generic_type.h"
#include "data_io.h"
#include "data_frame.h"
#include "mem_data_frame.h"
#include "mem_vector.h"
#include "vector_vector.h"
#include "factor.h"
#include "sparse_matrix_format.h"

using namespace fm;

/*
 * This class parses a line into an edge (source, destination).
 */
class edge_parser: public line_parser
{
public:
	size_t parse(const std::vector<std::string> &lines, data_frame &df) const;
	size_t get_num_cols() const {
		return 2;
	}

	std::string get_col_name(off_t idx) const {
		if (idx == 0)
			return "source";
		else
			return "dest";
	}

	const scalar_type &get_col_type(off_t idx) const {
		return get_scalar_type<fg::vertex_id_t>();
	}
};

size_t edge_parser::parse(const std::vector<std::string> &lines,
		data_frame &df) const
{
	type_mem_vector<fg::vertex_id_t>::ptr froms
		= type_mem_vector<fg::vertex_id_t>::create(lines.size());
	type_mem_vector<fg::vertex_id_t>::ptr tos
		= type_mem_vector<fg::vertex_id_t>::create(lines.size());
	for (size_t i = 0; i < lines.size(); i++) {
		const char *line = lines[i].c_str();
		int len = strlen(line);
		const char *first = line;
		// Skip space
		for (; isspace(*first); first++);
		// Make sure we get a number.
		if (!isdigit(*first)) {
			BOOST_LOG_TRIVIAL(error)
				<< std::string("the first entry isn't a number: ") + first;
			continue;
		}
		long from = atol(first);
		assert(from >= 0 && from < fg::MAX_VERTEX_ID);

		const char *second = first;
		// Go to the end of the first number.
		for (; isdigit(*second); second++);
		if (second - line == len) {
			BOOST_LOG_TRIVIAL(error)
				<< std::string("there isn't second entry: ") + line;
			continue;
		}
		// Skip space between two numbers.
		for (; isspace(*second); second++);
		// Make sure we get a number.
		if (!isdigit(*second)) {
			BOOST_LOG_TRIVIAL(error)
				<< std::string("the second entry isn't a number: ") + second;
			continue;
		}
		long to = atol(second);
		assert(to >= 0 && to < fg::MAX_VERTEX_ID);

		froms->set(i, from);
		tos->set(i, to);
	}

	df.get_vec(0)->append(*froms);
	df.get_vec(1)->append(*tos);
	return froms->get_length();
}

/*
 * This applies to a vector of values corresponding to the same key,
 * and generates an adjacency list.
 */
class adj_apply_operate: public gr_apply_operate<data_frame>
{
	fg::edge_type etype;
public:
	adj_apply_operate(fg::edge_type etype) {
		this->etype = etype;
	}

	void run(const void *key, const data_frame &val, mem_vector &out) const;

	const scalar_type &get_key_type() const {
		return get_scalar_type<fg::vertex_id_t>();
	}

	const scalar_type &get_output_type() const {
		return get_scalar_type<char>();
	}

	size_t get_num_out_eles() const {
		return 1;
	}
};

void adj_apply_operate::run(const void *key, const data_frame &val,
		mem_vector &out) const
{
	fg::vertex_id_t vid = *(const fg::vertex_id_t *) key;
	if (vid == fg::INVALID_VERTEX_ID) {
		out.resize(0);
		return;
	}

	// Right now, we just assume there aren't attributes.
	size_t edge_data_size = 0;
	assert(val.get_num_vecs() == 2);

	assert(out.is_type<char>());
	const type_mem_vector<fg::vertex_id_t> *vec;
	if (etype == fg::edge_type::OUT_EDGE)
		vec = (const type_mem_vector<fg::vertex_id_t> *) &val.get_vec_ref("dest");
	else
		vec = (const type_mem_vector<fg::vertex_id_t> *) &val.get_vec_ref("source");

	// I added an invalid edge for each vertex.
	// The invalid edge is the maximal integer.
	fg::vsize_t num_edges = vec->get_length() - 1;
	// TODO we actually don't need to alloate memory multiple times.
	std::unique_ptr<fg::vertex_id_t[]> edge_buf
		= std::unique_ptr<fg::vertex_id_t[]>(new fg::vertex_id_t[num_edges]);
	size_t edge_idx = 0;
	for (size_t i = 0; i < vec->get_length(); i++) {
		if (vec->get(i) == fg::INVALID_VERTEX_ID)
			continue;
		edge_buf[edge_idx++] = vec->get(i);
	}
	assert(edge_idx == num_edges);
	std::sort(edge_buf.get(), edge_buf.get() + num_edges);
	size_t size = fg::ext_mem_undirected_vertex::num_edges2vsize(num_edges,
			edge_data_size);
	out.resize(size);

	// Even if we generate a directed, we still can use undirected vertex to
	// store one type of edges of a vertex.
	fg::in_mem_undirected_vertex<> v(vid, edge_data_size > 0);
	for (size_t i = 0; i < num_edges; i++)
		v.add_edge(fg::edge<>(vid, edge_buf[i]));
	fg::ext_mem_undirected_vertex::serialize(v, out.get_raw_arr(), size, etype);
}

class count_agg_operate: public agg_operate
{
public:
	virtual void run(size_t num_eles, const void *in, void *output) const {
		fg::vsize_t *out = (fg::vsize_t *) output;
		out[0] = num_eles - 1;
	}

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<fg::vertex_id_t>();
	}

	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<fg::vsize_t>();
	}
};

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

void verify_2d_matrix(const std::string &mat_file, const std::string &mat_idx_file)
{
	SpM_2d_index::ptr idx = SpM_2d_index::load(mat_idx_file);
	SpM_2d_storage::ptr mat = SpM_2d_storage::load(mat_file, idx);
	mat->verify();
}

void create_2d_matrix(vector_vector::ptr out_adjs, vector_vector::ptr in_adjs,
		size_t num_vertices, const std::string &mat_name)
{
	std::string mat_file = mat_name + ".mat";
	std::string mat_idx_file = mat_name + ".mat_idx";
	std::string t_mat_file = mat_name + "_t.mat";
	std::string t_mat_idx_file = mat_name + "_t.mat_idx";

	// Construct 2D partitioning of the adjacency matrix.
	size_t block_height
		= ((size_t) std::numeric_limits<unsigned short>::max()) + 1;
	block_2d_size block_size(block_height, block_height);
	factor f(ceil(((double) num_vertices) / block_size.get_num_rows()));
	factor_vector::ptr labels = factor_vector::create(f, num_vertices);
	labels->set_data(set_2d_label_operate(block_size));
	vector_vector::ptr res = out_adjs->groupby(*labels,
			part_2d_apply_operate(block_size, num_vertices));

	matrix_header mheader(matrix_type::SPARSE, 0, num_vertices, num_vertices,
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
	verify_2d_matrix(mat_file, mat_idx_file);

	// Construct the transpose of the adjacency matrix.
	f_2d = fopen(t_mat_file.c_str(), "w");
	if (f_2d == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format("open %1%: %2%")
			% t_mat_file % strerror(errno);
		return;
	}
	fwrite(&mheader, sizeof(mheader), 1, f_2d);
	res = in_adjs->groupby(*labels, part_2d_apply_operate(block_size,
				num_vertices));
	ret = mem_vector::cast(res->cat())->export2(f_2d);
	assert(ret);
	fclose(f_2d);

	// Construct the index file of the adjacency matrix.
	assert(offsets.size() == res->get_num_vecs() + 1);
	off = sizeof(mheader);
	for (size_t i = 0; i < res->get_num_vecs(); i++) {
		offsets[i] = off;
		off += res->get_length(i);
	}
	offsets[res->get_num_vecs()] = off;
	mindex = SpM_2d_index::create(mheader, offsets);
	mindex->dump(t_mat_idx_file);
	verify_2d_matrix(t_mat_file, t_mat_idx_file);
}

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "el2al edge_file graph_file matrix_file\n");
		return -1;
	}

	std::string file_name = argv[1];
	std::string adj_file = std::string(argv[2]) + ".adj";
	std::string index_file = std::string(argv[2]) + ".index";
	std::string mat_name = std::string(argv[3]);
	std::vector<std::string> files;
	files.push_back(file_name);

	edge_parser parser;
	data_frame::ptr df = read_lines(files, parser);
	size_t num_edges = df->get_num_entries();
	printf("There are %ld edges\n", num_edges);
	fg::vertex_id_t max_vid = 0;
	for (size_t i = 0; i < df->get_num_vecs(); i++) {
		type_mem_vector<fg::vertex_id_t>::ptr vec
			= type_mem_vector<fg::vertex_id_t>::cast(df->get_vec(i));
		max_vid = std::max(max_vid, vec->max());
	}
	printf("max id: %d\n", max_vid);

	vector::ptr seq_vec = create_vector<fg::vertex_id_t>(0, max_vid, 1);
	vector::ptr rep_vec = create_vector<fg::vertex_id_t>(max_vid + 1,
			fg::INVALID_VERTEX_ID);
	assert(seq_vec->get_length() == rep_vec->get_length());
	// I artificially add an invalid out-edge for each vertex, so it's
	// guaranteed that each vertex exists in the adjacency lists.
	mem_data_frame::ptr new_df = mem_data_frame::create();
	new_df->add_vec(parser.get_col_name(0), seq_vec);
	new_df->add_vec(parser.get_col_name(1), rep_vec);
	df->append(new_df);

	// I artificially add an invalid in-edge for each vertex.
	new_df = mem_data_frame::create();
	new_df->add_vec(parser.get_col_name(1), seq_vec);
	new_df->add_vec(parser.get_col_name(0), rep_vec);
	df->append(new_df);

	df->sort("dest");
	assert(df->is_sorted("dest"));
	type_mem_vector<fg::vsize_t>::ptr num_in_edges
		= type_mem_vector<fg::vsize_t>::cast(df->get_vec("dest")->groupby(
			count_agg_operate(), false)->get_vec("agg"));
	adj_apply_operate in_adj_op(fg::edge_type::IN_EDGE);
	vector_vector::ptr in_adjs = df->groupby("dest", in_adj_op);
	printf("There are %ld adjacency lists and they use %ld bytes in total\n",
			in_adjs->get_num_vecs(), in_adjs->get_tot_num_entries());

	df->sort("source");
	assert(df->is_sorted("source"));
	type_mem_vector<fg::vsize_t>::ptr num_out_edges
		= type_mem_vector<fg::vsize_t>::cast(df->get_vec("source")->groupby(
			count_agg_operate(), false)->get_vec("agg"));
	adj_apply_operate out_adj_op(fg::edge_type::OUT_EDGE);
	vector_vector::ptr out_adjs = df->groupby("source", out_adj_op);
	printf("There are %ld adjacency lists and they use %ld bytes in total\n",
			out_adjs->get_num_vecs(), out_adjs->get_tot_num_entries());

	size_t num_vertices = out_adjs->get_num_vecs();
	fg::graph_header header(fg::graph_type::DIRECTED, num_vertices, num_edges, 0);

	// Construct the vertex index.
	// The vectors that contains the numbers of edges have the length of #V + 1
	// because we add -1 to the edge lists artificially and the last entries
	// are the number of vertices.
	assert(num_in_edges->get_length() - 1 == num_vertices);
	assert(num_out_edges->get_length() - 1 == num_vertices);
	fg::cdirected_vertex_index::ptr vindex
		= fg::cdirected_vertex_index::construct(num_vertices,
				(const fg::vsize_t *) num_in_edges->get_raw_arr(),
				(const fg::vsize_t *) num_out_edges->get_raw_arr(),
				header);
	vindex->dump(index_file);

	// Construct the file for the adjacency list file.
	FILE *f_graph = fopen(adj_file.c_str(), "w");
	if (f_graph == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format("open %1%: %2%")
			% adj_file % strerror(errno);
		return -1;
	}
	fwrite(&header, sizeof(header), 1, f_graph);
	bool ret = mem_vector::cast(in_adjs->cat())->export2(f_graph);
	assert(ret);
	ret = mem_vector::cast(out_adjs->cat())->export2(f_graph);
	assert(ret);
	fclose(f_graph);

	create_2d_matrix(out_adjs, in_adjs, num_vertices, mat_name);
	return 0;
}
