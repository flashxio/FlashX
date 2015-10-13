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
#include "in_mem_storage.h"

#include "fm_utils.h"
#include "factor.h"
#include "generic_type.h"
#include "local_vec_store.h"
#include "data_frame.h"
#include "vector_vector.h"
#include "local_vv_store.h"
#include "mem_vv_store.h"
#include "EM_vector.h"

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
 * This applies to a vector of values corresponding to the same key,
 * and generates an adjacency list.
 */
class adj_apply_operate: public gr_apply_operate<sub_data_frame>
{
public:
	void run(const void *key, const sub_data_frame &val,
			local_vec_store &out) const;

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

void adj_apply_operate::run(const void *key, const sub_data_frame &val,
		local_vec_store &out) const
{
	fg::vertex_id_t vid = *(const fg::vertex_id_t *) key;
	if (vid == fg::INVALID_VERTEX_ID) {
		out.resize(0);
		return;
	}

	// Right now, we just assume there aren't attributes.
	size_t edge_data_size = 0;
	assert(val.size() == 2);

	assert(out.is_type<char>());
	// The data frame is sorted based on the first vector and now we need
	// to access the entries in the second vector.
	const local_vec_store &vec = *val[1];
	assert(vec.get_type() == get_scalar_type<fg::vertex_id_t>());

	// I added an invalid edge for each vertex.
	// The invalid edge is the maximal integer.
	fg::vsize_t num_edges = vec.get_length() - 1;
	// TODO we actually don't need to alloate memory multiple times.
	std::unique_ptr<fg::vertex_id_t[]> edge_buf
		= std::unique_ptr<fg::vertex_id_t[]>(new fg::vertex_id_t[num_edges]);
	size_t edge_idx = 0;
	for (size_t i = 0; i < vec.get_length(); i++) {
		if (vec.get<fg::vertex_id_t>(i) == fg::INVALID_VERTEX_ID
				// skip self-edges.
				|| (remove_selfe && vec.get<fg::vertex_id_t>(i) == vid))
			continue;
		edge_buf[edge_idx++] = vec.get<fg::vertex_id_t>(i);
	}
	assert(edge_idx <= num_edges);
	// If there are self-edges, edge_idx has the actual number of edges.
	num_edges = edge_idx;
	std::sort(edge_buf.get(), edge_buf.get() + num_edges);
	if (deduplicate) {
		fg::vertex_id_t *end = std::unique(edge_buf.get(),
				edge_buf.get() + num_edges);
		num_edges = end - edge_buf.get();
	}
	size_t size = fg::ext_mem_undirected_vertex::num_edges2vsize(num_edges,
			edge_data_size);
	out.resize(size);

	// Even if we generate a directed, we still can use undirected vertex to
	// store one type of edges of a vertex.
	fg::in_mem_undirected_vertex<> v(vid, edge_data_size > 0);
	for (size_t i = 0; i < num_edges; i++)
		v.add_edge(fg::edge<>(vid, edge_buf[i]));
	fg::ext_mem_undirected_vertex::serialize(v, out.get_raw_arr(), size,
			// The edge type here actually doesn't matter since it's
			// an undirected vertex.
			fg::edge_type::OUT_EDGE);
}

template<class AttrType>
class attr_adj_apply_operate: public gr_apply_operate<sub_data_frame>
{
	typedef std::pair<fg::vertex_id_t, AttrType> edge_type;

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
	void run(const void *key, const sub_data_frame &val,
			local_vec_store &out) const;

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

template<class AttrType>
void attr_adj_apply_operate<AttrType>::run(const void *key,
		const sub_data_frame &val, local_vec_store &out) const
{
	fg::vertex_id_t vid = *(const fg::vertex_id_t *) key;
	if (vid == fg::INVALID_VERTEX_ID) {
		out.resize(0);
		return;
	}

	assert(val.size() == 3);
	assert(out.is_type<char>());
	// The data frame is sorted based on the first vector and now we need
	// to access the entries in the second vector.
	const local_vec_store &vec = *val[1];
	const local_vec_store &attr_vec = *val[2];
	assert(vec.get_type() == get_scalar_type<fg::vertex_id_t>());

	// I added an invalid edge for each vertex.
	// The invalid edge is the maximal integer.
	fg::vsize_t num_edges = vec.get_length() - 1;
	// TODO we actually don't need to alloate memory multiple times.
	std::unique_ptr<edge_type[]> edge_buf
		= std::unique_ptr<edge_type[]>(new edge_type[num_edges]);
	size_t edge_idx = 0;
	for (size_t i = 0; i < vec.get_length(); i++) {
		if (vec.get<fg::vertex_id_t>(i) == fg::INVALID_VERTEX_ID
				// skip self-edges.
				|| (remove_selfe && vec.get<fg::vertex_id_t>(i) == vid))
			continue;
		edge_buf[edge_idx].first = vec.get<fg::vertex_id_t>(i);
		edge_buf[edge_idx].second = attr_vec.get<AttrType>(i);
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
	size_t edge_data_size = val[2]->get_entry_size();
	size_t size = fg::ext_mem_undirected_vertex::num_edges2vsize(num_edges,
			edge_data_size);
	out.resize(size);

	// Even if we generate a directed, we still can use undirected vertex to
	// store one type of edges of a vertex.
	fg::in_mem_undirected_vertex<AttrType> v(vid, edge_data_size > 0);
	for (size_t i = 0; i < num_edges; i++)
		v.add_edge(fg::edge<AttrType>(vid, edge_buf[i].first,
					edge_buf[i].second));
	fg::ext_mem_undirected_vertex::serialize(v, out.get_raw_arr(), size,
			// The edge type here actually doesn't matter since it's
			// an undirected vertex.
			fg::edge_type::OUT_EDGE);
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

}

vector_vector::ptr create_1d_matrix(data_frame::ptr df)
{
	if (df->get_num_vecs() < 2) {
		BOOST_LOG_TRIVIAL(error)
			<< "The data frame needs to contain at least 2 vectors";
		return vector_vector::ptr();
	}

	struct timeval start, end;
	gettimeofday(&start, NULL);
	std::string sort_vec_name = df->get_vec_name(0);
	data_frame::const_ptr sorted_df = df->sort(sort_vec_name);
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to sort the edge list\n",
			time_diff(start, end));
	gettimeofday(&start, NULL);
	assert(sorted_df->is_sorted(sort_vec_name));
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to test if the edge list is sorted\n",
			time_diff(start, end));

	gettimeofday(&start, NULL);
	vector_vector::ptr ret;
	if (df->get_num_vecs() == 2)
		ret = sorted_df->groupby(sort_vec_name, adj_apply_operate());
	// Instead of giving the real data type, we give a type that indicates
	// the size of the edge data size. Actually, we don't interpret data type
	// here. Only the data size matters.
	else if (df->get_vec(2)->get_entry_size() == 4)
		ret = sorted_df->groupby(sort_vec_name, attr_adj_apply_operate<unit4>());
	else if (df->get_vec(2)->get_entry_size() == 8)
		ret = sorted_df->groupby(sort_vec_name, attr_adj_apply_operate<unit8>());
	else {
		BOOST_LOG_TRIVIAL(error)
			<< "The edge attribute has an unsupported type";
		return vector_vector::ptr();
	}
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to groupby the edge list.\n",
			time_diff(start, end));
	return ret;
}

static std::pair<fg::vertex_index::ptr, detail::vec_store::ptr> create_fg_directed_graph(
		const std::string &graph_name, data_frame::ptr df)
{
	struct timeval start, end;
	// Leave the space for graph header.
	detail::vec_store::ptr graph_data = detail::vec_store::create(
			fg::graph_header::get_header_size(),
			get_scalar_type<char>(), df->get_vec(0)->is_in_mem());
	size_t edge_data_size = 0;
	if (df->get_num_vecs() == 3) {
		auto vec = df->get_vec("attr");
		assert(vec);
		edge_data_size = vec->get_entry_size();
	}

	/*
	 * Construct the in-edge adjacency lists.
	 * All edges share the same destination vertex should be stored together.
	 */

	data_frame::ptr tmp = data_frame::create();
	tmp->add_vec("dest", df->get_vec("dest"));
	tmp->add_vec("source", df->get_vec("source"));
	if (df->get_num_vecs() == 3) {
		auto vec = df->get_vec("attr");
		assert(vec);
		tmp->add_vec("attr", vec);
	}
	df = tmp;
	vector_vector::ptr in_adjs = fm::create_1d_matrix(df);
	size_t num_vertices = in_adjs->get_num_vecs();
	printf("There are %ld in-edge adjacency lists and they use %ld bytes in total\n",
			in_adjs->get_num_vecs(), in_adjs->get_tot_num_entries());
	gettimeofday(&start, NULL);
	// Get the number of in-edges for each vertex.
	detail::smp_vec_store::ptr num_in_edges = detail::smp_vec_store::create(
			num_vertices, get_scalar_type<fg::vsize_t>());
	for (size_t i = 0; i < num_vertices; i++) {
		num_in_edges->set<fg::vsize_t>(i,
				fg::ext_mem_undirected_vertex::vsize2num_edges(
					in_adjs->get_length(i), edge_data_size));
	}
	size_t num_edges = vector::create(num_in_edges)->sum<fg::vsize_t>();
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to get #in-edges\n", time_diff(start, end));
	// Move in-edge adjacency lists to the final image.
	const detail::vv_store &in_adj_store
		= dynamic_cast<const detail::vv_store &>(in_adjs->get_data());
	gettimeofday(&start, NULL);
	graph_data->append(in_adj_store.get_data());
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to append in-edge adjacency list\n",
			time_diff(start, end));
	in_adjs = NULL;

	/*
	 * Construct out-edge adjacency lists.
	 * All edges share the same source vertex should be stored together.
	 */

	tmp = data_frame::create();
	tmp->add_vec("source", df->get_vec("source"));
	tmp->add_vec("dest", df->get_vec("dest"));
	if (df->get_num_vecs() == 3) {
		auto vec = df->get_vec("attr");
		assert(vec);
		tmp->add_vec("attr", vec);
	}
	df = tmp;
	vector_vector::ptr out_adjs = create_1d_matrix(df);
	printf("There are %ld out-edge adjacency lists and they use %ld bytes in total\n",
			out_adjs->get_num_vecs(), out_adjs->get_tot_num_entries());
	assert(out_adjs->get_num_vecs() == num_vertices);
	// Get the number of out-edge for each vertex.
	gettimeofday(&start, NULL);
	detail::smp_vec_store::ptr num_out_edges = detail::smp_vec_store::create(
			num_vertices, get_scalar_type<fg::vsize_t>());
	for (size_t i = 0; i < num_vertices; i++) {
		num_out_edges->set<fg::vsize_t>(i,
				fg::ext_mem_undirected_vertex::vsize2num_edges(
					out_adjs->get_length(i), edge_data_size));
	}
	assert(vector::create(num_out_edges)->sum<fg::vsize_t>() == num_edges);
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to get #out-edges\n", time_diff(start, end));
	printf("There are %ld edges\n", num_edges);
	// Move out-edge adjacency lists to the final image.
	const detail::vv_store &out_adj_store
		= dynamic_cast<const detail::vv_store &>(out_adjs->get_data());
	gettimeofday(&start, NULL);
	graph_data->append(out_adj_store.get_data());
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to append out-edge adjacency list\n",
			time_diff(start, end));
	out_adjs = NULL;

	// Construct the graph header.
	gettimeofday(&start, NULL);
	fg::graph_header header(fg::graph_type::DIRECTED, num_vertices, num_edges,
			edge_data_size);
	local_vec_store::ptr header_store(new local_buf_vec_store(0,
			fg::graph_header::get_header_size(), get_scalar_type<char>(), -1));
	memcpy(header_store->get_raw_arr(), &header,
			fg::graph_header::get_header_size());
	graph_data->set_portion(header_store, 0);
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to construct the graph header\n",
			time_diff(start, end));

	// Construct the vertex index.
	// The vectors that contains the numbers of edges have the length of #V + 1
	// because we add -1 to the edge lists artificially and the last entries
	// are the number of vertices.
	printf("create the vertex index image\n");
	gettimeofday(&start, NULL);
	fg::cdirected_vertex_index::ptr vindex
		= fg::cdirected_vertex_index::construct(num_vertices,
				(const fg::vsize_t *) num_in_edges->get_raw_arr(),
				(const fg::vsize_t *) num_out_edges->get_raw_arr(),
				header);
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to construct the graph index\n",
			time_diff(start, end));
	return std::pair<fg::vertex_index::ptr, detail::vec_store::ptr>(vindex,
			graph_data);
}

static std::pair<fg::vertex_index::ptr, detail::vec_store::ptr> create_fg_undirected_graph(
		const std::string &graph_name, data_frame::ptr df)
{
	struct timeval start, end;
	// Leave the space for graph header.
	detail::vec_store::ptr graph_data = detail::vec_store::create(0,
			get_scalar_type<char>(), df->get_vec(0)->is_in_mem());
	size_t edge_data_size = 0;
	if (df->get_num_vecs() == 3) {
		auto vec = df->get_vec("attr");
		assert(vec);
		edge_data_size = vec->get_entry_size();
	}

	vector_vector::ptr adjs = create_1d_matrix(df);
	printf("There are %ld vertices and they use %ld bytes in total\n",
			adjs->get_num_vecs(), adjs->get_tot_num_entries());

	gettimeofday(&start, NULL);
	size_t num_vertices = adjs->get_num_vecs();
	detail::smp_vec_store::ptr num_out_edges = detail::smp_vec_store::create(
			num_vertices, get_scalar_type<fg::vsize_t>());
	size_t num_edges = 0;
	for (size_t i = 0; i < num_vertices; i++) {
		size_t local_num_edges = fg::ext_mem_undirected_vertex::vsize2num_edges(
					adjs->get_length(i), edge_data_size);
		num_out_edges->set<fg::vsize_t>(i, local_num_edges);
		num_edges += local_num_edges;
	}
	assert(num_edges % 2 == 0);
	num_edges /= 2;
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to get #edges\n", time_diff(start, end));
	printf("There are %ld edges\n", num_edges);

	printf("create the graph image\n");
	gettimeofday(&start, NULL);
	fg::graph_header header(fg::graph_type::UNDIRECTED, num_vertices, num_edges,
			edge_data_size);
	local_vec_store::ptr header_store(new local_buf_vec_store(0,
			fg::graph_header::get_header_size(), get_scalar_type<char>(), -1));
	memcpy(header_store->get_raw_arr(), &header,
			fg::graph_header::get_header_size());
	graph_data->append(*header_store);

	const detail::vv_store &adj_store
		= dynamic_cast<const detail::vv_store &>(adjs->get_data());
	graph_data->append(adj_store.get_data());
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to append the adjacency list\n",
			time_diff(start, end));

	// Construct the vertex index.
	// The vectors that contains the numbers of edges have the length of #V + 1
	// because we add -1 to the edge lists artificially and the last entries
	// are the number of vertices.
	printf("create the vertex index image\n");
	gettimeofday(&start, NULL);
	fg::cundirected_vertex_index::ptr vindex
		= fg::cundirected_vertex_index::construct(num_vertices,
				(const fg::vsize_t *) num_out_edges->get_raw_arr(), header);
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to construct the graph index\n",
			time_diff(start, end));
	return std::pair<fg::vertex_index::ptr, detail::vec_store::ptr>(vindex,
			graph_data);
}

fg::FG_graph::ptr create_fg_graph(const std::string &graph_name,
		data_frame::ptr df, bool directed)
{
	std::pair<fg::vertex_index::ptr, detail::vec_store::ptr> res;
	if (directed)
		res = create_fg_directed_graph(graph_name, df);
	else
		res = create_fg_undirected_graph(graph_name, df);

	if (res.second->is_in_mem()) {
		fg::in_mem_graph::ptr graph = fg::in_mem_graph::create(graph_name,
				detail::smp_vec_store::cast(res.second)->get_raw_data(),
				res.second->get_length());
		return fg::FG_graph::create(graph, res.first, graph_name, NULL);
	}
	else {
		detail::EM_vec_store::ptr graph_data = detail::EM_vec_store::cast(res.second);
		detail::EM_vec_store::ptr index_vec = detail::EM_vec_store::create(0,
				get_scalar_type<char>());
		local_cref_vec_store index_store((const char *) res.first.get(), 0,
					res.first->get_index_size(), get_scalar_type<char>(), -1);
		index_vec->append(index_store);
		std::string graph_file_name = graph_name + ".adj";
		bool ret = graph_data->set_persistent(graph_file_name);
		assert(ret);
		std::string index_file_name = graph_name + ".index";
		ret = index_vec->set_persistent(index_file_name);
		assert(ret);
		return fg::FG_graph::create(graph_file_name, index_file_name, NULL);
	}
}

class set_2d_label_operate: public type_set_vec_operate<factor_value_t>
{
	block_2d_size block_size;
public:
	set_2d_label_operate(const block_2d_size &_size): block_size(_size) {
	}

	virtual void set(factor_value_t *arr, size_t num_eles, off_t start_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = (start_idx + i) / block_size.get_num_rows();
	}
};

class part_2d_apply_operate: public gr_apply_operate<local_vv_store>
{
	// The row length (aka. the total number of columns) of the matrix.
	size_t row_len;
	size_t nz_size;
	block_2d_size block_size;
public:
	part_2d_apply_operate(const block_2d_size &_size,
			size_t row_len, size_t nz_size): block_size(_size) {
		this->row_len = row_len;
		this->nz_size = nz_size;
	}

	void run(const void *key, const local_vv_store &val,
			local_vec_store &out) const;

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

namespace
{

struct coo_less
{
	bool operator()(const coo_nz_t &a, const coo_nz_t &b) {
		if (a.first == b.first)
			return a.second < b.second;
		else
			return a.first < b.first;
	}
};

class block_nz_data
{
	size_t entry_size;
	std::vector<char> data;
	size_t num_entries;
public:
	block_nz_data(size_t entry_size) {
		this->entry_size = entry_size;
		num_entries = 0;
	}

	void clear() {
		num_entries = 0;
		data.clear();
	}

	void push_back(char *new_entry) {
		assert(entry_size > 0);
		if (data.empty())
			data.resize(16 * entry_size);
		// If full
		else if (data.size() == entry_size * num_entries)
			data.resize(data.size() * 2);
		memcpy(&data[num_entries * entry_size], new_entry, entry_size);
		num_entries++;
	}

	void append(const char *new_entries, size_t num) {
		assert(entry_size > 0);
		if (data.empty())
			data.resize(num_entries * entry_size);
		else if (data.size() < (this->num_entries + num) * entry_size)
			data.resize((this->num_entries + num) * entry_size);
		memcpy(&data[num_entries * entry_size], new_entries, entry_size * num);
		num_entries += num;
	}

	const char *get_data() const {
		if (entry_size == 0)
			return NULL;
		else
			return data.data();
	}

	size_t get_size() const {
		return num_entries;
	}

	bool is_empty() const {
		if (entry_size == 0)
			return true;
		else
			return data.empty();
	}
};

}

void part_2d_apply_operate::run(const void *key, const local_vv_store &val,
		local_vec_store &out) const
{
	size_t block_height = block_size.get_num_rows();
	size_t block_width = block_size.get_num_cols();
	size_t num_blocks = ceil(((double) row_len) / block_width);
	factor_value_t block_row_id = *(const factor_value_t *) key;
	size_t tot_num_non_zeros = 0;
	size_t max_row_parts = 0;
	const fg::ext_mem_undirected_vertex *first_v
		= (const fg::ext_mem_undirected_vertex *) val.get_raw_arr(0);
	fg::vertex_id_t start_vid = first_v->get_id();
	std::vector<std::vector<const fg::ext_mem_undirected_vertex *> > edge_dist_map(
			num_blocks);
	for (size_t i = 0; i < val.get_num_vecs(); i++) {
		const fg::ext_mem_undirected_vertex *v
			= (const fg::ext_mem_undirected_vertex *) val.get_raw_arr(i);
		assert(val.get_length(i) == v->get_size());
		assert(v->get_id() / block_height == (size_t) block_row_id);
		tot_num_non_zeros += v->get_num_edges();
		// I definitely over estimate the number of row parts.
		// If a row doesn't have many non-zero entries, I assume that
		// the non-zero entries distribute evenly across all row parts.
		max_row_parts += std::min(num_blocks, v->get_num_edges());

		// Fill the edge distribution map.
		for (size_t i = 0; i < v->get_num_edges(); i++) {
			size_t vector_idx = v->get_neighbor(i) / block_width;
			if (edge_dist_map[vector_idx].empty()
					|| edge_dist_map[vector_idx].back() != v)
				edge_dist_map[vector_idx].push_back(v);
		}
	}

	// Containers of non-zero values.
	block_nz_data data(nz_size);
	block_nz_data single_nz_data(nz_size);

	std::vector<size_t> neigh_idxs(val.get_num_vecs());
	// The maximal size of a block.
	size_t max_block_size
		// Even if a block is empty, its header still exists. The size is
		// accurate.
		= sizeof(sparse_block_2d) * num_blocks
		// Each block has an empty row part in the end to indicate the end
		// of the block.
		+ sparse_row_part::get_size(0) * num_blocks
		// The size for row part headers is highly over estimated.
		+ sizeof(sparse_row_part) * max_row_parts
		// The size is accurate.
		+ sparse_row_part::get_col_entry_size() * tot_num_non_zeros
		// The size of the non-zero entries.
		+ nz_size * tot_num_non_zeros;
	size_t max_idx_size = max_block_size - nz_size * tot_num_non_zeros;
	out.resize(max_block_size);
	size_t curr_size = 0;
	// The maximal size of a row part.
	size_t max_row_size = sparse_row_part::get_size(block_width);
	std::unique_ptr<char[]> buf = std::unique_ptr<char[]>(new char[max_row_size]);
	std::vector<fg::vertex_id_t> local_neighs;
	size_t num_non_zeros = 0;
	// Iterate columns. Actually it strides instead of iterating all columns.
	for (size_t col_idx = 0; col_idx < row_len; col_idx += block_width) {
		sparse_block_2d *block
			= new (out.get_raw_arr() + curr_size) sparse_block_2d(
					block_row_id, col_idx / block_width);
		data.clear();
		single_nz_data.clear();
		const std::vector<const fg::ext_mem_undirected_vertex *> &v_ptrs
			= edge_dist_map[col_idx / block_width];
		std::vector<coo_nz_t> single_nnz;
		// Iterate the vectors in the vector_vector one by one.
		for (size_t i = 0; i < v_ptrs.size(); i++) {
			const fg::ext_mem_undirected_vertex *v = v_ptrs[i];
			fg::vertex_id_t row_idx = v->get_id() - start_vid;
			assert(v->get_edge_data_size() == nz_size);
			assert(neigh_idxs[row_idx] < v->get_num_edges());
			assert(v->get_neighbor(neigh_idxs[row_idx]) >= col_idx
					&& v->get_neighbor(neigh_idxs[row_idx]) < col_idx + block_width);

			size_t idx = neigh_idxs[row_idx];
			local_neighs.clear();
			for (; idx < v->get_num_edges()
					&& v->get_neighbor(idx) < col_idx + block_width; idx++)
				local_neighs.push_back(v->get_neighbor(idx));
			// Get all edge data.
			if (nz_size > 0 && local_neighs.size() == 1)
				single_nz_data.push_back(v->get_raw_edge_data(neigh_idxs[row_idx]));
			else if (nz_size > 0) {
				size_t tmp_idx = neigh_idxs[row_idx];
				for (; tmp_idx < v->get_num_edges()
						&& v->get_neighbor(tmp_idx) < col_idx + block_width;
						tmp_idx++)
					data.push_back(v->get_raw_edge_data(tmp_idx));
				assert(tmp_idx == idx);
			}
			size_t local_nnz = local_neighs.size();
			assert(local_nnz <= block_width);
			neigh_idxs[row_idx] = idx;

			if (local_neighs.size() > 1) {
				sparse_row_part *part = new (buf.get()) sparse_row_part(row_idx);
				rp_edge_iterator edge_it = part->get_edge_iterator();
				for (size_t k = 0; k < local_neighs.size(); k++)
					edge_it.append(block_size, local_neighs[k]);
				// These two parts don't count the space used by non-zero entries.
				assert(block->get_size(0) + sparse_row_part::get_size(local_nnz)
						// So is here.
						<= max_idx_size - curr_size);
				block->append(*part, sparse_row_part::get_size(local_nnz));
			}
			else if (local_neighs.size() == 1)
				single_nnz.push_back(coo_nz_t(row_idx, local_neighs[0]));
		}
		if (!single_nnz.empty())
			block->add_coo(single_nnz, block_size);
		if (!single_nz_data.is_empty())
			data.append(single_nz_data.get_data(), single_nz_data.get_size());
		// After we finish adding rows to a block, we need to finalize it.
		block->finalize(data.get_data(), data.get_size() * nz_size);
		if (!block->is_empty()) {
			curr_size += block->get_size(nz_size);
			block->verify(block_size);
		}
		num_non_zeros += block->get_nnz();
	}
	assert(tot_num_non_zeros == num_non_zeros);
	// If the block row doesn't have any non-zero entries, let's insert an
	// empty block row, so the matrix index can work correctly.
	if (curr_size == 0) {
		curr_size = sizeof(sparse_block_2d);
		const sparse_block_2d *block
			= (const sparse_block_2d *) out.get_raw_arr();
		assert(block->is_empty());
	}
	out.resize(curr_size);

#ifdef VERIFY_ENABLED
	std::vector<coo_nz_t> all_coos;
	block_row_iterator br_it((const sparse_block_2d *) out.get_raw_arr(),
			(const sparse_block_2d *) (out.get_raw_arr() + curr_size));
	while (br_it.has_next()) {
		const sparse_block_2d &block = br_it.next();
		std::vector<coo_nz_t> coos = block.get_non_zeros(block_size);
		all_coos.insert(all_coos.end(), coos.begin(), coos.end());
	}
	assert(tot_num_non_zeros == all_coos.size());
	std::sort(all_coos.begin(), all_coos.end(), coo_less());

	auto coo_it = all_coos.begin();
	for (size_t i = 0; i < val.get_num_vecs(); i++) {
		const fg::ext_mem_undirected_vertex *v
			= (const fg::ext_mem_undirected_vertex *) val.get_raw_arr(i);
		bool correct = true;
		for (size_t j = 0; j < v->get_num_edges(); j++) {
			if (v->get_id() != coo_it->first)
				correct = false;
			if (v->get_neighbor(j) != coo_it->second)
				correct = false;
			coo_it++;
		}
		assert(correct);
	}
#endif
}

std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr> create_2d_matrix(
		vector_vector::ptr adjs, const block_2d_size &block_size,
		const scalar_type *entry_type)
{
	size_t entry_size = 0;
	if (entry_type)
		entry_size = entry_type->get_size();
	size_t num_rows = adjs->get_num_vecs();
	factor f(ceil(((double) num_rows) / block_size.get_num_rows()));
	factor_vector::ptr labels = factor_vector::create(f, num_rows,
			adjs->is_in_mem(), set_2d_label_operate(block_size));
	printf("groupby multiple vectors in the vector vector\n");
	vector_vector::ptr res = adjs->groupby(*labels,
			part_2d_apply_operate(block_size, num_rows, entry_size));

	prim_type type = prim_type::P_BOOL;
	if (entry_type)
		type = entry_type->get_type();
	matrix_header mheader(matrix_type::SPARSE, entry_size, num_rows, num_rows,
			matrix_layout_t::L_ROW_2D, type, block_size);

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

std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr> create_2d_matrix(
		data_frame::ptr df, const block_2d_size &block_size,
		const scalar_type *entry_type)
{
	vector_vector::ptr vv = create_1d_matrix(df);
	if (vv == NULL)
		return std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr>();
	else
		return create_2d_matrix(vv, block_size, entry_type);
}

void export_2d_matrix(vector_vector::ptr adjs, const block_2d_size &block_size,
		const scalar_type *entry_type, const std::string &mat_file,
		const std::string &mat_idx_file, bool to_safs)
{
	size_t entry_size = 0;
	if (entry_type)
		entry_size = entry_type->get_size();
	size_t num_rows = adjs->get_num_vecs();
	factor f(ceil(((double) num_rows) / block_size.get_num_rows()));
	factor_vector::ptr labels = factor_vector::create(f, num_rows,
			adjs->is_in_mem(), set_2d_label_operate(block_size));
	vector_vector::ptr res = adjs->groupby(*labels,
			part_2d_apply_operate(block_size, num_rows, entry_size));

	prim_type type = prim_type::P_BOOL;
	if (entry_type)
		type = entry_type->get_type();
	matrix_header mheader(matrix_type::SPARSE, entry_size, num_rows, num_rows,
			matrix_layout_t::L_ROW_2D, type, block_size);
	if (!to_safs) {
		FILE *f_2d = fopen(mat_file.c_str(), "w");
		if (f_2d == NULL) {
			BOOST_LOG_TRIVIAL(error) << boost::format("open %1%: %2%")
				% mat_file % strerror(errno);
			return;
		}
		fwrite(&mheader, sizeof(mheader), 1, f_2d);
		bool ret = res->cat()->export2(f_2d);
		assert(ret);
		fclose(f_2d);
	}
	else {
		detail::EM_vec_store::ptr vec = detail::EM_vec_store::create(0,
				get_scalar_type<char>());
		local_cref_vec_store header_store((const char *) &mheader,
				0, sizeof(mheader), get_scalar_type<char>(), -1);
		vec->append(header_store);
		vec->append(dynamic_cast<const detail::vv_store &>(
					res->get_data()).get_data());
		vec->set_persistent(mat_file);
	}

	// Construct the index file of the adjacency matrix.
	std::vector<off_t> offsets(res->get_num_vecs() + 1);
	off_t off = sizeof(mheader);
	for (size_t i = 0; i < res->get_num_vecs(); i++) {
		offsets[i] = off;
		off += res->get_length(i);
	}
	offsets[res->get_num_vecs()] = off;
	SpM_2d_index::ptr mindex = SpM_2d_index::create(mheader, offsets);
	if (!to_safs)
		mindex->dump(mat_idx_file);
	else
		mindex->safs_dump(mat_idx_file);
}

}
