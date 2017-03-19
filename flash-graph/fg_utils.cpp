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

#include "fg_utils.h"
#include "in_mem_storage.h"

#include "factor.h"
#include "generic_type.h"
#include "local_vec_store.h"
#include "data_frame.h"
#include "vector_vector.h"
#include "local_vv_store.h"
#include "mem_vv_store.h"
#include "EM_vv_store.h"
#include "EM_vector.h"
#include "fg_utils.h"
#include "data_io.h"

using namespace fm;

namespace fg
{

size_t get_out_size(fg::vertex_index::ptr vindex)
{
	if (vindex->is_compressed()) {
		fg::vsize_t num_vertices = vindex->get_num_vertices();
		if (vindex->get_graph_header().is_directed_graph()) {
			fg::in_mem_cdirected_vertex_index::ptr dindex
				= fg::in_mem_cdirected_vertex_index::create(*vindex);
			fg::directed_vertex_entry dentry = dindex->get_vertex(num_vertices - 1);
			return dentry.get_out_off() + dindex->get_out_size(
					num_vertices - 1) - vindex->get_out_part_loc();
		}
		else {
			fg::in_mem_cundirected_vertex_index::ptr uindex
				= fg::in_mem_cundirected_vertex_index::create(*vindex);
			fg::vertex_offset off = uindex->get_vertex(num_vertices - 1);
			return off.get_off() + uindex->get_size(
					num_vertices - 1) - vindex->get_header_size();
		}
	}
	else {
		if (vindex->get_graph_header().is_directed_graph()) {
			fg::directed_vertex_index::ptr dindex
				= fg::directed_vertex_index::cast(vindex);
			return dindex->get_graph_size() - vindex->get_out_part_loc();
		}
		else {
			fg::undirected_vertex_index::ptr uindex
				= fg::undirected_vertex_index::cast(vindex);
			return uindex->get_graph_size() - vindex->get_header_size();
		}
	}
}

void init_out_offs(fg::vertex_index::ptr vindex, std::vector<off_t> &out_offs)
{
	size_t num_vertices = vindex->get_num_vertices();
	assert(num_vertices + 1 == out_offs.size());
	if (vindex->is_compressed()) {
		out_offs[0] = get_out_off(vindex);
		if (vindex->get_graph_header().is_directed_graph()) {
			fg::in_mem_cdirected_vertex_index::ptr dindex
				= fg::in_mem_cdirected_vertex_index::create(*vindex);
			for (size_t i = 1; i <= num_vertices; i++)
				out_offs[i] = out_offs[i - 1] + dindex->get_out_size(i - 1);
		}
		else {
			fg::in_mem_cundirected_vertex_index::ptr uindex
				= fg::in_mem_cundirected_vertex_index::create(*vindex);
			for (size_t i = 1; i <= num_vertices; i++)
				out_offs[i] = out_offs[i - 1] + uindex->get_size(i - 1);
		}
		assert((size_t) out_offs[num_vertices]
				== get_out_size(vindex) + out_offs[0]);
	}
	else {
		if (vindex->get_graph_header().is_directed_graph()) {
			off_t out_part_loc = vindex->get_out_part_loc();
			fg::directed_vertex_index::ptr dindex
				= fg::directed_vertex_index::cast(vindex);
			for (size_t i = 0; i < num_vertices; i++)
				out_offs[i] = dindex->get_vertex(i).get_out_off();
			out_offs[num_vertices] = get_out_size(vindex) + out_part_loc;
		}
		else {
			fg::undirected_vertex_index::ptr uindex
				= fg::undirected_vertex_index::cast(vindex);
			for (size_t i = 0; i < num_vertices; i++)
				out_offs[i] = uindex->get_vertex(i).get_off();
			out_offs[num_vertices]
				= get_out_size(vindex) + vindex->get_header_size();
		}
	}
	for (size_t i = 1; i <= num_vertices; i++)
		assert(out_offs[i] > out_offs[i - 1]);
}

void init_in_offs(fg::vertex_index::ptr vindex, std::vector<off_t> &in_offs)
{
	size_t num_vertices = vindex->get_num_vertices();
	assert(num_vertices + 1 == in_offs.size());
	assert(vindex->get_graph_header().is_directed_graph());
	if (vindex->is_compressed()) {
		fg::in_mem_cdirected_vertex_index::ptr dindex
			= fg::in_mem_cdirected_vertex_index::create(*vindex);
		in_offs[0] = get_in_off(vindex);
		for (size_t i = 1; i <= num_vertices; i++)
			in_offs[i] = in_offs[i - 1] + dindex->get_in_size(i - 1);
		assert((size_t) in_offs[num_vertices]
				== get_in_size(vindex) + in_offs[0]);
	}
	else {
		fg::directed_vertex_index::ptr dindex
			= fg::directed_vertex_index::cast(vindex);
		for (size_t i = 0; i < num_vertices; i++)
			in_offs[i] = dindex->get_vertex(i).get_in_off();
		in_offs[num_vertices] = get_in_size(vindex) + vindex->get_header_size();
	}
	for (size_t i = 1; i <= num_vertices; i++)
		assert(in_offs[i] > in_offs[i - 1]);
}

edge_list::ptr edge_list::create(data_frame::ptr df, bool directed)
{
	if (df->get_num_vecs() < 2) {
		BOOST_LOG_TRIVIAL(error)
			<< "The data frame needs to contain at least 2 vectors";
		return edge_list::ptr();
	}

	fg::vertex_id_t max_vid;
	{
		vector::ptr vec = vector::create(df->get_vec(0));
		assert(vec->get_type() == get_scalar_type<fg::vertex_id_t>());
		max_vid = vec->max<fg::vertex_id_t>();
		vec = vector::create(df->get_vec(1));
		assert(vec->get_type() == get_scalar_type<fg::vertex_id_t>());
		max_vid = std::max(max_vid, vec->max<fg::vertex_id_t>());
	}
	printf("max id: %d\n", max_vid);

	detail::vec_store::ptr seq_vec = detail::create_seq_vec_store<fg::vertex_id_t>(
			0, max_vid, 1);
	detail::vec_store::ptr rep_vec = detail::create_rep_vec_store<fg::vertex_id_t>(
			max_vid + 1, fg::INVALID_VERTEX_ID);
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
		return edge_list::ptr();
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
	return edge_list::ptr(new edge_list(df, directed));
}

edge_list::ptr edge_list::sort_source() const
{
	std::string sort_name = df->get_vec_name(0);
	return edge_list::ptr(new edge_list(df->sort(sort_name), directed));
}

vector_vector::ptr edge_list::groupby_source(
		const gr_apply_operate<sub_data_frame> &op) const
{
	std::string name = df->get_vec_name(0);
	return df->groupby(name, op);
}

size_t edge_list::get_attr_size() const
{
	if (get_num_vecs() > 2) {
		auto vec = df->get_vec(2);
		if (vec)
			return vec->get_entry_size();
		else
			return 0;
	}
	else
		return 0;
}

edge_list::ptr edge_list::reverse_edge() const
{
	std::vector<off_t> vec_idxs(df->get_num_vecs());
	vec_idxs[0] = 1;
	vec_idxs[1] = 0;
	for (size_t i = 2; i < vec_idxs.size(); i++)
		vec_idxs[i] = i;
	data_frame::const_ptr new_df = df->shuffle_vecs(vec_idxs);
	return edge_list::ptr(new edge_list(new_df, directed));
}

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
	std::vector<size_t> max_col_idxs;
public:
	adj_apply_operate() {
		int num_threads = detail::mem_thread_pool::get_global_num_threads();
		max_col_idxs.resize(num_threads);
	}

	size_t get_max_col_idx() const {
		size_t max = max_col_idxs[0];
		for (size_t i = 1; i < max_col_idxs.size(); i++)
			max = std::max(max, max_col_idxs[i]);
		return max;
	}

	virtual bool ignore_key(const void *key) const {
		fg::vertex_id_t vid = *(const fg::vertex_id_t *) key;
		return vid == fg::INVALID_VERTEX_ID;
	}

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
	size_t max_col_idx = 0;
	for (size_t i = 0; i < vec.get_length(); i++) {
		if (vec.get<fg::vertex_id_t>(i) != fg::INVALID_VERTEX_ID)
			max_col_idx = std::max(max_col_idx,
					(size_t) vec.get<fg::vertex_id_t>(i));
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

	// Here is the max column index I have found so far.
	int thread_id = detail::mem_thread_pool::get_curr_thread_id();
	adj_apply_operate *mutable_this = const_cast<adj_apply_operate *>(this);
	mutable_this->max_col_idxs[thread_id] = std::max(max_col_idx,
			mutable_this->max_col_idxs[thread_id]);

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
	std::vector<size_t> max_col_idxs;
public:
	attr_adj_apply_operate() {
		int num_threads = detail::mem_thread_pool::get_global_num_threads();
		max_col_idxs.resize(num_threads);
	}

	size_t get_max_col_idx() const {
		size_t max = max_col_idxs[0];
		for (size_t i = 1; i < max_col_idxs.size(); i++)
			max = std::max(max, max_col_idxs[i]);
		return max;
	}

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
	size_t max_col_idx = 0;
	for (size_t i = 0; i < vec.get_length(); i++) {
		if (vec.get<fg::vertex_id_t>(i) != fg::INVALID_VERTEX_ID)
			max_col_idx = std::max(max_col_idx,
					(size_t) vec.get<fg::vertex_id_t>(i));
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

	// Here is the max column index I have found so far.
	int thread_id = detail::mem_thread_pool::get_curr_thread_id();
	attr_adj_apply_operate<AttrType> *mutable_this
		= const_cast<attr_adj_apply_operate<AttrType> *>(this);
	mutable_this->max_col_idxs[thread_id] = std::max(max_col_idx,
			mutable_this->max_col_idxs[thread_id]);

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

std::pair<vector_vector::ptr, size_t> create_1d_matrix(edge_list::ptr el)
{
	struct timeval start, end;
	gettimeofday(&start, NULL);
	edge_list::const_ptr sorted_el = el->sort_source();
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to sort the edge list\n",
			time_diff(start, end));

	gettimeofday(&start, NULL);
	vector_vector::ptr ret;
	size_t max_col_idx = 0;
	if (!el->has_attr()) {
		std::unique_ptr<adj_apply_operate> op(new adj_apply_operate());
		ret = sorted_el->groupby_source(*op);
		max_col_idx = op->get_max_col_idx();
	}
	// Instead of giving the real data type, we give a type that indicates
	// the size of the edge data size. Actually, we don't interpret data type
	// here. Only the data size matters.
	else if (el->get_attr_size() == 4) {
		std::unique_ptr<attr_adj_apply_operate<unit4> > op(
				new attr_adj_apply_operate<unit4>());
		ret = sorted_el->groupby_source(*op);
		max_col_idx = op->get_max_col_idx();
	}
	else if (el->get_attr_size() == 8) {
		std::unique_ptr<attr_adj_apply_operate<unit8> > op(
				new attr_adj_apply_operate<unit8>());
		ret = sorted_el->groupby_source(*op);
		max_col_idx = op->get_max_col_idx();
	}
	else {
		BOOST_LOG_TRIVIAL(error)
			<< "The edge attribute has an unsupported type";
		return std::pair<vector_vector::ptr, size_t>();
	}
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to groupby the edge list.\n",
			time_diff(start, end));
	return std::pair<vector_vector::ptr, size_t>(ret, max_col_idx + 1);
}

static std::pair<fg::vertex_index::ptr, detail::vec_store::ptr> create_fg_directed_graph(
		const std::string &graph_name, edge_list::ptr el)
{
	struct timeval start, end;
	// Leave the space for graph header.
	// TODO I should make this a NUMA vector.
	detail::vec_store::ptr graph_data = detail::vec_store::create(
			fg::graph_header::get_header_size(),
			get_scalar_type<char>(), -1, el->is_in_mem());
	size_t edge_data_size = el->get_attr_size();

	/*
	 * Construct the in-edge adjacency lists.
	 * All edges share the same destination vertex should be stored together.
	 */

	auto oned_mat = create_1d_matrix(el->reverse_edge());
	vector_vector::ptr in_adjs = oned_mat.first;
	size_t num_vertices = in_adjs->get_num_vecs();
	// A graph is stored in a square matrix, so the number of vertices should
	// >= the number of columns.
	assert(num_vertices >= oned_mat.second);
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

	oned_mat = create_1d_matrix(el);
	vector_vector::ptr out_adjs = oned_mat.first;
	printf("There are %ld out-edge adjacency lists and they use %ld bytes in total\n",
			out_adjs->get_num_vecs(), out_adjs->get_tot_num_entries());
	assert(num_vertices >= oned_mat.second);
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
	printf("#out edges: %d, #in edges: %ld\n",
			vector::create(num_out_edges)->sum<fg::vsize_t>(), num_edges);
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
		const std::string &graph_name, edge_list::ptr el)
{
	struct timeval start, end;
	// Leave the space for graph header.
	// TODO I should make this a NUMA vector.
	detail::vec_store::ptr graph_data = detail::vec_store::create(0,
			get_scalar_type<char>(), -1, el->is_in_mem());
	size_t edge_data_size = el->get_attr_size();

	auto oned_mat = create_1d_matrix(el);
	vector_vector::ptr adjs = oned_mat.first;
	printf("There are %ld vertices and they use %ld bytes in total\n",
			adjs->get_num_vecs(), adjs->get_tot_num_entries());

	gettimeofday(&start, NULL);
	size_t num_vertices = adjs->get_num_vecs();
	assert(num_vertices >= oned_mat.second);
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

fg::FG_graph::ptr construct_FG_graph(
		const std::pair<fg::vertex_index::ptr, detail::vec_store::ptr> &g,
		const std::string &graph_name)
{
	if (g.second->is_in_mem()) {
		fg::in_mem_graph::ptr graph = fg::in_mem_graph::create(graph_name,
				detail::smp_vec_store::cast(g.second)->get_raw_data(),
				g.second->get_length());
		return fg::FG_graph::create(graph, g.first, graph_name, NULL);
	}
	else {
		detail::EM_vec_store::ptr graph_data = detail::EM_vec_store::cast(
				g.second);
		detail::EM_vec_store::ptr index_vec = detail::EM_vec_store::create(0,
				get_scalar_type<char>());
		local_cref_vec_store index_store((const char *) g.first.get(), 0,
					g.first->get_index_size(), get_scalar_type<char>(), -1);
		index_vec->append(index_store);
		std::string graph_file_name = graph_name + ".adj";
		bool ret = graph_data->set_persistent(graph_file_name);
		if (!ret) {
			fprintf(stderr, "can't make the graph file persistent in SAFS\n");
			return fg::FG_graph::ptr();
		}
		std::string index_file_name = graph_name + ".index";
		ret = index_vec->set_persistent(index_file_name);
		if (!ret) {
			fprintf(stderr, "can't make the index file persistent in SAFS\n");
			return fg::FG_graph::ptr();
		}
		return fg::FG_graph::create(graph_file_name, index_file_name, NULL);
	}
}

fg::FG_graph::ptr create_fg_graph(const std::string &graph_name,
		edge_list::ptr el)
{
	std::pair<fg::vertex_index::ptr, detail::vec_store::ptr> res;
	if (el->is_directed())
		res = create_fg_directed_graph(graph_name, el);
	else
		res = create_fg_undirected_graph(graph_name, el);
	return construct_FG_graph(res, graph_name);
}

static vector_vector::ptr conv_fg2vv(fg::FG_graph::ptr graph, bool is_out_edge)
{
	auto vindex = graph->get_index_data();
	std::vector<off_t> offs(vindex->get_num_vertices() + 1);
	if (is_out_edge)
		init_out_offs(vindex, offs);
	else
		init_in_offs(vindex, offs);

	if (graph->is_in_mem()) {
		auto graph_data = graph->get_graph_data();
		size_t len = graph_data->get_data().get_length();
		detail::smp_vec_store::ptr mem_vec = detail::smp_vec_store::create(
				len, get_scalar_type<char>());
		graph_data->get_data().copy_to(mem_vec->get_raw_arr(), len, 0);
		return vector_vector::create(detail::mem_vv_store::create(offs, mem_vec));
	}
	else {
		safs::file_io_factory::shared_ptr factory = graph->get_graph_io_factory(
				safs::REMOTE_ACCESS);
		if (factory == NULL) {
			fprintf(stderr, "can't access the graph\n");
			return vector_vector::ptr();
		}

		detail::EM_object::file_holder::ptr holder
			= detail::EM_object::file_holder::create(factory->get_name());
		detail::EM_object::io_set::ptr ios(new detail::EM_object::io_set(factory));
		size_t file_size = factory->get_file_size();
		detail::EM_vec_store::ptr vec = detail::EM_vec_store::create(holder, ios,
				file_size, get_scalar_type<char>());
		return vector_vector::create(detail::EM_vv_store::create(offs, vec));
	}
}

namespace
{

typedef std::unordered_map<fg::vertex_id_t, fg::vertex_id_t> vertex_map_t;

class subgraph_apply: public gr_apply_operate<local_vv_store>
{
	const vertex_map_t &vmap;
	bool compact;
	// TODO I should avoid using a global variable.
	std::atomic<size_t> tot_num_edges;
public:
	subgraph_apply(const vertex_map_t &_vmap, bool compact): vmap(_vmap) {
		this->compact = compact;
		tot_num_edges = 0;
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

	size_t get_tot_num_edges() const {
		return tot_num_edges.load();
	}
};

void subgraph_apply::run(const void *key, const local_vv_store &val,
		local_vec_store &out) const
{
	factor_value_t vid = *(const factor_value_t *) key;
	const fg::ext_mem_undirected_vertex *v
		= (const fg::ext_mem_undirected_vertex *) val.get_raw_arr(0);
	assert(vid == v->get_id());

	// If we don't need the vertex, we can jump out now.
	auto it = vmap.find(vid);
	if (it == vmap.end() && compact) {
		out.resize(0);
		return;
	}
	// In this case, I need to add an empty vertex.
	else if (it == vmap.end() && !compact) {
		out.resize(fg::ext_mem_undirected_vertex::get_header_size());
		new (out.get_raw_arr()) fg::ext_mem_undirected_vertex(vid, 0, 0);
		return;
	}
	fg::vertex_id_t new_vid = it->second;

	// Get all the neighbors we need for the new vertex.
	std::vector<fg::vertex_id_t> neighs;
	for (size_t i = 0; i < v->get_num_edges(); i++) {
		auto it = vmap.find(v->get_neighbor(i));
		if (it != vmap.end())
			neighs.push_back(it->second);
	}

	// Construct the new vertex.
	size_t vsize = fg::ext_mem_undirected_vertex::num_edges2vsize(
			neighs.size(), 0);
	out.resize(vsize);
	fg::ext_mem_undirected_vertex *out_v
		= new (out.get_raw_arr()) fg::ext_mem_undirected_vertex(new_vid,
				neighs.size(), 0);
	for (size_t i = 0; i < neighs.size(); i++)
		out_v->set_neighbor(i, neighs[i]);

	const_cast<subgraph_apply *>(this)->tot_num_edges += neighs.size();
}

class set_subgraph_label_operate: public type_set_vec_operate<factor_value_t>
{
public:
	virtual void set(factor_value_t *arr, size_t num_eles,
			off_t start_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = start_idx + i;
	}
};

}

fg::FG_graph::ptr fetch_subgraph(fg::FG_graph::ptr graph,
		const std::vector<fg::vertex_id_t> &vertices,
		const std::string &graph_name, bool compact)
{
	const size_t attr_size = 0;
	// map from the original vertex id to the new id.
	std::unordered_map<fg::vertex_id_t, fg::vertex_id_t> vmap;
	if (compact) {
		for (size_t i = 0; i < vertices.size(); i++)
			vmap.insert(std::pair<fg::vertex_id_t, fg::vertex_id_t>(
						vertices[i], i));
	}
	else {
		for (size_t i = 0; i < vertices.size(); i++)
			vmap.insert(std::pair<fg::vertex_id_t, fg::vertex_id_t>(vertices[i],
						vertices[i]));
	}

	const fg::graph_header &old_header = graph->get_graph_header();
	factor f(old_header.get_num_vertices());
	factor_vector::ptr labels = factor_vector::create(f,
			old_header.get_num_vertices(), -1, true,
			set_subgraph_label_operate());

	detail::vec_store::ptr graph_data = detail::vec_store::create(
			fg::graph_header::get_header_size(),
			get_scalar_type<char>(), -1, graph->is_in_mem());
	fg::vertex_index::ptr vindex;

	size_t num_vertices = compact ? vmap.size() : old_header.get_num_vertices();
	if (old_header.is_directed_graph()) {
		vector_vector::ptr res, adjs;
		size_t num_edges;

		// Work on the in-edge lists.
		subgraph_apply in_apply(vmap, compact);
		adjs = conv_fg2vv(graph, false);
		res = adjs->groupby(*labels, in_apply);
		num_edges = in_apply.get_tot_num_edges();
		assert(res->get_num_vecs() == num_vertices);
		// Move in-edge adjacency lists to the final image.
		const detail::vv_store &in_adj_store
			= dynamic_cast<const detail::vv_store &>(res->get_data());
		graph_data->append(in_adj_store.get_data());
		// Get the number of edges of each vertex.
		std::vector<fg::vsize_t> num_in_edges(num_vertices);
		for (size_t i = 0; i < res->get_num_vecs(); i++)
			num_in_edges[i] = fg::ext_mem_undirected_vertex::vsize2num_edges(
					res->get_length(i), attr_size);

		// Work on the out-edge lists.
		subgraph_apply out_apply(vmap, compact);
		adjs = conv_fg2vv(graph, true);
		res = adjs->groupby(*labels, out_apply);
		assert(res->get_num_vecs() == num_vertices);
		assert(num_edges == out_apply.get_tot_num_edges());
		// Move out-edge adjacency lists to the final image.
		const detail::vv_store &out_adj_store
			= dynamic_cast<const detail::vv_store &>(res->get_data());
		graph_data->append(out_adj_store.get_data());
		// Get the number of edges of each vertex.
		std::vector<fg::vsize_t> num_out_edges(num_vertices);
		for (size_t i = 0; i < res->get_num_vecs(); i++)
			num_out_edges[i] = fg::ext_mem_undirected_vertex::vsize2num_edges(
					res->get_length(i), attr_size);

		// We need to free the space for these data structures.
		// They are very large for large graphs.
		vmap.clear();
		labels = NULL;

		fg::graph_header header(old_header.get_graph_type(), num_vertices,
				num_edges, attr_size);
		vindex = fg::cdirected_vertex_index::construct(num_vertices,
				num_in_edges.data(), num_out_edges.data(), header);

		// Construct the graph header.
		local_vec_store::ptr header_store(new local_buf_vec_store(0,
					fg::graph_header::get_header_size(), get_scalar_type<char>(), -1));
		memcpy(header_store->get_raw_arr(), &header,
				fg::graph_header::get_header_size());
		graph_data->set_portion(header_store, 0);
	}
	else {
		size_t num_edges;

		subgraph_apply apply(vmap, compact);
		vector_vector::ptr adjs = conv_fg2vv(graph, true);
		vector_vector::ptr res = adjs->groupby(*labels, apply);
		num_edges = apply.get_tot_num_edges();
		assert(res->get_num_vecs() == num_vertices);
		const detail::vv_store &adj_store
			= dynamic_cast<const detail::vv_store &>(res->get_data());
		graph_data->append(adj_store.get_data());
		// Get the number of edges of each vertex.
		std::vector<fg::vsize_t> num_vedges(num_vertices);
		for (size_t i = 0; i < res->get_num_vecs(); i++)
			num_vedges[i] = fg::ext_mem_undirected_vertex::vsize2num_edges(
					res->get_length(i), attr_size);

		// We need to free the space for these data structures.
		// They are very large for large graphs.
		vmap.clear();
		labels = NULL;

		fg::graph_header header(old_header.get_graph_type(), num_vertices,
				num_edges, attr_size);
		vindex = fg::cundirected_vertex_index::construct(num_vertices,
				num_vedges.data(), header);

		// Construct the graph header.
		local_vec_store::ptr header_store(new local_buf_vec_store(0,
					fg::graph_header::get_header_size(), get_scalar_type<char>(), -1));
		memcpy(header_store->get_raw_arr(), &header,
				fg::graph_header::get_header_size());
		graph_data->set_portion(header_store, 0);
	}

	return construct_FG_graph(
			std::pair<fg::vertex_index::ptr, detail::vec_store::ptr>(vindex,
				graph_data), graph_name);
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

namespace
{

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

}

class part_2d_apply_operate: public gr_apply_operate<local_vv_store>
{
	// The row length (aka. the total number of columns) of the matrix.
	size_t row_len;
	size_t nz_size;
	block_2d_size block_size;

	size_t collect_block_info(size_t block_row_id, const local_vv_store &val,
			std::vector<block_info> &infos) const;
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

static void add_nz(block_pointers &ps,
		const fg::ext_mem_undirected_vertex &v,
		const std::vector<fg::vertex_id_t> &neighs,
		block_2d_size block_size, size_t nz_size,
		char *buf)
{
	// There is only one edge from the vertex in the block
	// We store it in the COO region.
	if (neighs.size() == 1) {
		// Add index
		*ps.coos = local_coo_t(v.get_id() & block_size.get_nrow_mask(),
				neighs[0] & block_size.get_ncol_mask());
		// Add the non-zero value.
		if (nz_size > 0)
			memcpy(ps.coo_vals, v.get_raw_edge_data(0), nz_size);

		// move to the next one.
		ps.coos++;
		if (nz_size > 0)
			ps.coo_vals += nz_size;
	}
	// There are more than one edge from the vertex.
	// We store it in SCSR.
	else if (neighs.size() > 1) {
		// Add the index
		size_t row_idx = v.get_id() & block_size.get_nrow_mask();
		sparse_row_part *part = new (buf) sparse_row_part(row_idx);
		rp_edge_iterator edge_it = part->get_edge_iterator();
		for (size_t k = 0; k < neighs.size(); k++)
			edge_it.append(block_size, neighs[k]);
		ps.block->append(*part, sparse_row_part::get_size(neighs.size()));
		// Add the non-zero values.
		if (nz_size > 0)
			memcpy(ps.nz_vals, v.get_raw_edge_data(0), nz_size * neighs.size());

		// Move to the next one.
		ps.nz_vals += nz_size * neighs.size();
	}
}

size_t part_2d_apply_operate::collect_block_info(size_t block_row_id,
		const local_vv_store &val, std::vector<block_info> &block_infos) const
{
	// We calculate the number of bytes in each block accurately.
	for (size_t i = 0; i < val.get_num_vecs(); i++) {
		const fg::ext_mem_undirected_vertex *v
			= (const fg::ext_mem_undirected_vertex *) val.get_raw_arr(i);
		assert(val.get_length(i) == v->get_size());
		assert((v->get_id() >> block_size.get_nrow_log()) == block_row_id);

		// The id of the current block.
		size_t curr_bid = 0;
		// The number of edges in the block.
		size_t num_edges_block = 0;
		for (size_t j = 0; j < v->get_num_edges(); j++) {
			fg::vertex_id_t id = v->get_neighbor(j);
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

void part_2d_apply_operate::run(const void *key, const local_vv_store &val,
		local_vec_store &out) const
{
	factor_value_t block_row_id = *(const factor_value_t *) key;

	size_t num_blocks = div_ceil(row_len, block_size.get_num_cols());
	std::vector<block_info> block_infos(num_blocks);
	size_t tot_num_bytes = collect_block_info(block_row_id, val, block_infos);

	// If the block row doesn't have any non-zero entries, let's insert an
	// empty block row, so the matrix index can work correctly.
	if (tot_num_bytes == 0) {
		tot_num_bytes = sizeof(sparse_block_2d);
		out.resize(tot_num_bytes);
		sparse_block_2d *block = new (out.get_raw_arr()) sparse_block_2d(
					block_row_id, 0);
		assert(block->is_empty());
		return;
	}

	out.resize(tot_num_bytes);

	// Get the location where we can fill data to.
	std::vector<block_pointers> blocks(num_blocks);
	size_t curr_size = 0;
	for (size_t i = 0; i < block_infos.size(); i++) {
		// If the block doesn't have non-zero entries, we will skip it.
		if (block_infos[i].nnz == 0)
			continue;

		sparse_block_2d *block
			= new (out.get_raw_arr() + curr_size) sparse_block_2d(
					block_row_id, i, block_infos[i].nnz, block_infos[i].nrow,
					block_infos[i].num_coos);
		curr_size += block->get_size(nz_size);
		blocks[i].coos = block->get_coo_start();
		blocks[i].coo_vals = block->get_coo_val_start(nz_size);
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
	std::vector<fg::vertex_id_t> neighs;
	for (size_t i = 0; i < val.get_num_vecs(); i++) {
		const fg::ext_mem_undirected_vertex *v
			= (const fg::ext_mem_undirected_vertex *) val.get_raw_arr(i);

		size_t curr_bid = 0;
		neighs.clear();
		for (size_t j = 0; j < v->get_num_edges(); j++) {
			fg::vertex_id_t id = v->get_neighbor(j);
			size_t block_id = id >> block_size.get_ncol_log();
			if (curr_bid == block_id)
				neighs.push_back(id);
			else {
				add_nz(blocks[curr_bid], *v, neighs, block_size, nz_size,
						buf.get());
				curr_bid = block_id;
				neighs.clear();
				neighs.push_back(id);
			}
		}

		add_nz(blocks[curr_bid], *v, neighs, block_size, nz_size, buf.get());
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

std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr> create_2d_matrix(
		vector_vector::ptr adjs, size_t num_cols,
		const block_2d_size &block_size, const scalar_type *entry_type)
{
	size_t entry_size = 0;
	if (entry_type)
		entry_size = entry_type->get_size();
	size_t num_rows = adjs->get_num_vecs();
	factor f(ceil(((double) num_rows) / block_size.get_num_rows()));
	// TODO I should make this a NUMA vector.
	factor_vector::ptr labels = factor_vector::create(f, num_rows, -1,
			adjs->is_in_mem(), set_2d_label_operate(block_size));
	printf("groupby multiple vectors in the vector vector\n");
	struct timeval start, end;
	gettimeofday(&start, NULL);
	vector_vector::ptr res = adjs->groupby(*labels,
			part_2d_apply_operate(block_size, num_cols, entry_size));
	gettimeofday(&end, NULL);
	printf("groupby takes %f seconds\n", time_diff(start, end));

	prim_type type = prim_type::P_BOOL;
	if (entry_type)
		type = entry_type->get_type();
	matrix_header mheader(matrix_type::SPARSE, entry_size, num_rows,
			num_cols, matrix_layout_t::L_ROW_2D, type, block_size);

	// Construct the index file of the adjacency matrix.
	std::vector<off_t> offsets(res->get_num_vecs() + 1);
	off_t off = sizeof(mheader);
	for (size_t i = 0; i < res->get_num_vecs(); i++) {
		offsets[i] = off;
		off += res->get_length(i);
	}
	offsets[res->get_num_vecs()] = off;
	SpM_2d_index::ptr idx = SpM_2d_index::create(mheader, offsets);

	mheader.verify();
	return std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr>(
			idx, SpM_2d_storage::create(mheader, *res, idx));
}

std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr> create_2d_matrix(
		edge_list::ptr el, const block_2d_size &block_size,
		const scalar_type *entry_type)
{
	auto ret = create_1d_matrix(el);
	if (ret.first == NULL)
		return std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr>();
	else
		return create_2d_matrix(ret.first, ret.second, block_size, entry_type);
}

void export_2d_matrix(vector_vector::ptr adjs, size_t num_cols,
		const block_2d_size &block_size, const scalar_type *entry_type,
		const std::string &mat_file, const std::string &mat_idx_file,
		bool to_safs)
{
	size_t entry_size = 0;
	if (entry_type)
		entry_size = entry_type->get_size();
	size_t num_rows = adjs->get_num_vecs();
	factor f(ceil(((double) num_rows) / block_size.get_num_rows()));
	// TODO I should make this a NUMA vector.
	factor_vector::ptr labels = factor_vector::create(f, num_rows, -1,
			true, set_2d_label_operate(block_size));
	struct timeval start, end;
	gettimeofday(&start, NULL);

	part_2d_apply_operate op(block_size, num_cols, entry_size);

	// Allocate memory to store the groupby result.
	const scalar_type &out_type = op.get_output_type();
	// We use the storage size of the data frame to approximate the storage size
	// for the groupby result.
	size_t num_bytes = sizeof(matrix_header) + adjs->get_data().get_num_bytes();
	detail::vec_store::ptr store = detail::vec_store::create(0, out_type, -1,
			adjs->is_in_mem());
	store->reserve(num_bytes / out_type.get_size());

	// Reserve the space for the matrix header.
	local_buf_vec_store::ptr header_store(new local_buf_vec_store(0,
				sizeof(matrix_header), get_scalar_type<char>(), -1));
	store->append(*header_store);

	vector_vector::ptr res = adjs->groupby(*labels, op, store);
	gettimeofday(&end, NULL);
	printf("groupby takes %f seconds\n", time_diff(start, end));

	prim_type type = prim_type::P_BOOL;
	if (entry_type)
		type = entry_type->get_type();
	// Save the header.
	matrix_header *mheader = new(header_store->get_raw_arr()) matrix_header(
			matrix_type::SPARSE, entry_size, num_rows, num_cols,
			matrix_layout_t::L_ROW_2D, type, block_size);
	store->set_portion(header_store, 0);

	// Save the groupby result to disks.
	if (!to_safs) {
		FILE *f_2d = fopen(mat_file.c_str(), "w");
		if (f_2d == NULL) {
			BOOST_LOG_TRIVIAL(error) << boost::format("open %1%: %2%")
				% mat_file % strerror(errno);
			return;
		}
		bool ret = res->cat()->export2(f_2d);
		assert(ret);
		fclose(f_2d);
	}
	else {
		detail::EM_vec_store::ptr em_store
			= std::dynamic_pointer_cast<detail::EM_vec_store>(store);
		assert(em_store);
		em_store->set_persistent(mat_file);
	}

	// Construct the index file of the adjacency matrix.
	std::vector<off_t> offsets(res->get_num_vecs() + 1);
	off_t off = sizeof(matrix_header);
	for (size_t i = 0; i < res->get_num_vecs(); i++) {
		offsets[i] = off;
		off += res->get_length(i);
	}
	offsets[res->get_num_vecs()] = off;
	SpM_2d_index::ptr mindex = SpM_2d_index::create(*mheader, offsets);
	if (!to_safs)
		mindex->dump(mat_idx_file);
	else
		mindex->safs_dump(mat_idx_file);
}

static void print_vertex(const ext_mem_undirected_vertex &v, bool directed,
		const std::string &delim, const scalar_type *edge_data_type, FILE *f)
{
	for (size_t i = 0; i < v.get_num_edges(); i++) {
		// For undirected vertices, we only need to print the first half.
		if (!directed && v.get_neighbor(i) > v.get_id())
			break;
		std::string str = std::to_string(v.get_id()) + delim
			+ std::to_string(v.get_neighbor(i));
		if (v.has_edge_data() && edge_data_type)
			str = str + delim + edge_data_type->conv2str(
					v.get_raw_edge_data(i), 1, "");
		fprintf(f, "%s\n", str.c_str());
	}
}

// We are going to print the graph in a file. We can assume the graph is
// small enough to be stored in memory.
void print_graph_el(FG_graph::ptr fg, const std::string &delim,
		const std::string &edge_data_type, FILE *f)
{
	if (!edge_data_type.empty() && !valid_ele_type(edge_data_type)) {
		BOOST_LOG_TRIVIAL(error) << "unknown edge data type";
		return;
	}
	const scalar_type *type = NULL;
	if (!edge_data_type.empty())
		type = &get_ele_type(edge_data_type);

	// If the graph isn't stored in memory, we need to load it to memory first.
	if (!fg->is_in_mem())
		fg = FG_graph::create(fg->get_graph_data(), fg->get_index_data(),
				"tmp", fg->get_configs());

	vector_vector::ptr vv = conv_fg2vv(fg, true);
	detail::mem_vv_store::const_ptr store
		= std::dynamic_pointer_cast<const detail::mem_vv_store>(
				vv->get_raw_store());
	if (store == NULL) {
		BOOST_LOG_TRIVIAL(error) << "can't load the graph to memory";
		return;
	}
	bool directed = fg->is_directed();
	for (size_t i = 0; i < vv->get_num_vecs(); i++) {
		const ext_mem_undirected_vertex *v
			= reinterpret_cast<const ext_mem_undirected_vertex *>(
					store->get_raw_arr(i));
		print_vertex(*v, directed, delim, type, f);
	}
}

}
