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
#include <vector>

#include "FGlib.h"
#include "in_mem_storage.h"

#include "factor.h"
#include "mem_vec_store.h"
#include "vector_vector.h"
#include "local_vv_store.h"
#include "fm_utils.h"
#include "sparse_matrix.h"

using namespace fm;

typedef off_t row_ptr_t;
typedef int col_idx_t;

std::vector<row_ptr_t> read_offs(const std::string &file)
{
	safs::native_file f(file);
	size_t size = f.get_size();
	assert(size % sizeof(row_ptr_t) == 0);
	std::vector<row_ptr_t> offs(size / sizeof(row_ptr_t));

	FILE *fd = fopen(file.c_str(), "r");
	assert(fd);
	size_t ret = fread(offs.data(), offs.size() * sizeof(row_ptr_t), 1, fd);
	assert(ret == 1);
	return offs;
}

detail::mem_vec_store::ptr read_col_idxs(const std::string &file)
{
	safs::native_file f(file);
	size_t size = f.get_size();
	assert(size % sizeof(col_idx_t) == 0);

	const scalar_type &type = get_scalar_type<col_idx_t>();
	size_t buf_size = ROUNDUP(size, PAGE_SIZE);
	detail::mem_vec_store::ptr vec = detail::mem_vec_store::create(
			buf_size / type.get_size(), -1, type);
	char *data = vec->get_raw_arr();

	int fd = open(file.c_str(), O_RDONLY | O_DIRECT);
	assert(fd >= 0);
	size_t expected_size = size;
	while (expected_size > 0) {
		ssize_t ret = read(fd, data, buf_size);
		assert(ret >= 0);
		assert((size_t) ret <= expected_size);
		data += ret;
		buf_size -= ret;
		expected_size -= ret;
	}
	vec->resize(size / type.get_size());
	return vec;
}

class csr2fg_apply: public gr_apply_operate<local_vv_store>
{
public:
	void run(const void *key, const local_vv_store &val,
			local_vec_store &out) const {
		fg::vertex_id_t id = *(const fg::vertex_id_t *) key;
		assert(val.get_num_vecs() == 1);
		const col_idx_t *data = (const col_idx_t *) val.get_raw_arr(0);
		size_t num_edges = val.get_length(0);

		size_t size = fg::ext_mem_undirected_vertex::num_edges2vsize(num_edges, 0);
		out.resize(size);
		fg::ext_mem_undirected_vertex *v
			= new (out.get_raw_arr()) fg::ext_mem_undirected_vertex(id,
					num_edges, 0);

		for (size_t i = 0; i < num_edges; i++)
			v->set_neighbor(i, data[i]);
	}

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

int main(int argc, char *argv[])
{
	if (argc < 5) {
		fprintf(stderr,
				"csr2fg conf_file row_ptr_file col_file name\n");
		return -1;
	}

	std::string conf_file = argv[1];
	std::string row_ptr_file = argv[2];
	std::string col_file = argv[3];
	std::string adj_file = std::string(argv[4]) + ".adj";
	std::string index_file = std::string(argv[4]) + ".index";
	size_t edge_data_size = 0;

	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	{
		auto offs = read_offs(row_ptr_file);
		auto col_idxs = read_col_idxs(col_file);
		printf("%ld offs, first: %ld, last: %ld, %ld non-zero\n",
				offs.size(), offs.front(), offs.back(), col_idxs->get_length());
		// Offsets in vv store is in bytes.
		for (size_t i = 0; i < offs.size(); i++)
			offs[i] *= sizeof(col_idx_t);
		vector_vector::ptr vv = vector_vector::create(detail::vv_store::create(
					offs, col_idxs));
		size_t num_vertices = vv->get_num_vecs();
		factor_vector::ptr labels = factor_vector::create(factor(num_vertices),
				detail::create_seq_vec_store<factor_value_t>(0, num_vertices - 1, 1));
		vector_vector::ptr adjs = vv->groupby(*labels, csr2fg_apply());

		detail::smp_vec_store::ptr num_out_edges = detail::smp_vec_store::create(
				num_vertices, get_scalar_type<fg::vsize_t>());
		size_t num_edges = 0;
		for (size_t i = 0; i < num_vertices; i++) {
			size_t local_num_edges = vv->get_length(i);
			num_out_edges->set<fg::vsize_t>(i, local_num_edges);
			num_edges += local_num_edges;
		}
		assert(num_edges % 2 == 0);
		num_edges /= 2;
		printf("There are %ld edges\n", num_edges);

		printf("create the graph image\n");
		fg::graph_header header(fg::graph_type::UNDIRECTED, num_vertices, num_edges,
				edge_data_size);
		local_vec_store::ptr header_store(new local_buf_vec_store(0,
					fg::graph_header::get_header_size(), get_scalar_type<char>(), -1));
		memcpy(header_store->get_raw_arr(), &header,
				fg::graph_header::get_header_size());
		detail::vec_store::ptr graph_data = detail::vec_store::create(0,
				get_scalar_type<char>(), -1, true);
		const detail::vv_store &adj_store
			= dynamic_cast<const detail::vv_store &>(adjs->get_data());
		graph_data->reserve(header_store->get_length() + adj_store.get_num_bytes());
		graph_data->append(*header_store);
		graph_data->append(adj_store.get_data());

		// Construct the vertex index.
		// The vectors that contains the numbers of edges have the length of #V + 1
		// because we add -1 to the edge lists artificially and the last entries
		// are the number of vertices.
		printf("create the vertex index image\n");
		fg::cundirected_vertex_index::ptr vindex
			= fg::cundirected_vertex_index::construct(num_vertices,
					(const fg::vsize_t *) num_out_edges->get_raw_arr(), header);

		fg::FG_graph::ptr graph = construct_FG_graph(
				std::pair<fg::vertex_index::ptr, detail::vec_store::ptr>(vindex,
					graph_data), "");
		if (graph && graph->get_index_data())
			graph->get_index_data()->dump(index_file);
		if (graph && graph->get_graph_data())
			graph->get_graph_data()->dump(adj_file);
	}

	destroy_flash_matrix();
}
