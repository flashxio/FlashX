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
#include "fg_utils.h"

#include "factor.h"
#include "mem_vec_store.h"
#include "vector_vector.h"
#include "local_vv_store.h"
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

vector_vector::ptr read_csr(const std::string &row_ptr_file,
		const std::string &col_file)
{
	// Store the CSR data in the vv store.
	auto offs = read_offs(row_ptr_file);
	auto col_idxs = read_col_idxs(col_file);
	printf("%ld offs, first: %ld, last: %ld, %ld non-zero\n",
			offs.size(), offs.front(), offs.back(), col_idxs->get_length());
	// Offsets in vv store is in bytes.
	for (size_t i = 0; i < offs.size(); i++)
		offs[i] *= sizeof(col_idx_t);
	return vector_vector::create(detail::vv_store::create(offs, col_idxs));
}

int main(int argc, char *argv[])
{
	if (argc < 5) {
		fprintf(stderr,
				"csr2fg conf_file name row_ptr_file col_file [row_ptr_file col_file]\n");
		return -1;
	}

	std::string conf_file = argv[1];
	std::string adj_file = std::string(argv[2]) + ".adj";
	std::string index_file = std::string(argv[2]) + ".index";
	std::string row_ptr_file1 = argv[3];
	std::string col_file1 = argv[4];
	std::string row_ptr_file2, col_file2;
	fg::graph_type gtype = fg::graph_type::UNDIRECTED;
	if (argc == 7) {
		row_ptr_file2 = argv[5];
		col_file2 = argv[6];
		gtype = fg::graph_type::DIRECTED;
	}
	size_t edge_data_size = 0;

	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	{
		vector_vector::ptr vv = read_csr(row_ptr_file1, col_file1);

		// Count the statitistics of the graph.
		size_t num_vertices = vv->get_num_vecs();
		size_t tot_edges = 0;
		detail::smp_vec_store::ptr num_edges = detail::smp_vec_store::create(
				num_vertices, get_scalar_type<fg::vsize_t>());
		for (size_t i = 0; i < num_vertices; i++) {
			size_t local_num_edges = vv->get_length(i);
			num_edges->set<fg::vsize_t>(i, local_num_edges);
			tot_edges += local_num_edges;
		}
		// For undirected graphs.
		if (gtype == fg::graph_type::UNDIRECTED) {
			assert(tot_edges % 2 == 0);
			tot_edges /= 2;
		}
		printf("There are %ld edges\n", tot_edges);

		factor_vector::ptr labels = factor_vector::create(factor(num_vertices),
				detail::create_seq_vec_store<factor_value_t>(0, num_vertices - 1, 1));
		csr2fg_apply op;
		const scalar_type &type = op.get_output_type();

		// Construct the graph header.
		fg::graph_header header(gtype, num_vertices, tot_edges, edge_data_size);
		local_vec_store::ptr header_store(new local_buf_vec_store(0,
					fg::graph_header::get_header_size(), get_scalar_type<char>(), -1));
		memcpy(header_store->get_raw_arr(), &header,
				fg::graph_header::get_header_size());

		// Prepare the storage for the groupby result.
		size_t num_bytes = sizeof(matrix_header) + vv->get_data().get_num_bytes()
			+ sizeof(fg::ext_mem_undirected_vertex) * num_vertices;
		// The number of in-edges is the same as the number of out-edges.
		if (gtype == fg::graph_type::DIRECTED)
			num_bytes += vv->get_data().get_num_bytes();
		detail::vec_store::ptr graph_data = detail::vec_store::create(0, type, -1,
				vv->is_in_mem());
		graph_data->reserve(num_bytes / type.get_size());
		graph_data->append(*header_store);

		printf("create the graph image\n");
		vv->groupby(*labels, op, graph_data);
		vv = NULL;

		detail::smp_vec_store::ptr num_out_edges;
		if (gtype == fg::graph_type::DIRECTED) {
			printf("construct out-edges\n");
			vector_vector::ptr out_vv = read_csr(row_ptr_file2, col_file2);
			assert(out_vv->get_num_vecs() == num_vertices);
			out_vv->groupby(*labels, op, graph_data);

			num_out_edges = detail::smp_vec_store::create(num_vertices,
					get_scalar_type<fg::vsize_t>());

			size_t tot_out_edges = 0;
			for (size_t i = 0; i < num_vertices; i++) {
				size_t local_num_out_edges = out_vv->get_length(i);
				num_out_edges->set<fg::vsize_t>(i, local_num_out_edges);
				tot_out_edges += local_num_out_edges;
			}
			assert(tot_edges == tot_out_edges);
		}

		// Construct the vertex index.
		// The vectors that contains the numbers of edges have the length of #V + 1
		// because we add -1 to the edge lists artificially and the last entries
		// are the number of vertices.
		printf("create the vertex index image\n");
		fg::vertex_index::ptr vindex;
		if (gtype == fg::graph_type::DIRECTED)
			vindex = fg::cdirected_vertex_index::construct(num_vertices,
					(const fg::vsize_t *) num_edges->get_raw_arr(),
					(const fg::vsize_t *) num_out_edges->get_raw_arr(), header);
		else
			vindex = fg::cundirected_vertex_index::construct(num_vertices,
					(const fg::vsize_t *) num_edges->get_raw_arr(), header);

		fg::FG_graph::ptr graph = fg::construct_FG_graph(
				std::pair<fg::vertex_index::ptr, detail::vec_store::ptr>(vindex,
					graph_data), "");
		if (graph && graph->get_index_data())
			graph->get_index_data()->dump(index_file);
		if (graph && graph->get_graph_data())
			graph->get_graph_data()->dump(adj_file);
	}

	destroy_flash_matrix();
}
