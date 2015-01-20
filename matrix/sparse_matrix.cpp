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

#include "sparse_matrix.h"
#include "matrix_worker_thread.h"
#include "matrix_io.h"

static const int NUM_WORKERS = 32;
static const int NUM_NODES = 4;

void row_compute_task::run()
{
	char *buf_p = buf;
	ext_mem_undirected_vertex *v = (ext_mem_undirected_vertex *) buf_p;
	for (size_t i = 0; i < num_rows; i++) {
		size_t vsize = v->get_size();
		assert(buf_size >= vsize);
		buf_size -= vsize;
		buf_p += vsize;
		run_on_row(*v);
		v = (ext_mem_undirected_vertex *) buf_p;
	}
}

/*
 * Sparse square symmetric matrix.
 */
class sparse_sym_matrix: public sparse_matrix
{
	// This works like the index of the sparse matrix.
	std::vector<row_block> blocks;

	sparse_sym_matrix(file_io_factory::shared_ptr factory,
			size_t nrows): sparse_matrix(factory, nrows, nrows, true,
				part_dim_t::ROW) {
	}
public:
	static ptr create(FG_graph::ptr);

	virtual void compute(task_creator::ptr creator) const;
};

sparse_matrix::ptr sparse_sym_matrix::create(FG_graph::ptr fg)
{
	// Initialize vertex index.
	vertex_index::ptr index = fg->get_index_data();
	if (index == NULL)
		index = vertex_index::safs_load(fg->get_index_file());
	assert(!index->get_graph_header().is_directed_graph()
			&& !index->is_compressed());
	default_vertex_index::ptr uindex = default_vertex_index::cast(index);

	vsize_t num_vertices = uindex->get_num_vertices();
	sparse_sym_matrix *m = new sparse_sym_matrix(create_io_factory(
				fg->get_graph_file(), REMOTE_ACCESS), num_vertices);

	// Generate the matrix index from the vertex index.
	for (size_t i = 0; i < num_vertices; i += ROW_BLOCK_SIZE) {
		ext_mem_vertex_info info = uindex->get_vertex_info(i);
		m->blocks.emplace_back(info.get_off());
	}
	m->blocks.emplace_back(uindex->get_graph_size());

	return sparse_matrix::ptr(m);
}

void sparse_sym_matrix::compute(task_creator::ptr creator) const
{
	int num_workers = NUM_WORKERS;
	int num_nodes = NUM_NODES;
	std::vector<matrix_worker_thread::ptr> workers(num_workers);
	for (int i = 0; i < num_workers; i++) {
		int node_id = i % num_nodes;
		matrix_worker_thread::ptr t = matrix_worker_thread::ptr(
				matrix_worker_thread::create(node_id, get_io_factory(),
					row_io_generator::create(blocks, get_num_rows(),
						get_num_cols(), get_file_id(), i, num_workers),
					creator));
		t->start();
		workers[i] = t;
	}
	for (int i = 0; i < num_workers; i++)
		workers[i]->join();
}

sparse_matrix::ptr sparse_matrix::create(FG_graph::ptr fg)
{
	return sparse_matrix::ptr();
}
