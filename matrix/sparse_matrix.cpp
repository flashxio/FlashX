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
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include "sparse_matrix.h"
#include "matrix_worker_thread.h"
#include "matrix_io.h"
#include "matrix_config.h"

matrix_config matrix_conf;

void row_compute_task::run(char *buf, size_t size)
{
	assert(this->buf == buf);
	assert(this->buf_size == size);

	char *buf_p = buf + (io.get_loc().get_offset() - off);
	fg::ext_mem_undirected_vertex *v = (fg::ext_mem_undirected_vertex *) buf_p;
	for (size_t i = 0; i < io.get_num_rows(); i++) {
		size_t vsize = v->get_size();
		assert(buf_size >= vsize);
		buf_size -= vsize;
		buf_p += vsize;
		run_on_row(*v);
		v = (fg::ext_mem_undirected_vertex *) buf_p;
	}
}

/*
 * Sparse square symmetric matrix. It is partitioned in rows.
 */
class sparse_sym_matrix: public sparse_matrix
{
	// This works like the index of the sparse matrix.
	std::vector<row_block> blocks;

	sparse_sym_matrix(safs::file_io_factory::shared_ptr factory,
			size_t nrows): sparse_matrix(factory, nrows, nrows, true,
				part_dim_t::ROW) {
	}
public:
	static ptr create(fg::FG_graph::ptr);

	// Nothing should happen for a symmetric matrix.
	virtual void transpose() {
	}

	virtual void compute(task_creator::ptr creator) const;
};

sparse_matrix::ptr sparse_sym_matrix::create(fg::FG_graph::ptr fg)
{
	// Initialize vertex index.
	fg::vertex_index::ptr index = fg->get_index_data();
	assert(index != NULL);
	assert(!index->get_graph_header().is_directed_graph()
			&& !index->is_compressed());
	fg::default_vertex_index::ptr uindex = fg::default_vertex_index::cast(index);

	fg::vsize_t num_vertices = uindex->get_num_vertices();
	sparse_sym_matrix *m = new sparse_sym_matrix(fg->get_graph_io_factory(
				safs::REMOTE_ACCESS), num_vertices);

	// Generate the matrix index from the vertex index.
	for (size_t i = 0; i < num_vertices; i += matrix_conf.get_row_block_size()) {
		fg::ext_mem_vertex_info info = uindex->get_vertex_info(i);
		m->blocks.emplace_back(info.get_off());
	}
	m->blocks.emplace_back(uindex->get_graph_size());

	return sparse_matrix::ptr(m);
}

void sparse_sym_matrix::compute(task_creator::ptr creator) const
{
	int num_workers = matrix_conf.get_num_threads();
	int num_nodes = safs::params.get_num_nodes();
	std::vector<matrix_worker_thread::ptr> workers(num_workers);
	std::vector<matrix_io_generator::ptr> io_gens(num_workers);
	for (int i = 0; i < num_workers; i++) {
		io_gens[i] = matrix_io_generator::create(blocks, get_num_rows(),
				get_num_cols(), get_file_id(), i, num_workers);
	}

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif
	for (int i = 0; i < num_workers; i++) {
		int node_id = i % num_nodes;
		matrix_worker_thread::ptr t = matrix_worker_thread::create(i, node_id,
				get_io_factory(), io_gens, creator);
		t->start();
		workers[i] = t;
	}
	for (int i = 0; i < num_workers; i++)
		workers[i]->join();
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
}

/*
 * Sparse asymmetric square matrix. It is partitioned in rows.
 */
class sparse_asym_matrix: public sparse_matrix
{
	// These work like the index of the sparse matrix.
	// out_blocks index the original matrix.
	std::vector<row_block> out_blocks;
	// in_blocks index the transpose of the matrix.
	std::vector<row_block> in_blocks;
	bool transposed;

	sparse_asym_matrix(safs::file_io_factory::shared_ptr factory,
			size_t nrows): sparse_matrix(factory, nrows, nrows, false,
				part_dim_t::ROW) {
		transposed = false;
	}
public:
	static ptr create(fg::FG_graph::ptr);

	virtual void transpose() {
		transposed = !transposed;
	}

	virtual void compute(task_creator::ptr creator) const;
};

sparse_matrix::ptr sparse_asym_matrix::create(fg::FG_graph::ptr fg)
{
	// Initialize vertex index.
	fg::vertex_index::ptr index = fg->get_index_data();
	assert(index != NULL);
	assert(index->get_graph_header().is_directed_graph()
			&& !index->is_compressed());
	fg::directed_vertex_index::ptr dindex = fg::directed_vertex_index::cast(index);

	fg::vsize_t num_vertices = dindex->get_num_vertices();
	sparse_asym_matrix *m = new sparse_asym_matrix(fg->get_graph_io_factory(
				safs::REMOTE_ACCESS), num_vertices);

	// Generate the matrix index from the vertex index.
	for (size_t i = 0; i < num_vertices; i += matrix_conf.get_row_block_size()) {
		fg::ext_mem_vertex_info info = dindex->get_vertex_info_out(i);
		m->out_blocks.emplace_back(info.get_off());

		info = dindex->get_vertex_info_in(i);
		m->in_blocks.emplace_back(info.get_off());
	}
	fg::ext_mem_vertex_info info = dindex->get_vertex_info_out(num_vertices - 1);
	m->out_blocks.emplace_back(info.get_off() + info.get_size());
	info = dindex->get_vertex_info_in(num_vertices - 1);
	m->in_blocks.emplace_back(info.get_off() + info.get_size());

	return sparse_matrix::ptr(m);
}

void sparse_asym_matrix::compute(task_creator::ptr creator) const
{
	int num_workers = matrix_conf.get_num_threads();
	int num_nodes = safs::params.get_num_nodes();
	std::vector<matrix_worker_thread::ptr> workers(num_workers);
	std::vector<matrix_io_generator::ptr> io_gens(num_workers);
	for (int i = 0; i < num_workers; i++) {
		io_gens[i] = matrix_io_generator::create(
				transposed ? in_blocks : out_blocks, get_num_rows(),
				get_num_cols(), get_file_id(), i, num_workers);
	}
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif
	for (int i = 0; i < num_workers; i++) {
		int node_id = i % num_nodes;
		matrix_worker_thread::ptr t = matrix_worker_thread::create(i, node_id,
				get_io_factory(), io_gens, creator);
		t->start();
		workers[i] = t;
	}
	for (int i = 0; i < num_workers; i++)
		workers[i]->join();
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
}

sparse_matrix::ptr sparse_matrix::create(fg::FG_graph::ptr fg)
{
	const fg::graph_header &header = fg->get_graph_header();
	if (header.is_directed_graph())
		return sparse_asym_matrix::create(fg);
	else
		return sparse_sym_matrix::create(fg);
}

void init_flash_matrix(config_map::ptr configs)
{
	safs::init_io_system(configs);
	matrix_conf.init(configs);
}

void destroy_flash_matrix()
{
	safs::destroy_io_system();
}
