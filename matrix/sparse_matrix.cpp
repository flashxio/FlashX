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

namespace fm
{

matrix_config matrix_conf;

void sparse_matrix::compute(task_creator::ptr creator) const
{
	int num_workers = matrix_conf.get_num_threads();
	int num_nodes = safs::params.get_num_nodes();
	std::vector<matrix_worker_thread::ptr> workers(num_workers);
	std::vector<matrix_io_generator::ptr> io_gens(num_workers);
	init_io_gens(io_gens);
#ifdef PROFILER
	if (!fg::graph_conf.get_prof_file().empty())
		ProfilerStart(fg::graph_conf.get_prof_file().c_str());
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
	if (!fg::graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
}

///////////// The code for sparse matrix of the FlashGraph format //////////////

void fg_row_compute_task::run(char *buf, size_t size)
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
class fg_sparse_sym_matrix: public sparse_matrix
{
	// This works like the index of the sparse matrix.
	std::vector<row_block> blocks;
	safs::file_io_factory::shared_ptr factory;

	fg_sparse_sym_matrix(safs::file_io_factory::shared_ptr factory,
			size_t nrows): sparse_matrix(nrows, true) {
		this->factory = factory;
	}
public:
	static ptr create(fg::FG_graph::ptr);

	// Nothing should happen for a symmetric matrix.
	virtual void transpose() {
	}

	virtual safs::file_io_factory::shared_ptr get_io_factory() const {
		return factory;
	}

	virtual void init_io_gens(
			std::vector<matrix_io_generator::ptr> &io_gens) const;
};

sparse_matrix::ptr fg_sparse_sym_matrix::create(fg::FG_graph::ptr fg)
{
	// Initialize vertex index.
	fg::vertex_index::ptr index = fg->get_index_data();
	assert(index != NULL);
	assert(!index->get_graph_header().is_directed_graph());

	fg::vsize_t num_vertices = index->get_num_vertices();
	fg_sparse_sym_matrix *m = new fg_sparse_sym_matrix(fg->get_graph_io_factory(
				safs::REMOTE_ACCESS), num_vertices);

	// Generate the matrix index from the vertex index.
	if (index->is_compressed()) {
		fg::in_mem_cundirected_vertex_index::ptr uindex
			= fg::in_mem_cundirected_vertex_index::create(*index);
		for (size_t i = 0; i < num_vertices;
				i += matrix_conf.get_row_block_size()) {
			fg::vertex_offset off = uindex->get_vertex(i);
			m->blocks.emplace_back(off.get_off());
		}
		size_t graph_size = uindex->get_vertex(num_vertices - 1).get_off()
			+ uindex->get_size(num_vertices - 1);
		m->blocks.emplace_back(graph_size);
	}
	else {
		fg::undirected_vertex_index::ptr uindex
			= fg::undirected_vertex_index::cast(index);
		for (size_t i = 0; i < num_vertices;
				i += matrix_conf.get_row_block_size()) {
			fg::ext_mem_vertex_info info = uindex->get_vertex_info(i);
			m->blocks.emplace_back(info.get_off());
		}
		m->blocks.emplace_back(uindex->get_graph_size());
	}

	return sparse_matrix::ptr(m);
}

void fg_sparse_sym_matrix::init_io_gens(
		std::vector<matrix_io_generator::ptr> &io_gens) const
{
	for (size_t i = 0; i < io_gens.size(); i++) {
		row_block_mapper mapper(blocks, i, io_gens.size(),
				matrix_conf.get_rb_io_size());
		io_gens[i] = matrix_io_generator::create(blocks, get_num_rows(),
				get_num_cols(), factory->get_file_id(), mapper);
	}
}

/*
 * Sparse asymmetric square matrix. It is partitioned in rows.
 */
class fg_sparse_asym_matrix: public sparse_matrix
{
	// These work like the index of the sparse matrix.
	// out_blocks index the original matrix.
	std::vector<row_block> out_blocks;
	// in_blocks index the transpose of the matrix.
	std::vector<row_block> in_blocks;
	safs::file_io_factory::shared_ptr factory;
	bool transposed;

	fg_sparse_asym_matrix(safs::file_io_factory::shared_ptr factory,
			size_t nrows): sparse_matrix(nrows, false) {
		transposed = false;
		this->factory = factory;
	}
public:
	static ptr create(fg::FG_graph::ptr);

	virtual safs::file_io_factory::shared_ptr get_io_factory() const {
		return factory;
	}

	virtual void transpose() {
		transposed = !transposed;
	}

	virtual void init_io_gens(
			std::vector<matrix_io_generator::ptr> &io_gens) const;
};

sparse_matrix::ptr fg_sparse_asym_matrix::create(fg::FG_graph::ptr fg)
{
	// Initialize vertex index.
	fg::vertex_index::ptr index = fg->get_index_data();
	assert(index != NULL);
	assert(index->get_graph_header().is_directed_graph());

	fg::vsize_t num_vertices = index->get_num_vertices();
	fg_sparse_asym_matrix *m = new fg_sparse_asym_matrix(fg->get_graph_io_factory(
				safs::REMOTE_ACCESS), num_vertices);

	if (index->is_compressed()) {
		fg::in_mem_cdirected_vertex_index::ptr dindex
			= fg::in_mem_cdirected_vertex_index::create(*index);
		for (size_t i = 0; i < num_vertices;
				i += matrix_conf.get_row_block_size()) {
			fg::directed_vertex_entry ventry = dindex->get_vertex(i);
			m->out_blocks.emplace_back(ventry.get_out_off());
			m->in_blocks.emplace_back(ventry.get_in_off());
		}
		fg::directed_vertex_entry ventry = dindex->get_vertex(num_vertices - 1);
		m->out_blocks.emplace_back(ventry.get_out_off()
				+ dindex->get_out_size(num_vertices - 1));
		m->in_blocks.emplace_back(ventry.get_in_off()
				+ dindex->get_in_size(num_vertices - 1));
	}
	else {
		fg::directed_vertex_index::ptr dindex
			= fg::directed_vertex_index::cast(index);
		// Generate the matrix index from the vertex index.
		for (size_t i = 0; i < num_vertices;
				i += matrix_conf.get_row_block_size()) {
			fg::ext_mem_vertex_info info = dindex->get_vertex_info_out(i);
			m->out_blocks.emplace_back(info.get_off());

			info = dindex->get_vertex_info_in(i);
			m->in_blocks.emplace_back(info.get_off());
		}
		fg::ext_mem_vertex_info info
			= dindex->get_vertex_info_out(num_vertices - 1);
		m->out_blocks.emplace_back(info.get_off() + info.get_size());
		info = dindex->get_vertex_info_in(num_vertices - 1);
		m->in_blocks.emplace_back(info.get_off() + info.get_size());
	}

	return sparse_matrix::ptr(m);
}

void fg_sparse_asym_matrix::init_io_gens(
		std::vector<matrix_io_generator::ptr> &io_gens) const
{
	for (size_t i = 0; i < io_gens.size(); i++) {
		row_block_mapper mapper(transposed ? in_blocks : out_blocks,
				i, io_gens.size(), matrix_conf.get_rb_io_size());
		io_gens[i] = matrix_io_generator::create(
				transposed ? in_blocks : out_blocks, get_num_rows(),
				get_num_cols(), factory->get_file_id(), mapper);
	}
}

sparse_matrix::ptr sparse_matrix::create(fg::FG_graph::ptr fg)
{
	const fg::graph_header &header = fg->get_graph_header();
	if (header.is_directed_graph())
		return fg_sparse_asym_matrix::create(fg);
	else
		return fg_sparse_sym_matrix::create(fg);
}

/////////////// The code for native 2D-partitioned sparse matrix ///////////////

class block_sparse_matrix: public sparse_matrix
{
	SpM_2d_index::ptr index;
	SpM_2d_storage::ptr mat;
	safs::file_io_factory::shared_ptr factory;
public:
	block_sparse_matrix(SpM_2d_index::ptr index,
			SpM_2d_storage::ptr mat): sparse_matrix(
				index->get_header().get_num_rows(),
				index->get_header().get_num_cols(),
				true, index->get_header().get_2d_block_size()) {
		this->index = index;
		this->mat = mat;
		factory = mat->create_io_factory();
	}

	virtual safs::file_io_factory::shared_ptr get_io_factory() const {
		return factory;
	}

	// Nothing should happen for a symmetric matrix.
	virtual void transpose() {
	}

	virtual void init_io_gens(
			std::vector<matrix_io_generator::ptr> &io_gens) const {
		for (size_t i = 0; i < io_gens.size(); i++) {
			row_block_mapper mapper(*index, i, io_gens.size(), 1);
			io_gens[i] = matrix_io_generator::create(index,
					factory->get_file_id(), mapper);
		}
	}
};

class block_sparse_asym_matrix: public sparse_matrix
{
	block_sparse_matrix::ptr mat;
	block_sparse_matrix::ptr t_mat;
public:
	block_sparse_asym_matrix(SpM_2d_index::ptr index, SpM_2d_storage::ptr mat,
			SpM_2d_index::ptr t_index, SpM_2d_storage::ptr t_mat): sparse_matrix(
				index->get_header().get_num_rows(),
				index->get_header().get_num_cols(),
				false, index->get_header().get_2d_block_size()) {
		this->mat = block_sparse_matrix::ptr(new block_sparse_matrix(index, mat));
		this->t_mat = block_sparse_matrix::ptr(new block_sparse_matrix(t_index,
					t_mat));
	}

	virtual safs::file_io_factory::shared_ptr get_io_factory() const {
		return mat->get_io_factory();
	}

	// Nothing should happen for a symmetric matrix.
	virtual void transpose() {
		block_sparse_matrix::ptr tmp = mat;
		mat = t_mat;
		t_mat = tmp;
	}

	virtual void init_io_gens(
			std::vector<matrix_io_generator::ptr> &io_gens) const {
		mat->init_io_gens(io_gens);
	}
};

sparse_matrix::ptr sparse_matrix::create(SpM_2d_index::ptr index,
		SpM_2d_storage::ptr mat)
{
	return sparse_matrix::ptr(new block_sparse_matrix(index, mat));
}

sparse_matrix::ptr sparse_matrix::create(SpM_2d_index::ptr index,
		SpM_2d_storage::ptr mat, SpM_2d_index::ptr t_index,
		SpM_2d_storage::ptr t_mat)
{
	return sparse_matrix::ptr(new block_sparse_asym_matrix(index, mat,
				t_index, t_mat));
}

static std::atomic<size_t> init_count;

void init_flash_matrix(config_map::ptr configs)
{
	size_t count = init_count.fetch_add(1);
	if (count == 0) {
		matrix_conf.init(configs);
		try {
			safs::init_io_system(configs);
		} catch (std::exception &e) {
			// If SAFS fails to initialize, we should remove the count
			// increase at the beginning of the function.
			init_count--;
			throw e;
		}
	}
}

void destroy_flash_matrix()
{
	if (init_count.fetch_sub(1) == 1)
		safs::destroy_io_system();
}

}
