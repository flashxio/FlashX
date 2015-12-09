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
#include "matrix_io.h"
#include "matrix_config.h"
#include "hilbert_curve.h"
#include "mem_worker_thread.h"

namespace fm
{

matrix_config matrix_conf;

/*
 * This processes the blocks in the their original order.
 * It can process an arbitrary number of blocks.
 */
class seq_exec_order: public block_exec_order
{
public:
	virtual bool is_valid_size(size_t height, size_t width) const {
		return true;
	}

	virtual bool exec(block_compute_task &task,
			std::vector<const sparse_block_2d *> &blocks) const {
		BOOST_FOREACH(const sparse_block_2d *b, blocks)
			// We allow null blocks.
			if (b)
				task.run_on_block(*b);
		return true;
	}
};

/*
 * This class processes the blocks in the hilbert order.
 * The hilbert order gives us very high CPU cache hits regardless of
 * the CPU cache size. It's kind of like cache oblivious algorithms,
 * but this is not completely cache oblivious because it still relies on
 * a right block size to generate a reasonable amount of cache hits.
 */
class hilbert_exec_order: public block_exec_order
{
	struct coordinate_order
	{
		std::pair<off_t, off_t> coo;
		size_t order;
		bool operator<(const coordinate_order &o) const {
			return this->order < o.order;
		}
	};

	size_t n;
	std::vector<coordinate_order> hilbert_orders;
public:
	hilbert_exec_order(size_t n);

	virtual bool is_valid_size(size_t height, size_t width) const {
		return n == height && n == width;
	}

	virtual bool exec(block_compute_task &task,
			std::vector<const sparse_block_2d *> &blocks) const;
};

hilbert_exec_order::hilbert_exec_order(size_t n)
{
	this->n = n;
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			coordinate_order o;
			o.coo.first = i;
			o.coo.second = j;
			o.order = hilbert_xy2d(n, i, j);
			hilbert_orders.push_back(o);
		}
	}
	assert(hilbert_orders.size() == n * n);
	// Order the coordinates ascendingly.
	std::sort(hilbert_orders.begin(), hilbert_orders.end());
}

/*
 * Process the blocks in the hilbert order.
 * Here I assume all blocks are in a square so there are n^2 blocks
 * in the vector. Furthermore, the blocks are organized in the row major
 * in the original matrix.
 */
bool hilbert_exec_order::exec(block_compute_task &task,
		std::vector<const sparse_block_2d *> &blocks) const
{
	if (blocks.size() != hilbert_orders.size()) {
		BOOST_LOG_TRIVIAL(error) << "The number of blocks need to be n^2";
		return false;
	}

	BOOST_FOREACH(coordinate_order o, hilbert_orders) {
		size_t row_idx = o.coo.first;
		size_t col_idx = o.coo.second;
		size_t idx = row_idx * n + col_idx;
		// We allow null blocks. There must be n^2 blocks in the vector.
		// If there are empty blocks in the square, we have to put a null
		// in the place to indicate the empty block.
		if (blocks[idx])
			task.run_on_block(*blocks[idx]);
	}
	return true;
}

namespace
{
struct buf_deleter {
	void operator()(char *buf) const {
		free(buf);
	}
};
}

block_compute_task::block_compute_task(const matrix_io &_io,
		const sparse_matrix &mat, block_exec_order::ptr order): io(
			_io), block_size(mat.get_block_size())
{
	this->entry_size = mat.get_entry_size();
	size_t num_block_rows
		= ceil(((double) io.get_num_rows()) / block_size.get_num_rows());
	if (order->is_valid_size(num_block_rows, num_block_rows))
		this->exec_order = order;
	else
		this->exec_order = block_exec_order::ptr(new seq_exec_order());

	off_t orig_off = io.get_loc().get_offset();
	off = ROUND_PAGE(orig_off);
	real_io_size = ROUNDUP_PAGE(orig_off - off + io.get_size());
	buf = detail::local_mem_buffer::get_irreg();
	// If there isn't a buffer available in the local thread or the local buffer
	// is smaller than required, we allocate a new buffer.
	// The smaller buffer will be deallocated if it exists.
	if (buf.second == NULL || buf.first < real_io_size) {
		std::shared_ptr<char> tmp((char *) valloc(real_io_size), buf_deleter());
		buf = detail::local_mem_buffer::irreg_buf_t(real_io_size, tmp);
	}

	// The last entry in the vector indicates the end of the last block row.
	block_rows.resize(num_block_rows + 1);
	assert(io.get_top_left().get_row_idx() % block_size.get_num_rows() == 0);
	// The first block row.
	off_t block_row_idx
		= io.get_top_left().get_row_idx() / block_size.get_num_rows();
	std::vector<off_t> block_row_idxs(num_block_rows + 1);
	for (size_t i = 0; i < block_row_idxs.size(); i++)
		block_row_idxs[i] = block_row_idx + i;
	std::vector<off_t> block_row_offs;
	mat.get_block_row_offs(block_row_idxs, block_row_offs);
	assert(io.get_loc().get_offset() == block_row_offs[0]);
	for (size_t i = 0; i < num_block_rows; i++) {
		// The offset of the block row in the buffer.
		off_t local_off = block_row_offs[i] - off;
		block_rows[i] = buf.second.get() + local_off;
	}
	block_rows[num_block_rows]
		= buf.second.get() + block_row_offs[num_block_rows] - off;
}

block_compute_task::~block_compute_task()
{
	detail::local_mem_buffer::cache_irreg(buf);
}

/*
 * A block compute task processes data in multiple block rows.
 * It's up to us in what order we should process the blocks in these block rows.
 */
void block_compute_task::run(char *buf, size_t size)
{
	off_t orig_off = io.get_loc().get_offset();
	off_t local_off = orig_off - ROUND_PAGE(orig_off);
	assert(local_off + io.get_size() <= size);
	size_t block_row_start
		= io.get_top_left().get_row_idx() / block_size.get_num_rows();
	size_t num_block_rows
		= ceil(((double) io.get_num_rows()) / block_size.get_num_rows());
	assert(io.get_top_left().get_col_idx() == 0);
	size_t block_col_start = 0;
	size_t num_block_cols
		= ceil(((double) io.get_num_cols()) / block_size.get_num_cols());

	// We access data in super blocks.
	// A super block is a set of blocks organized in a square or
	// a nearly square.
	// The last entry in the vector indicates the end of the last block row.
	size_t num_blocks = block_rows.size() - 1;
	std::vector<block_row_iterator> its;
	for (size_t i = 0; i < num_blocks; i++) {
		its.emplace_back((const sparse_block_2d *) block_rows[i],
				(const sparse_block_2d *) block_rows[i + 1]);
	}
	std::vector<const sparse_block_2d *> blocks(num_blocks * num_blocks);
	// The column index of a super block (in blocks).
	size_t sb_col_idx = 0;
	bool has_blocks;
	do {
		has_blocks = false;
		// Get a super block.
		for (size_t i = 0; i < num_blocks; i++) {
			for (size_t j = 0; j < num_blocks; j++) {
				size_t idx = i * num_blocks + j;
				// If the block row doesn't have blocks left, we should fill
				// the corresponding location with NULL.
				if (!its[i].has_next()) {
					blocks[idx] = NULL;
					continue;
				}
				const sparse_block_2d &b = its[i].get_curr();
				assert(b.get_block_row_idx() >= block_row_start
						&& b.get_block_row_idx() < block_row_start + num_block_rows);
				assert(b.get_block_col_idx() >= block_col_start
						&& b.get_block_col_idx() < block_col_start + num_block_cols);
				assert(b.get_block_col_idx() >= sb_col_idx + j);
				if (b.get_block_col_idx() == sb_col_idx + j) {
					blocks[idx] = &b;
					its[i].next(entry_size);
				}
				else
					blocks[idx] = NULL;
			}
			// As long as there is a block left in a block row, we need to
			// go through the process again.
			has_blocks |= its[i].has_next();
		}
		sb_col_idx += num_blocks;
		exec_order->exec(*this, blocks);
	} while (has_blocks);
	// We have process the entire block rows.
	notify_complete();
}

block_spmm_task::block_spmm_task(const detail::mem_matrix_store &_input,
		detail::matrix_store &_output, const matrix_io &io,
		const sparse_matrix &mat, block_exec_order::ptr order): block_compute_task(
			io, mat, order), input(_input), output(_output)
{
	entry_size = mat.get_entry_size();
	// We have to make sure the task processes the entire block rows.
	assert(io.get_num_cols() == mat.get_num_cols());
}

const char *block_spmm_task::get_in_rows(size_t start_row, size_t num_rows)
{
	const char *ret = input.get_rows(start_row, start_row + num_rows);
	assert(ret);
	return ret;
}

char *block_spmm_task::get_out_rows(size_t start_row, size_t num_rows)
{
	// out_part only needs to be initialized once because a task only
	// runs on certain block rows. 
	size_t block_row_start = get_io().get_top_left().get_row_idx();
	size_t block_num_rows = std::min(get_io().get_num_rows(),
			output.get_num_rows() - block_row_start);
	if (out_part == NULL) {
		size_t out_part_size = output.get_portion_size().first;
		size_t out_part_id = block_row_start / out_part_size;
		// It's guaranteed that all output rows are stored contiguously together.
		assert((block_row_start + block_num_rows - 1) / out_part_size
				== out_part_id);

		// We maintain a local buffer for the corresponding part of
		// the output matrix.
		out_part = detail::local_row_matrix_store::ptr(
				new detail::local_buf_row_matrix_store(block_row_start, 0,
					block_num_rows, output.get_num_cols(), output.get_type(),
					// we allocate the buffer in the local node.
					-1));
		out_part->reset_data();
	}

	// Get the contiguous rows in the input and output matrices.
	size_t local_start = start_row - out_part->get_global_start_row();
	size_t local_end = std::min(local_start + num_rows,
			out_part->get_num_rows());
	return out_part->get_rows(local_start, local_end);
}

void block_spmm_task::notify_complete()
{
	// It's possible that the entire block row is empty. In this case,
	// we didn't create out_part for the output portion. We need to reset
	// the data in the portion.
	if (out_part == NULL) {
		size_t block_row_start = get_io().get_top_left().get_row_idx();
		size_t block_num_rows = std::min(get_io().get_num_rows(),
				output.get_num_rows() - block_row_start);
		if (output.is_in_mem()) {
			detail::local_matrix_store::ptr tmp = output.get_portion(
					block_row_start, 0, block_num_rows, output.get_num_cols());
			assert(tmp);
			tmp->reset_data();
		}
		else {
			detail::local_matrix_store::ptr out_part(
					new detail::local_buf_row_matrix_store(block_row_start, 0,
						block_num_rows, output.get_num_cols(), output.get_type(),
						// we allocate the buffer in the local node.
						-1));
			out_part->reset_data();
			output.write_portion_async(out_part, block_row_start, 0);
		}
	}
	else {
		if (output.is_in_mem())
			output.get_portion(out_part->get_global_start_row(),
					out_part->get_global_start_col(), out_part->get_num_rows(),
					out_part->get_num_cols())->copy_from(*out_part);
		else
			output.write_portion_async(out_part,
					out_part->get_global_start_row(),
					out_part->get_global_start_col());
	}
}

/*
 * This is shared by all threads.
 */
class spm_dispatcher: public detail::task_dispatcher
{
	std::vector<matrix_io_generator::ptr> io_gens;
	std::vector<size_t> steal_io_ids;
	detail::EM_object::io_set::ptr ios;
	task_creator::ptr tcreator;
public:
	spm_dispatcher(const std::vector<matrix_io_generator::ptr> &io_gens,
			detail::EM_object::io_set::ptr ios, task_creator::ptr tcreator) {
		this->io_gens = io_gens;
		this->steal_io_ids.resize(io_gens.size());
		this->ios = ios;
		this->tcreator = tcreator;
	}

	virtual bool issue_task();
};

bool spm_dispatcher::issue_task()
{
	detail::pool_task_thread *curr
		= dynamic_cast<detail::pool_task_thread *>(thread::get_curr_thread());
	int thread_id = curr->get_pool_thread_id();

	matrix_io_generator::ptr this_io_gen = io_gens[thread_id];
	matrix_io mio;
	if (this_io_gen->has_next_io())
		mio = this_io_gen->get_next_io();
	else {
		size_t steal_io_id = steal_io_ids[thread_id];
		for (size_t i = 0; i < io_gens.size(); i++) {
			if (io_gens[steal_io_id]->has_next_io()) {
				mio = io_gens[steal_io_id]->steal_io();
				if (mio.is_valid())
					break;
			}
			steal_io_id = (steal_io_id + 1) % io_gens.size();
		}
		steal_io_ids[thread_id] = steal_io_id;
	}

	if (!mio.is_valid())
		return false;

	compute_task::ptr task = tcreator->create(mio);
	safs::io_request req = task->get_request();
	safs::io_interface &io = ios->get_curr_io();
	detail::portion_callback &cb = static_cast<detail::portion_callback &>(
			io.get_callback());
	cb.add(req, task);
	io.access(&req, 1);
	return true;
}

void sparse_matrix::compute(task_creator::ptr creator,
		size_t num_block_rows) const
{
	// We might have kept the memory buffers to avoid the overhead of memory
	// allocation. We should delete them all before running SpMM.
	detail::local_mem_buffer::clear_bufs();
	detail::mem_thread_pool::ptr threads
		= detail::mem_thread_pool::get_global_mem_threads();
	int num_workers = threads->get_num_threads();
	std::vector<matrix_io_generator::ptr> io_gens(num_workers);
	init_io_gens(num_block_rows, io_gens);
#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStart(matrix_conf.get_prof_file().c_str());
#endif
	if (ios == NULL) {
		sparse_matrix *mutable_this = const_cast<sparse_matrix *>(this);
		mutable_this->ios = detail::EM_object::io_set::ptr(
				new detail::EM_object::io_set(get_io_factory()));
	}
	spm_dispatcher::ptr dispatcher(new spm_dispatcher(io_gens, ios, creator));
	for (int i = 0; i < num_workers; i++) {
		detail::io_worker_task *task = new detail::io_worker_task(dispatcher);
		const detail::EM_object *sp_obj = this;
		task->register_EM_obj(const_cast<detail::EM_object *>(sp_obj));
		std::vector<detail::EM_object *> em_objs = creator->get_EM_objs();
		for (size_t i = 0; i < em_objs.size(); i++)
			task->register_EM_obj(em_objs[i]);
		threads->process_task(i % threads->get_num_nodes(), task);
	}
	threads->wait4complete();
#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	detail::local_mem_buffer::clear_bufs();
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
	block_2d_size block_size;
	// This works like the index of the sparse matrix.
	std::vector<row_block> blocks;
	safs::file_io_factory::shared_ptr factory;

	fg_sparse_sym_matrix(safs::file_io_factory::shared_ptr factory,
			size_t nrows, const scalar_type *entry_type): sparse_matrix(
				nrows, entry_type, true) {
		this->factory = factory;
	}
public:
	static ptr create(fg::FG_graph::ptr, const scalar_type *entry_type);

	// Nothing should happen for a symmetric matrix.
	virtual sparse_matrix::ptr transpose() const {
		return sparse_matrix::ptr(new fg_sparse_sym_matrix(*this));
	}

	virtual safs::file_io_factory::shared_ptr get_io_factory() const {
		return factory;
	}

	virtual void init_io_gens(size_t num_block_rows,
			std::vector<matrix_io_generator::ptr> &io_gens) const;

	virtual const block_2d_size &get_block_size() const {
		return block_size;
	}

	virtual void get_block_row_offs(const std::vector<off_t> &block_row_idxs,
			std::vector<off_t> &offs) const {
		throw unsupported_exception(
				"get_block_row_offs isn't supported in fg_sparse_sym_matrix");
	}

	virtual block_exec_order::ptr get_multiply_order(
			size_t num_block_rows, size_t num_block_cols) const {
		return block_exec_order::ptr(new seq_exec_order());
	}
};

sparse_matrix::ptr fg_sparse_sym_matrix::create(fg::FG_graph::ptr fg,
		const scalar_type *entry_type)
{
	// Initialize vertex index.
	fg::vertex_index::ptr index = fg->get_index_data();
	assert(index != NULL);
	assert(!index->get_graph_header().is_directed_graph());

	if (entry_type)
		assert(entry_type->get_size()
				== (size_t) index->get_graph_header().get_edge_data_size());
	fg::vsize_t num_vertices = index->get_num_vertices();
	fg_sparse_sym_matrix *m = new fg_sparse_sym_matrix(fg->get_graph_io_factory(
				safs::REMOTE_ACCESS), num_vertices, entry_type);

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

void fg_sparse_sym_matrix::init_io_gens(size_t num_block_rows,
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
	block_2d_size block_size;
	// These work like the index of the sparse matrix.
	// out_blocks index the original matrix.
	std::shared_ptr<std::vector<row_block> > out_blocks;
	// in_blocks index the transpose of the matrix.
	std::shared_ptr<std::vector<row_block> > in_blocks;
	safs::file_io_factory::shared_ptr factory;

	fg_sparse_asym_matrix(safs::file_io_factory::shared_ptr factory,
			size_t nrows, const scalar_type *entry_type): sparse_matrix(
				nrows, entry_type, false) {
		this->factory = factory;
		out_blocks = std::shared_ptr<std::vector<row_block> >(
				new std::vector<row_block>());
		in_blocks = std::shared_ptr<std::vector<row_block> >(
				new std::vector<row_block>());
	}
public:
	static ptr create(fg::FG_graph::ptr, const scalar_type *entry_type);

	virtual safs::file_io_factory::shared_ptr get_io_factory() const {
		return factory;
	}

	virtual sparse_matrix::ptr transpose() const {
		fg_sparse_asym_matrix *ret = new fg_sparse_asym_matrix(*this);
		ret->sparse_matrix::_transpose();
		ret->out_blocks = this->in_blocks;
		ret->in_blocks = this->out_blocks;
		return sparse_matrix::ptr(ret);
	}

	virtual void init_io_gens(size_t num_block_rows,
			std::vector<matrix_io_generator::ptr> &io_gens) const;

	virtual const block_2d_size &get_block_size() const {
		return block_size;
	}

	virtual void get_block_row_offs(const std::vector<off_t> &block_row_idxs,
			std::vector<off_t> &offs) const {
		throw unsupported_exception(
				"get_block_row_offs isn't supported in fg_sparse_asym_matrix");
	}

	virtual block_exec_order::ptr get_multiply_order(
			size_t num_block_rows, size_t num_block_cols) const {
		return block_exec_order::ptr(new seq_exec_order());
	}
};

sparse_matrix::ptr fg_sparse_asym_matrix::create(fg::FG_graph::ptr fg,
		const scalar_type *entry_type)
{
	// Initialize vertex index.
	fg::vertex_index::ptr index = fg->get_index_data();
	assert(index != NULL);
	assert(index->get_graph_header().is_directed_graph());

	if (entry_type)
		assert(entry_type->get_size()
				== (size_t) index->get_graph_header().get_edge_data_size());
	fg::vsize_t num_vertices = index->get_num_vertices();
	fg_sparse_asym_matrix *m = new fg_sparse_asym_matrix(fg->get_graph_io_factory(
				safs::REMOTE_ACCESS), num_vertices, entry_type);

	if (index->is_compressed()) {
		fg::in_mem_cdirected_vertex_index::ptr dindex
			= fg::in_mem_cdirected_vertex_index::create(*index);
		for (size_t i = 0; i < num_vertices;
				i += matrix_conf.get_row_block_size()) {
			fg::directed_vertex_entry ventry = dindex->get_vertex(i);
			m->out_blocks->emplace_back(ventry.get_out_off());
			m->in_blocks->emplace_back(ventry.get_in_off());
		}
		fg::directed_vertex_entry ventry = dindex->get_vertex(num_vertices - 1);
		m->out_blocks->emplace_back(ventry.get_out_off()
				+ dindex->get_out_size(num_vertices - 1));
		m->in_blocks->emplace_back(ventry.get_in_off()
				+ dindex->get_in_size(num_vertices - 1));
	}
	else {
		fg::directed_vertex_index::ptr dindex
			= fg::directed_vertex_index::cast(index);
		// Generate the matrix index from the vertex index.
		for (size_t i = 0; i < num_vertices;
				i += matrix_conf.get_row_block_size()) {
			fg::ext_mem_vertex_info info = dindex->get_vertex_info_out(i);
			m->out_blocks->emplace_back(info.get_off());

			info = dindex->get_vertex_info_in(i);
			m->in_blocks->emplace_back(info.get_off());
		}
		fg::ext_mem_vertex_info info
			= dindex->get_vertex_info_out(num_vertices - 1);
		m->out_blocks->emplace_back(info.get_off() + info.get_size());
		info = dindex->get_vertex_info_in(num_vertices - 1);
		m->in_blocks->emplace_back(info.get_off() + info.get_size());
	}

	return sparse_matrix::ptr(m);
}

void fg_sparse_asym_matrix::init_io_gens(size_t num_block_rows,
		std::vector<matrix_io_generator::ptr> &io_gens) const
{
	for (size_t i = 0; i < io_gens.size(); i++) {
		row_block_mapper mapper(*out_blocks, i, io_gens.size(),
				matrix_conf.get_rb_io_size());
		io_gens[i] = matrix_io_generator::create(*out_blocks, get_num_rows(),
				get_num_cols(), factory->get_file_id(), mapper);
	}
}

sparse_matrix::ptr sparse_matrix::create(fg::FG_graph::ptr fg,
		const scalar_type *entry_type)
{
	const fg::graph_header &header = fg->get_graph_header();
	if (header.is_directed_graph())
		return fg_sparse_asym_matrix::create(fg, entry_type);
	else
		return fg_sparse_sym_matrix::create(fg, entry_type);
}

/////////////// The code for native 2D-partitioned sparse matrix ///////////////

class block_sparse_matrix: public sparse_matrix
{
	// If the matrix is stored in the native format and is partitioned
	// in two dimensions, we need to know the block size.
	block_2d_size block_size;
	SpM_2d_index::ptr index;
	safs::file_io_factory::shared_ptr factory;
public:
	block_sparse_matrix(SpM_2d_index::ptr index,
			SpM_2d_storage::ptr mat): sparse_matrix(
				index->get_header().get_num_rows(),
				index->get_header().get_num_cols(),
				&index->get_header().get_data_type(), true), block_size(
				index->get_header().get_2d_block_size()) {
		this->index = index;
		factory = mat->create_io_factory();
	}

	block_sparse_matrix(SpM_2d_index::ptr index,
			safs::file_io_factory::shared_ptr factory): sparse_matrix(
				index->get_header().get_num_rows(),
				index->get_header().get_num_cols(),
				&index->get_header().get_data_type(), true), block_size(
				index->get_header().get_2d_block_size()) {
		this->index = index;
		this->factory = factory;
	}

	virtual safs::file_io_factory::shared_ptr get_io_factory() const {
		return factory;
	}

	// Nothing should happen for a symmetric matrix.
	virtual sparse_matrix::ptr transpose() const {
		return sparse_matrix::ptr(new block_sparse_matrix(*this));
	}

	virtual void init_io_gens(size_t num_block_rows,
			std::vector<matrix_io_generator::ptr> &io_gens) const {
		for (size_t i = 0; i < io_gens.size(); i++) {
			row_block_mapper mapper(*index, i, io_gens.size(), num_block_rows);
			io_gens[i] = matrix_io_generator::create(index,
					factory->get_file_id(), mapper);
		}
	}

	virtual const block_2d_size &get_block_size() const {
		return block_size;
	}

	virtual void get_block_row_offs(const std::vector<off_t> &block_row_idxs,
			std::vector<off_t> &offs) const {
		offs.resize(block_row_idxs.size());
		for (size_t i = 0; i < block_row_idxs.size(); i++)
			offs[i] = index->get_block_row_off(block_row_idxs[i]);
	}

	virtual block_exec_order::ptr get_multiply_order(
			size_t num_block_rows, size_t num_block_cols) const {
		if (num_block_rows != num_block_cols) {
			BOOST_LOG_TRIVIAL(error) << "hilbert order requires a square.";
			return block_exec_order::ptr();
		}
		double log2_nbr = log2(num_block_rows);
		if (log2_nbr != floor(log2_nbr)) {
			BOOST_LOG_TRIVIAL(error)
				<< "hilbert order requires a dimension of 2^n";
			return block_exec_order::ptr();
		}
		if (matrix_conf.use_hilbert_order())
			return block_exec_order::ptr(new hilbert_exec_order(num_block_rows));
		else
			return block_exec_order::ptr(new seq_exec_order());
	}
};

class block_sparse_asym_matrix: public sparse_matrix
{
	// If the matrix is stored in the native format and is partitioned
	// in two dimensions, we need to know the block size.
	block_2d_size block_size;
	block_sparse_matrix::ptr mat;
	block_sparse_matrix::ptr t_mat;
public:
	block_sparse_asym_matrix(SpM_2d_index::ptr index, SpM_2d_storage::ptr mat,
			SpM_2d_index::ptr t_index, SpM_2d_storage::ptr t_mat): sparse_matrix(
				index->get_header().get_num_rows(),
				index->get_header().get_num_cols(),
				&index->get_header().get_data_type(), false), block_size(
				index->get_header().get_2d_block_size()) {
		this->mat = block_sparse_matrix::ptr(new block_sparse_matrix(index, mat));
		this->t_mat = block_sparse_matrix::ptr(new block_sparse_matrix(t_index,
					t_mat));
	}

	block_sparse_asym_matrix(SpM_2d_index::ptr index,
			safs::file_io_factory::shared_ptr mat_io_fac,
			SpM_2d_index::ptr t_index,
			safs::file_io_factory::shared_ptr t_mat_io_fac): sparse_matrix(
				index->get_header().get_num_rows(),
				index->get_header().get_num_cols(),
				&index->get_header().get_data_type(), false), block_size(
				index->get_header().get_2d_block_size()) {
		this->mat = block_sparse_matrix::ptr(new block_sparse_matrix(index, mat_io_fac));
		this->t_mat = block_sparse_matrix::ptr(new block_sparse_matrix(t_index,
					t_mat_io_fac));
	}

	virtual safs::file_io_factory::shared_ptr get_io_factory() const {
		return mat->get_io_factory();
	}

	// Nothing should happen for a symmetric matrix.
	virtual sparse_matrix::ptr transpose() const {
		block_sparse_asym_matrix *ret = new block_sparse_asym_matrix(*this);
		ret->mat = this->t_mat;
		ret->t_mat = this->mat;
		ret->sparse_matrix::_transpose();
		return sparse_matrix::ptr(ret);
	}

	virtual void init_io_gens(size_t num_block_rows,
			std::vector<matrix_io_generator::ptr> &io_gens) const {
		mat->init_io_gens(num_block_rows, io_gens);
	}

	virtual const block_2d_size &get_block_size() const {
		return block_size;
	}

	virtual void get_block_row_offs(const std::vector<off_t> &block_row_idxs,
			std::vector<off_t> &offs) const {
		mat->get_block_row_offs(block_row_idxs, offs);
	}

	virtual block_exec_order::ptr get_multiply_order(
			size_t num_block_rows, size_t num_block_cols) const {
		return mat->get_multiply_order(num_block_rows, num_block_cols);
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

sparse_matrix::ptr sparse_matrix::create(SpM_2d_index::ptr index,
			safs::file_io_factory::shared_ptr mat_io_fac)
{
	return sparse_matrix::ptr(new block_sparse_matrix(index, mat_io_fac));
}

sparse_matrix::ptr sparse_matrix::create(SpM_2d_index::ptr index,
			safs::file_io_factory::shared_ptr mat_io_fac,
			SpM_2d_index::ptr t_index,
			safs::file_io_factory::shared_ptr t_mat_io_fac)
{
	return sparse_matrix::ptr(new block_sparse_asym_matrix(index, mat_io_fac,
				t_index, t_mat_io_fac));
}

bool sparse_matrix::multiply(detail::matrix_store::const_ptr in,
		detail::matrix_store::ptr out, task_creator::ptr create) const
{
	if (in->get_num_rows() != ncols
			|| in->get_num_cols() != out->get_num_cols()
			|| out->get_num_rows() != this->get_num_rows()) {
		BOOST_LOG_TRIVIAL(error) <<
			"the input and output matrix have incompatible dimensions";
		return false;
	}

	dense_matrix::ptr in_tmp;
	if (in->store_layout() != matrix_layout_t::L_ROW) {
		in_tmp = dense_matrix::create(in);
		in_tmp = in_tmp->conv2(matrix_layout_t::L_ROW);
	}
	if (!in->is_in_mem()) {
		if (in_tmp == NULL)
			in_tmp = dense_matrix::create(in);
		in_tmp = in_tmp->conv_store(true, matrix_conf.get_num_nodes());
	}

	size_t sblock_size = cal_super_block_size(get_block_size(),
			in->get_entry_size() * in->get_num_cols());
	detail::mem_matrix_store::const_ptr mem_in;
	if (in_tmp) {
		in_tmp->materialize_self();
		mem_in = detail::mem_matrix_store::cast(in_tmp->get_raw_store());
	}
	else
		mem_in = detail::mem_matrix_store::cast(in);
	bool ret = create->set_data(mem_in, out);
	if (ret)
		compute(create, sblock_size);
	return ret;
}

std::vector<safs::io_interface::ptr> sparse_matrix::create_ios() const
{
	std::vector<safs::io_interface::ptr> ret(1);
	assert(ios);
	ret[0] = ios->create_io();
	return ret;
}

static std::atomic<size_t> init_count;

void init_flash_matrix(config_map::ptr configs)
{
	size_t count = init_count.fetch_add(1);
	if (count == 0) {
		if (configs) {
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
		size_t num_nodes = matrix_conf.get_num_nodes();
		size_t num_threads = matrix_conf.get_num_DM_threads();
		detail::local_mem_buffer::init();
		detail::mem_thread_pool::init_global_mem_threads(num_nodes,
				num_threads / num_nodes);
	}
}

void destroy_flash_matrix()
{
	if (init_count.fetch_sub(1) == 1) {
		safs::destroy_io_system();
		detail::local_mem_buffer::destroy();
		// TODO I should also destroy the worker thread here.
	}
}

}
