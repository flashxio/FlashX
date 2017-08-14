#include "vertex.h"
#include "FGlib.h"

#include "sparse_matrix.h"
#include "generic_type.h"
#include "matrix_header.h"
#include "col_vec.h"
#include "EM_vector.h"
#include "EM_vv_store.h"

namespace fg
{

namespace detail
{

/*
 * This task performs computation on a sparse matrix in the FlashGraph format.
 */
class fg_row_compute_task: public fm::detail::compute_task
{
	fm::matrix_io io;
	off_t off;
	char *buf;
	size_t buf_size;
protected:
	// rows in the input matrix.
	fm::detail::row_portions::ptr in_row_portions;
	// rows in the output matrix.
	fm::detail::local_matrix_store::ptr out_rows;
	char *raw_out_rows;
public:
	fg_row_compute_task(fm::detail::matrix_store &output, const fm::matrix_io &_io,
			fm::detail::row_portions::ptr in_row_portions): io(_io) {
		off_t orig_off = io.get_loc().get_offset();
		off = ROUND_PAGE(orig_off);
		buf_size = ROUNDUP_PAGE(orig_off - off + io.get_size());
		buf = (char *) valloc(buf_size);

		this->in_row_portions = in_row_portions;
		this->out_rows = output.get_portion(_io.get_top_left().get_row_idx(),
				0, _io.get_num_rows(), output.get_num_cols());
		assert(out_rows);
		this->raw_out_rows = out_rows->get_raw_arr();
		assert(raw_out_rows);
	}

	~fg_row_compute_task() {
		free(buf);
	}
	virtual void run(char *buf, size_t size);
	virtual void run_on_row(const fg::ext_mem_undirected_vertex &v) = 0;
	virtual safs::io_request get_request() const {
		return safs::io_request(buf, safs::data_loc_t(io.get_loc().get_file_id(),
					off), buf_size, READ);
	}
};

/*
 * This task performs sparse matrix dense matrix multiplication
 * in the FlashGraph format.
 * We implement this method for the sake of compatibility. It doesn't
 * run very fast.
 */
template<class DenseType, class SparseType, int ROW_WIDTH>
class fg_row_spmm_task: public fg_row_compute_task
{
	fm::detail::matrix_store &output;
public:
	fg_row_spmm_task(fm::detail::row_portions::ptr in_row_portions,
			fm::detail::matrix_store &_output, const fm::matrix_io &_io): fg_row_compute_task(
				_output, _io, in_row_portions), output(_output) {
		assert(in_row_portions->get_type() == fm::get_scalar_type<DenseType>());
		assert(output.get_type() == fm::get_scalar_type<DenseType>());
		assert(in_row_portions->get_num_cols() == output.get_num_cols());
		assert(in_row_portions->get_num_cols() == (size_t) ROW_WIDTH);
	}

	void run_on_row(const fg::ext_mem_undirected_vertex &v);
};

template<class DenseType, class SparseType, int ROW_WIDTH>
void fg_row_spmm_task<DenseType, SparseType, ROW_WIDTH>::run_on_row(
		const fg::ext_mem_undirected_vertex &v)
{
	DenseType res[ROW_WIDTH];
	for (size_t i = 0; i < (size_t) ROW_WIDTH; i++)
		res[i] = 0;

	bool has_val = v.has_edge_data();
	for (size_t i = 0; i < v.get_num_edges(); i++) {
		fg::vertex_id_t id = v.get_neighbor(i);
		SparseType data = 1;
		if (has_val)
			data = v.get_edge_data<SparseType>(i);
		const DenseType *row = (const DenseType *) in_row_portions->get_row(id);
		for (size_t j = 0; j < (size_t) ROW_WIDTH; j++)
			res[j] += row[j] * data;
	}
	size_t rel_row_idx = v.get_id() - out_rows->get_global_start_row();
	assert(rel_row_idx < out_rows->get_num_rows());
	char *row = raw_out_rows + rel_row_idx * ROW_WIDTH * sizeof(DenseType);
	memcpy(row, res, sizeof(DenseType) * ROW_WIDTH);
}

template<class DenseType, class SparseType>
class fg_row_spmm_task<DenseType, SparseType, 0>: public fg_row_compute_task
{
	fm::detail::matrix_store &output;
public:
	fg_row_spmm_task(fm::detail::row_portions::ptr in_row_portions,
			fm::detail::matrix_store &_output, const fm::matrix_io &_io): fg_row_compute_task(
				_output, _io, in_row_portions), output(_output) {
		assert(in_row_portions->get_type() == fm::get_scalar_type<DenseType>());
		assert(output.get_type() == fm::get_scalar_type<DenseType>());
		assert(in_row_portions->get_num_cols() == output.get_num_cols());
	}

	void run_on_row(const fg::ext_mem_undirected_vertex &v) {
		DenseType res[in_row_portions->get_num_cols()];
		for (size_t i = 0; i < in_row_portions->get_num_cols(); i++)
			res[i] = 0;

		bool has_val = v.has_edge_data();
		for (size_t i = 0; i < v.get_num_edges(); i++) {
			fg::vertex_id_t id = v.get_neighbor(i);
			SparseType data = 1;
			if (has_val)
				data = v.get_edge_data<SparseType>(i);
			const DenseType *row = (const DenseType *) in_row_portions->get_row(id);
			for (size_t j = 0; j < in_row_portions->get_num_cols(); j++)
				res[j] += row[j] * data;
		}
		size_t rel_row_idx = v.get_id() - out_rows->get_global_start_row();
		assert(rel_row_idx < out_rows->get_num_rows());
		char *row = raw_out_rows
			+ rel_row_idx * output.get_num_cols() * sizeof(DenseType);
		memcpy(row, res, sizeof(DenseType) * output.get_num_cols());
	}
};

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

template<class DenseType, class SparseType>
class spmm_creator: public fm::detail::task_creator
{
	fm::detail::row_portions::ptr in_row_portions;
	fm::detail::matrix_store::ptr output;
	fm::detail::EM_matrix_stream::ptr output_stream;
	const fm::sparse_matrix &mat;

	spmm_creator(const fm::sparse_matrix &_mat, size_t num_in_cols): mat(_mat) {
	}

public:
	static fm::detail::task_creator::ptr create(const fm::sparse_matrix &mat,
			size_t num_in_cols) {
		return fm::detail::task_creator::ptr(new spmm_creator<DenseType, SparseType>(
					mat, num_in_cols));
	}

	virtual void complete() {
		if (output_stream) {
			fm::detail::EM_matrix_store::ptr em_out
				= std::dynamic_pointer_cast<fm::detail::EM_matrix_store>(output);
			assert(em_out);
			em_out->wait4complete();
			assert(output_stream->is_complete());
		}
	}

	virtual bool set_data(fm::detail::matrix_store::const_ptr in,
			fm::detail::matrix_store::ptr out, const fm::block_2d_size &block_size) {
		if (in->get_type() != fm::get_scalar_type<DenseType>()
				|| out->get_type() != fm::get_scalar_type<DenseType>()) {
			BOOST_LOG_TRIVIAL(error) << "wrong matrix type in spmm creator";
			return false;
		}
		this->output = out;
		this->in_row_portions = fm::detail::row_portions::create(in, block_size.get_num_cols());
		if (in_row_portions == NULL)
			return false;
		if (!output->is_in_mem())
			output_stream = fm::detail::EM_matrix_stream::create(
					std::dynamic_pointer_cast<fm::detail::EM_matrix_store>(output));
		return true;
	}

	virtual std::vector<fm::detail::EM_object *> get_EM_objs() {
		std::vector<fm::detail::EM_object *> ret;
		if (!output->is_in_mem()) {
			const fm::detail::EM_object *obj
				= dynamic_cast<const fm::detail::EM_object *>(output.get());
			ret.push_back(const_cast<fm::detail::EM_object *>(obj));
		}
		return ret;
	}

	virtual fm::detail::compute_task::ptr create(const fm::matrix_io &io) const {
		assert(in_row_portions);
		switch (output->get_num_cols()) {
			case 1: return fm::detail::compute_task::ptr(new fg_row_spmm_task<DenseType,
							SparseType, 1>(in_row_portions, *output, io));
			case 2: return fm::detail::compute_task::ptr(new fg_row_spmm_task<DenseType,
							SparseType, 2>(in_row_portions, *output, io));
			case 4: return fm::detail::compute_task::ptr(new fg_row_spmm_task<DenseType,
							SparseType, 4>(in_row_portions, *output, io));
			case 8: return fm::detail::compute_task::ptr(new fg_row_spmm_task<DenseType,
							SparseType, 8>(in_row_portions, *output, io));
			case 16: return fm::detail::compute_task::ptr(new fg_row_spmm_task<DenseType,
							 SparseType, 16>(in_row_portions, *output, io));
			default: return fm::detail::compute_task::ptr(new fg_row_spmm_task<DenseType,
							 SparseType, 0>(in_row_portions, *output, io));
		}
	}
};

namespace
{

class get_diag_apply: public fm::arr_apply_operate
{
	const fm::scalar_type &entry_type;
	size_t block_size;
public:
	get_diag_apply(const fm::scalar_type &_type,
			size_t block_size): entry_type(_type) {
		this->block_size = block_size;
	}

	virtual void run(const fm::local_vec_store &in,
			fm::local_vec_store &out) const;
	virtual size_t get_num_out_eles(size_t num_input) const {
		return block_size;
	}

	virtual const fm::scalar_type &get_input_type() const {
		return fm::get_scalar_type<char>();
	}
	virtual const fm::scalar_type &get_output_type() const {
		return entry_type;
	}
};

void get_diag_apply::run(const fm::local_vec_store &in,
		fm::local_vec_store &out) const
{
	out.resize(block_size);
	out.reset_data();
	auto v = reinterpret_cast<const fg::ext_mem_undirected_vertex *>(
			in.get_raw_arr());
	auto first_vid = v->get_id();
	size_t off = 0;
	size_t num_vertices = 0;
	for (size_t i = 0; i < block_size; i++) {
		num_vertices++;
		auto vid = v->get_id();
		for (size_t j = 0; j < v->get_num_edges(); j++) {
			if (v->get_neighbor(j) == vid && v->get_edge_data_size() == 0) {
				assert(out.get_type() == fm::get_scalar_type<bool>());
				out.set<bool>(vid - first_vid, true);
				break;
			}
			else if (v->get_neighbor(j) == vid) {
				assert(entry_type.get_size() == v->get_edge_data_size());
				memcpy(out.get(vid - first_vid), v->get_raw_edge_data(j),
						entry_type.get_size());
				break;
			}
			// We don't need to search any more.
			else if (v->get_neighbor(j) > vid)
				break;
		}
		off += v->get_size();
		// If we have explored all vertices in the row block, we can
		// jump out now.
		if (off >= in.get_length())
			break;
		v = reinterpret_cast<const fg::ext_mem_undirected_vertex *>(
				in.get_raw_arr() + off);
	}
	out.resize(num_vertices);
}

}

static fm::col_vec::ptr get_diag(const std::vector<fm::row_block> &blocks,
		safs::file_io_factory::shared_ptr factory,
		const fm::scalar_type &entry_type, size_t block_size)
{
	fm::detail::vec_store::ptr vec = fm::detail::vec_store::create(factory,
			fm::get_scalar_type<char>());
	std::vector<off_t> offs(blocks.size());
	for (size_t i = 0; i < blocks.size(); i++)
		offs[i] = blocks[i].get_offset();
	fm::vector_vector::ptr vv = fm::vector_vector::create(
			fm::detail::vv_store::create(offs, vec));
	fm::vector_vector::ptr ret = vv->apply(get_diag_apply(entry_type,
				block_size));
	return fm::col_vec::create(ret->cat());
}

/*
 * Sparse square symmetric matrix. It is partitioned in rows.
 */
class fg_sparse_sym_matrix: public fm::sparse_matrix
{
	fm::block_2d_size block_size;
	// This works like the index of the sparse matrix.
	std::vector<fm::row_block> blocks;
	safs::file_io_factory::shared_ptr factory;

	fg_sparse_sym_matrix(safs::file_io_factory::shared_ptr factory,
			size_t nrows, const fm::scalar_type *entry_type): fm::sparse_matrix(
				nrows, entry_type, true) {
		this->factory = factory;
	}
public:
	static ptr create(fg::FG_graph::ptr, const fm::scalar_type *entry_type);

	// Nothing should happen for a symmetric matrix.
	virtual fm::sparse_matrix::ptr transpose() const {
		return fm::sparse_matrix::ptr(new fg_sparse_sym_matrix(*this));
	}

	virtual safs::file_io_factory::shared_ptr get_io_factory() const {
		return factory;
	}

	virtual fm::matrix_io_generator::ptr create_io_gen(
			const fm::detail::matrix_store &in) const {
		assert(!in.is_wide());
		size_t num_rows = in.get_portion_size().first;
		size_t num_brows = std::min((size_t) fm::matrix_conf.get_rb_io_size(),
				num_rows / fm::matrix_conf.get_row_block_size());
		fm::row_block_mapper mapper(blocks, num_brows);
		return fm::matrix_io_generator::create(blocks, get_num_rows(),
				get_num_cols(), factory->get_file_id(), mapper);
	}

	virtual const fm::block_2d_size &get_block_size() const {
		return block_size;
	}

	virtual void get_block_row_offs(const std::vector<off_t> &block_row_idxs,
			std::vector<off_t> &offs) const {
		throw unsupported_exception(
				"get_block_row_offs isn't supported in fg_sparse_sym_matrix");
	}

	virtual fm::detail::block_exec_order::ptr get_multiply_order(
			size_t num_block_rows, size_t num_block_cols) const {
		return fm::detail::block_exec_order::ptr(new fm::detail::seq_exec_order());
	}

	virtual fm::detail::task_creator::ptr get_multiply_creator(
			const fm::scalar_type &type, size_t num_in_cols) const {
		if (type == fm::get_scalar_type<int>())
			return spmm_creator<int, int>::create(*this, num_in_cols);
		else if (type == fm::get_scalar_type<long>())
			return spmm_creator<long, long>::create(*this, num_in_cols);
		else if (type == fm::get_scalar_type<size_t>())
			return spmm_creator<size_t, size_t>::create(*this, num_in_cols);
		else if (type == fm::get_scalar_type<float>())
			return spmm_creator<float, float>::create(*this, num_in_cols);
		else if (type == fm::get_scalar_type<double>())
			return spmm_creator<double, double>::create(*this, num_in_cols);
		else {
			BOOST_LOG_TRIVIAL(error) << "unsupported type";
			return fm::detail::task_creator::ptr();
		}
	}
	virtual fm::col_vec::ptr get_diag() const {
		return fg::detail::get_diag(blocks, factory, get_type(),
				fm::matrix_conf.get_row_block_size());
	}
};

fm::sparse_matrix::ptr fg_sparse_sym_matrix::create(fg::FG_graph::ptr fg,
		const fm::scalar_type *entry_type)
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
				i += fm::matrix_conf.get_row_block_size()) {
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
				i += fm::matrix_conf.get_row_block_size()) {
			fg::ext_mem_vertex_info info = uindex->get_vertex_info(i);
			m->blocks.emplace_back(info.get_off());
		}
		m->blocks.emplace_back(uindex->get_graph_size());
	}

	return fm::sparse_matrix::ptr(m);
}

/*
 * Sparse asymmetric square matrix. It is partitioned in rows.
 */
class fg_sparse_asym_matrix: public fm::sparse_matrix
{
	fm::block_2d_size block_size;
	// These work like the index of the sparse matrix.
	// out_blocks index the original matrix.
	std::shared_ptr<std::vector<fm::row_block> > out_blocks;
	// in_blocks index the transpose of the matrix.
	std::shared_ptr<std::vector<fm::row_block> > in_blocks;
	safs::file_io_factory::shared_ptr factory;

	fg_sparse_asym_matrix(safs::file_io_factory::shared_ptr factory,
			size_t nrows, const fm::scalar_type *entry_type): fm::sparse_matrix(
				nrows, entry_type, false) {
		this->factory = factory;
		out_blocks = std::shared_ptr<std::vector<fm::row_block> >(
				new std::vector<fm::row_block>());
		in_blocks = std::shared_ptr<std::vector<fm::row_block> >(
				new std::vector<fm::row_block>());
	}
public:
	static ptr create(fg::FG_graph::ptr, const fm::scalar_type *entry_type);

	virtual safs::file_io_factory::shared_ptr get_io_factory() const {
		return factory;
	}

	virtual fm::sparse_matrix::ptr transpose() const {
		fg_sparse_asym_matrix *ret = new fg_sparse_asym_matrix(*this);
		ret->fm::sparse_matrix::_transpose();
		ret->out_blocks = this->in_blocks;
		ret->in_blocks = this->out_blocks;
		return fm::sparse_matrix::ptr(ret);
	}

	virtual fm::matrix_io_generator::ptr create_io_gen(
			const fm::detail::matrix_store &in) const {
		if (in.is_wide()) {
			BOOST_LOG_TRIVIAL(error) << boost::format(
					"The input matrix has unexpected shape: %1% x %2%")
				% in.get_num_rows() % in.get_num_cols();
			return fm::matrix_io_generator::ptr();
		}
		size_t num_rows = in.get_portion_size().first;
		size_t num_brows = std::min((size_t) fm::matrix_conf.get_rb_io_size(),
				num_rows / fm::matrix_conf.get_row_block_size());
		fm::row_block_mapper mapper(*out_blocks, num_brows);
		return fm::matrix_io_generator::create(*out_blocks, get_num_rows(),
				get_num_cols(), factory->get_file_id(), mapper);
	}

	virtual const fm::block_2d_size &get_block_size() const {
		return block_size;
	}

	virtual void get_block_row_offs(const std::vector<off_t> &block_row_idxs,
			std::vector<off_t> &offs) const {
		throw unsupported_exception(
				"get_block_row_offs isn't supported in fg_sparse_asym_matrix");
	}

	virtual fm::detail::block_exec_order::ptr get_multiply_order(
			size_t num_block_rows, size_t num_block_cols) const {
		return fm::detail::block_exec_order::ptr(new fm::detail::seq_exec_order());
	}

	virtual fm::detail::task_creator::ptr get_multiply_creator(
			const fm::scalar_type &type, size_t num_in_cols) const {
		if (type == fm::get_scalar_type<int>())
			return spmm_creator<int, int>::create(*this, num_in_cols);
		else if (type == fm::get_scalar_type<long>())
			return spmm_creator<long, long>::create(*this, num_in_cols);
		else if (type == fm::get_scalar_type<size_t>())
			return spmm_creator<size_t, size_t>::create(*this, num_in_cols);
		else if (type == fm::get_scalar_type<float>())
			return spmm_creator<float, float>::create(*this, num_in_cols);
		else if (type == fm::get_scalar_type<double>())
			return spmm_creator<double, double>::create(*this, num_in_cols);
		else {
			BOOST_LOG_TRIVIAL(error) << "unsupported type";
			return fm::detail::task_creator::ptr();
		}
	}
	virtual fm::col_vec::ptr get_diag() const {
		return fg::detail::get_diag(*in_blocks, factory, get_type(),
				fm::matrix_conf.get_row_block_size());
	}
};

fm::sparse_matrix::ptr fg_sparse_asym_matrix::create(fg::FG_graph::ptr fg,
		const fm::scalar_type *entry_type)
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
				i += fm::matrix_conf.get_row_block_size()) {
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
				i += fm::matrix_conf.get_row_block_size()) {
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

	return fm::sparse_matrix::ptr(m);
}

}

fm::sparse_matrix::ptr create_sparse_matrix(fg::FG_graph::ptr fg,
		const fm::scalar_type *entry_type)
{
	const fg::graph_header &header = fg->get_graph_header();
	if (header.is_directed_graph())
		return detail::fg_sparse_asym_matrix::create(fg, entry_type);
	else
		return detail::fg_sparse_sym_matrix::create(fg, entry_type);
}

}
