#include <unordered_set>

#include "in_mem_storage.h"

#include "fm_utils.h"
#include "sparse_matrix.h"
#include "data_frame.h"

using namespace fm;

typedef std::pair<fg::vertex_id_t, fg::vertex_id_t> edge_t;

int num_nodes = 1;

struct hash_edge
{
	size_t operator()(const edge_t &e) const {
		return e.first + e.second;
	}
};

struct edge_equal
{
	bool operator()(const edge_t &e1, const edge_t &e2) const {
		return e1.first == e2.first && e1.second == e2.second;
	}
};

data_frame::ptr create_rand_el()
{
	int num_rows = 1024 * 16;
	int num_cols = 1024 * 16;
	std::unordered_set<edge_t, hash_edge, edge_equal> edges;
	for (size_t i = 0; i < 100000; i++) {
		edge_t e;
		e.first = random() % num_rows;
		e.second = random() % num_cols;
		edges.insert(e);
	}
	printf("There are %ld edges\n", edges.size());
	detail::smp_vec_store::ptr sources = detail::smp_vec_store::create(
			edges.size(), get_scalar_type<fg::vertex_id_t>());
	detail::smp_vec_store::ptr dests = detail::smp_vec_store::create(
			edges.size(), get_scalar_type<fg::vertex_id_t>());
	detail::smp_vec_store::ptr vals = detail::smp_vec_store::create(
			edges.size(), get_scalar_type<float>());
	size_t idx = 0;
	BOOST_FOREACH(edge_t e, edges) {
		sources->set<fg::vertex_id_t>(idx, e.first);
		dests->set<fg::vertex_id_t>(idx, e.second);
		vals->set<float>(idx, idx);
		idx++;
	}

	data_frame::ptr df = data_frame::create();
	df->add_vec("source", sources);
	df->add_vec("dest", dests);
	df->add_vec("attr", vals);
	return df;
}

void test_spmv_block(SpM_2d_index::ptr idx, SpM_2d_storage::ptr mat,
		const std::vector<size_t> &degrees)
{
	printf("test SpMV on 2D-partitioned matrix\n");
	sparse_matrix::ptr spm = sparse_matrix::create(idx, mat);
	size_t num_cols = spm->get_num_cols();
	size_t num_rows = spm->get_num_rows();
	printf("test_spmv: the sparse matrix has %ld rows and %ld cols\n",
			num_rows, num_cols);
	detail::NUMA_vec_store::ptr in = detail::NUMA_vec_store::create(num_cols,
			num_nodes, get_scalar_type<int>());
	for (size_t i = 0; i < num_cols; i++)
		in->set(i, 1);
	detail::NUMA_vec_store::ptr out = detail::NUMA_vec_store::create(num_rows,
			num_nodes, get_scalar_type<int>());
	spm->multiply<int>(*in, *out);
	assert(out->get_length() == num_rows);
	assert(degrees.size() == num_rows);
	for (size_t i = 0; i < num_rows; i++)
		assert(out->get<int>(i) == degrees[i]);
}

void verify_spmm(sparse_matrix::ptr spm, detail::mem_matrix_store::ptr in_mat,
		detail::mem_matrix_store::ptr out)
{
	size_t num_cols = spm->get_num_cols();
	size_t num_rows = spm->get_num_rows();
	for (size_t i = 0; i < in_mat->get_num_cols(); i++) {
		detail::NUMA_vec_store::ptr in_vec = detail::NUMA_vec_store::create(
				num_cols, num_nodes, get_scalar_type<int>());
		for (size_t j = 0; j < num_rows; j++)
			in_vec->set<int>(j, in_mat->get<int>(j, i));
		detail::NUMA_vec_store::ptr out_vec = detail::NUMA_vec_store::create(
				num_rows, num_nodes, get_scalar_type<int>());
		spm->multiply<int>(*in_vec, *out_vec);
		for (size_t j = 0; j < num_rows; j++)
			assert(out_vec->get<int>(j) == out->get<int>(j, i));
	}
}

void test_spmm_block(SpM_2d_index::ptr idx, SpM_2d_storage::ptr mat,
		const std::vector<size_t> &degrees)
{
	printf("test SpMM on 2D-partitioned matrix\n");
	size_t num_cols = idx->get_header().get_num_cols();
	size_t num_rows = idx->get_header().get_num_rows();
	detail::mem_matrix_store::ptr in_mat
		= detail::NUMA_row_tall_matrix_store::create(num_cols, 10, num_nodes,
				get_scalar_type<int>());
	int val = 0;
	for (size_t i = 0; i < in_mat->get_num_rows(); i++)
		for (size_t j = 0; j < in_mat->get_num_cols(); j++)
			in_mat->set(i, j, val++);
	sparse_matrix::ptr spm = sparse_matrix::create(idx, mat);

	detail::mem_matrix_store::ptr out
		= detail::NUMA_row_tall_matrix_store::create(num_rows, 10, num_nodes,
				get_scalar_type<int>());
	spm->multiply<int>(*in_mat, *out);
	verify_spmm(spm, in_mat, out);

	out = detail::NUMA_col_tall_matrix_store::create(num_rows, 10, num_nodes,
				get_scalar_type<int>());
	spm->multiply<int>(*in_mat, *out);
	verify_spmm(spm, in_mat, out);
}

void test_multiply_block()
{
	printf("Multiply on 2D-partitioned matrix\n");
	const block_2d_size block_size(1024, 1024);
	data_frame::ptr df = create_rand_el();

	vector::ptr vec = vector::create(df->get_vec(0));
	fg::vertex_id_t max_vid = vec->max<fg::vertex_id_t>();
	vec = vector::create(df->get_vec(1));
	max_vid = std::max(max_vid, vec->max<fg::vertex_id_t>());

	// I artificially add an invalid out-edge for each vertex, so it's
	// guaranteed that each vertex exists in the adjacency lists.
	detail::vec_store::ptr seq_vec = detail::create_vec_store<fg::vertex_id_t>(
			0, max_vid, 1);
	detail::vec_store::ptr rep_vec = detail::create_vec_store<fg::vertex_id_t>(
			max_vid + 1, fg::INVALID_VERTEX_ID);
	detail::vec_store::ptr val_vec = detail::create_vec_store<float>(
			max_vid + 1, 0);
	assert(seq_vec->get_length() == rep_vec->get_length());
	data_frame::ptr new_df = data_frame::create();
	new_df->add_vec(df->get_vec_name(0), seq_vec);
	new_df->add_vec(df->get_vec_name(1), rep_vec);
	new_df->add_vec(df->get_vec_name(2), val_vec);
	df->append(new_df);

	vector_vector::ptr adj = create_1d_matrix(df);
	std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr> mat
		= create_2d_matrix(adj, block_size, val_vec->get_entry_size());
	assert(mat.first);
	assert(mat.second);
	mat.second->verify();
	std::vector<size_t> degrees(adj->get_num_vecs());
	for (size_t i = 0; i < adj->get_num_vecs(); i++)
		degrees[i] = fg::ext_mem_undirected_vertex::vsize2num_edges(
				adj->get_length(i), val_vec->get_entry_size());

	test_spmv_block(mat.first, mat.second, degrees);
	test_spmm_block(mat.first, mat.second, degrees);
}

void test_spmv_fg(fg::FG_graph::ptr fg)
{
	printf("test SpMV on FlashGraph matrix\n");
	sparse_matrix::ptr spm = sparse_matrix::create(fg);
	size_t num_cols = spm->get_num_cols();
	size_t num_rows = spm->get_num_rows();
	detail::smp_vec_store::ptr in = detail::smp_vec_store::create(num_cols,
			get_scalar_type<int>());
	for (size_t i = 0; i < num_cols; i++)
		in->set(i, 1);
	detail::smp_vec_store::ptr out = detail::smp_vec_store::create(num_rows,
			get_scalar_type<int>());
	spm->multiply<int>(*in, *out);
	assert(out->get_length() == num_rows);
//	assert(degrees.size() == num_rows);
//	for (size_t i = 0; i < num_rows; i++)
//		assert(out->get<int>(i) == degrees[i]);
}

void test_spmm_fg(fg::FG_graph::ptr fg)
{
	printf("test SpMM on FlashGraph matrix\n");
	sparse_matrix::ptr spm = sparse_matrix::create(fg);
	size_t num_cols = spm->get_num_cols();
	size_t num_rows = spm->get_num_rows();
	detail::mem_matrix_store::ptr in_mat
		= detail::mem_row_matrix_store::create(num_cols, 10,
				get_scalar_type<int>());
	int val = 0;
	for (size_t i = 0; i < in_mat->get_num_rows(); i++)
		for (size_t j = 0; j < in_mat->get_num_cols(); j++)
			in_mat->set(i, j, val++);

	detail::mem_matrix_store::ptr out
		= detail::mem_row_matrix_store::create(num_rows, 10,
				get_scalar_type<int>());
	spm->multiply<int>(*in_mat, *out);
//	verify_spmm(spm, in_mat, out);
}

void test_multiply_fg()
{
	printf("Multiply on FlashGraph matrix\n");
	data_frame::ptr df = create_rand_el();
	vector::ptr vec = vector::create(df->get_vec(0));
	fg::vertex_id_t max_vid = vec->max<fg::vertex_id_t>();
	vec = vector::create(df->get_vec(1));
	max_vid = std::max(max_vid, vec->max<fg::vertex_id_t>());

	// I artificially add an invalid out-edge for each vertex, so it's
	// guaranteed that each vertex exists in the adjacency lists.
	detail::vec_store::ptr seq_vec = detail::create_vec_store<fg::vertex_id_t>(
			0, max_vid, 1);
	detail::vec_store::ptr rep_vec = detail::create_vec_store<fg::vertex_id_t>(
			max_vid + 1, fg::INVALID_VERTEX_ID);
	detail::vec_store::ptr val_vec = detail::create_vec_store<float>(
			max_vid + 1, 0);
	assert(seq_vec->get_length() == rep_vec->get_length());
	data_frame::ptr new_df = data_frame::create();
	new_df->add_vec(df->get_vec_name(0), seq_vec);
	new_df->add_vec(df->get_vec_name(1), rep_vec);
	new_df->add_vec(df->get_vec_name(2), val_vec);
	df->append(new_df);

	// I artificially add an invalid in-edge for each vertex.
	new_df = data_frame::create();
	new_df->add_vec(df->get_vec_name(1), seq_vec);
	new_df->add_vec(df->get_vec_name(0), rep_vec);
	new_df->add_vec(df->get_vec_name(2), val_vec);
	df->append(new_df);

	fg::FG_graph::ptr fg = create_fg_graph("test", df, true);
	test_spmv_fg(fg);
	test_spmm_fg(fg);
}

int main()
{
	test_multiply_fg();
	test_multiply_block();
}
