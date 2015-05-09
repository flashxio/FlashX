#include <unordered_set>

#include "fm_utils.h"
#include "sparse_matrix.h"
#include "mem_data_frame.h"

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
	fg::vertex_id_t max_vid = 0;
	for (size_t i = 0; i < 100000; i++) {
		edge_t e;
		e.first = random() % num_rows;
		e.second = random() % num_cols;
		max_vid = std::max(max_vid, e.first);
		max_vid = std::max(max_vid, e.second);
		edges.insert(e);
	}
	printf("There are %ld edges\n", edges.size());
	mem_vector::ptr sources = mem_vector::create(edges.size(),
			get_scalar_type<fg::vertex_id_t>());
	mem_vector::ptr dests = mem_vector::create(edges.size(),
			get_scalar_type<fg::vertex_id_t>());
	size_t idx = 0;
	BOOST_FOREACH(edge_t e, edges) {
		sources->set(idx, e.first);
		dests->set(idx, e.second);
		idx++;
	}
	mem_data_frame::ptr df = mem_data_frame::create();
	df->add_vec("source", sources);
	df->add_vec("dest", dests);

	vector::ptr seq_vec = create_vector<fg::vertex_id_t>(0, max_vid, 1);
	vector::ptr rep_vec = create_vector<fg::vertex_id_t>(max_vid + 1,
			fg::INVALID_VERTEX_ID);
	assert(seq_vec->get_length() == rep_vec->get_length());
	// I artificially add an invalid out-edge for each vertex, so it's
	// guaranteed that each vertex exists in the adjacency lists.
	mem_data_frame::ptr new_df = mem_data_frame::create();
	new_df->add_vec("source", seq_vec);
	new_df->add_vec("dest", rep_vec);
	df->append(new_df);
	return df;
}

void test_spmv(SpM_2d_index::ptr idx, SpM_2d_storage::ptr mat,
		const std::vector<size_t> &degrees)
{
	sparse_matrix::ptr spm = sparse_matrix::create(idx, mat);
	size_t num_cols = spm->get_num_cols();
	size_t num_rows = spm->get_num_rows();
	printf("test_spmv: the sparse matrix has %ld rows and %ld cols\n",
			num_rows, num_cols);
	NUMA_vector::ptr in = NUMA_vector::create(num_cols, get_scalar_type<int>());
	for (size_t i = 0; i < num_cols; i++)
		in->set(i, 1);
	NUMA_vector::ptr out = NUMA_vector::create(num_rows, get_scalar_type<int>());
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
		NUMA_vector::ptr in_vec = NUMA_vector::create(num_cols, num_nodes,
				get_scalar_type<int>());
		for (size_t j = 0; j < num_rows; j++)
			in_vec->set<int>(j, in_mat->get<int>(j, i));
		NUMA_vector::ptr out_vec = NUMA_vector::create(num_rows, num_nodes,
				get_scalar_type<int>());
		spm->multiply<int>(*in_vec, *out_vec);
		for (size_t j = 0; j < num_rows; j++)
			assert(out_vec->get<int>(j) == out->get<int>(j, i));
	}
}

void test_spmm(SpM_2d_index::ptr idx, SpM_2d_storage::ptr mat,
		const std::vector<size_t> &degrees)
{
	size_t num_cols = idx->get_header().get_num_cols();
	size_t num_rows = idx->get_header().get_num_rows();
	printf("test_spmm: the sparse matrix has %ld rows and %ld cols\n",
			num_rows, num_cols);
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

int main()
{

	const block_2d_size block_size(1024, 1024);
	data_frame::ptr df = create_rand_el();
	vector_vector::ptr adj = create_1d_matrix(df);
	std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr> mat
		= create_2d_matrix(adj, block_size);
	assert(mat.first);
	assert(mat.second);
	mat.second->verify();
	std::vector<size_t> degrees(adj->get_num_vecs());
	for (size_t i = 0; i < adj->get_num_vecs(); i++)
		degrees[i] = fg::ext_mem_undirected_vertex::vsize2num_edges(
				adj->get_length(i), 0);

	test_spmv(mat.first, mat.second, degrees);
	test_spmm(mat.first, mat.second, degrees);
}
