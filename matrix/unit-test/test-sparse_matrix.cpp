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

void print_cols(detail::mem_matrix_store::ptr store)
{
	dense_matrix::ptr mat = dense_matrix::create(store);
	if (mat->store_layout() == matrix_layout_t::L_ROW)
		mat = mat->conv2(matrix_layout_t::L_COL);
	mat->materialize_self();
	for (size_t i = 0; i < mat->get_num_cols(); i++) {
		vector::ptr col = mat->get_col(i);
		printf("%ld: %f\n", i, col->sum<float>());
	}
}

data_frame::ptr create_rand_el(bool with_attr)
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
	if (with_attr)
		df->add_vec("attr", vals);
	return df;
}

void test_spmm_block(SpM_2d_index::ptr idx, SpM_2d_storage::ptr mat,
		const std::vector<size_t> &degrees)
{
	printf("test SpMM on 2D-partitioned matrix\n");
	size_t num_cols = idx->get_header().get_num_cols();
	size_t num_rows = idx->get_header().get_num_rows();
	detail::mem_matrix_store::ptr in_mat
		= detail::NUMA_row_tall_matrix_store::create(num_cols, 10, num_nodes,
				get_scalar_type<float>());
	int val = 0;
	for (size_t i = 0; i < in_mat->get_num_rows(); i++)
		for (size_t j = 0; j < in_mat->get_num_cols(); j++)
			in_mat->set<float>(i, j, val++);
	sparse_matrix::ptr spm = sparse_matrix::create(idx, mat);

	detail::mem_matrix_store::ptr out
		= detail::NUMA_row_tall_matrix_store::create(num_rows, 10, num_nodes,
				get_scalar_type<float>());
	spm->multiply<float, float>(in_mat, out);
	print_cols(out);

	out = detail::NUMA_col_tall_matrix_store::create(num_rows, 10, num_nodes,
				get_scalar_type<float>());
	spm->multiply<float, float>(in_mat, out);
	print_cols(out);
}

void test_multiply_block(data_frame::ptr df)
{
	printf("Multiply on 2D-partitioned matrix\n");
	const block_2d_size block_size(1024, 1024);
	// Clone the original data frame.
	data_frame::ptr new_df = data_frame::create();
	for (size_t i = 0; i < df->get_num_vecs(); i++)
		new_df->add_vec(df->get_vec_name(i), df->get_vec(i)->deep_copy());
	df = new_df;

	// Get the max ID.
	vector::ptr vec = vector::create(df->get_vec(0));
	fg::vertex_id_t max_vid = vec->max<fg::vertex_id_t>();
	vec = vector::create(df->get_vec(1));
	max_vid = std::max(max_vid, vec->max<fg::vertex_id_t>());

	// I artificially add an invalid out-edge for each vertex, so it's
	// guaranteed that each vertex exists in the adjacency lists.
	detail::vec_store::ptr seq_vec = detail::create_seq_vec_store<fg::vertex_id_t>(
			0, max_vid, 1);
	detail::vec_store::ptr rep_vec = detail::create_rep_vec_store<fg::vertex_id_t>(
			max_vid + 1, fg::INVALID_VERTEX_ID);
	assert(seq_vec->get_length() == rep_vec->get_length());
	new_df = data_frame::create();
	new_df->add_vec(df->get_vec_name(0), seq_vec);
	new_df->add_vec(df->get_vec_name(1), rep_vec);
	detail::vec_store::ptr val_vec;
	size_t entry_size = 0;
	const scalar_type *entry_type = NULL;
	if (df->get_num_vecs() == 3) {
		val_vec = detail::create_rep_vec_store<float>(max_vid + 1, 0);
		new_df->add_vec(df->get_vec_name(2), val_vec);
		entry_size = val_vec->get_entry_size();
		entry_type = &val_vec->get_type();
	}
	df->append(new_df);

	auto oned_mat = create_1d_matrix(df);
	vector_vector::ptr adj = oned_mat.first;
	std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr> mat
		= create_2d_matrix(adj, oned_mat.second, block_size, entry_type);
	assert(mat.first);
	assert(mat.second);
	mat.second->verify();
	std::vector<size_t> degrees(adj->get_num_vecs());
	for (size_t i = 0; i < adj->get_num_vecs(); i++)
		degrees[i] = fg::ext_mem_undirected_vertex::vsize2num_edges(
				adj->get_length(i), entry_size);

	test_spmm_block(mat.first, mat.second, degrees);
}

void test_spmm_fg(fg::FG_graph::ptr fg)
{
	printf("test SpMM on FlashGraph matrix\n");
	sparse_matrix::ptr spm;
	if (fg->get_graph_header().has_edge_data()) {
		assert(fg->get_graph_header().get_edge_data_size() == sizeof(float));
		spm = sparse_matrix::create(fg, &get_scalar_type<float>());
	}
	else
		spm = sparse_matrix::create(fg, NULL);

	size_t num_cols = spm->get_num_cols();
	size_t num_rows = spm->get_num_rows();
	detail::mem_matrix_store::ptr in_mat
		= detail::mem_row_matrix_store::create(num_cols, 10,
				get_scalar_type<float>());
	int val = 0;
	for (size_t i = 0; i < in_mat->get_num_rows(); i++)
		for (size_t j = 0; j < in_mat->get_num_cols(); j++)
			in_mat->set<float>(i, j, val++);

	detail::mem_matrix_store::ptr out
		= detail::mem_row_matrix_store::create(num_rows, 10,
				get_scalar_type<float>());
	out->reset_data();
	spm->multiply<float, float>(in_mat, out);
	print_cols(out);
}

void test_multiply_fg(data_frame::ptr df)
{
	printf("Multiply on FlashGraph matrix\n");
	// Clone the original data frame.
	data_frame::ptr new_df = data_frame::create();
	for (size_t i = 0; i < df->get_num_vecs(); i++)
		new_df->add_vec(df->get_vec_name(i), df->get_vec(i)->deep_copy());
	df = new_df;

	// Get the max ID.
	vector::ptr vec = vector::create(df->get_vec(0));
	fg::vertex_id_t max_vid = vec->max<fg::vertex_id_t>();
	vec = vector::create(df->get_vec(1));
	max_vid = std::max(max_vid, vec->max<fg::vertex_id_t>());

	// I artificially add an invalid out-edge for each vertex, so it's
	// guaranteed that each vertex exists in the adjacency lists.
	detail::vec_store::ptr seq_vec = detail::create_seq_vec_store<fg::vertex_id_t>(
			0, max_vid, 1);
	detail::vec_store::ptr rep_vec = detail::create_rep_vec_store<fg::vertex_id_t>(
			max_vid + 1, fg::INVALID_VERTEX_ID);
	assert(seq_vec->get_length() == rep_vec->get_length());
	new_df = data_frame::create();
	new_df->add_vec(df->get_vec_name(0), seq_vec);
	new_df->add_vec(df->get_vec_name(1), rep_vec);
	detail::vec_store::ptr val_vec;
	if (df->get_num_vecs() == 3) {
		val_vec = detail::create_rep_vec_store<float>(max_vid + 1, 0);
		new_df->add_vec(df->get_vec_name(2), val_vec);
	}
	df->append(new_df);

	// I artificially add an invalid in-edge for each vertex.
	new_df = data_frame::create();
	new_df->add_vec(df->get_vec_name(1), seq_vec);
	new_df->add_vec(df->get_vec_name(0), rep_vec);
	if (df->get_num_vecs() == 3)
		new_df->add_vec(df->get_vec_name(2), val_vec);
	df->append(new_df);

	fg::FG_graph::ptr fg = create_fg_graph("test", df, true);
	test_spmm_fg(fg);
}

int main()
{
	data_frame::ptr df = create_rand_el(false);
	test_multiply_fg(df);
	test_multiply_block(df);

	df = create_rand_el(true);
	test_multiply_fg(df);
	test_multiply_block(df);
}
