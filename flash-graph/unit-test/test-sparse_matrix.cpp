#include <unordered_set>

#include "in_mem_storage.h"

#include "fg_utils.h"
#include "sparse_matrix.h"
#include "data_frame.h"
#include "col_vec.h"

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
	mat = mat->cast_ele_type(get_scalar_type<double>());
	dense_matrix::ptr sum = mat->col_sum();
	sum->materialize_self();
	detail::mem_matrix_store::const_ptr mem_sum
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				sum->get_raw_store());
	for (size_t i = 0; i < mem_sum->get_num_cols(); i++)
		printf("%ld: %f\n", i, mem_sum->get<double>(0, i));
}

fg::edge_list::ptr create_rand_el(bool with_attr, bool directed)
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
	// For an undirected graph, we need to store an edge twice.
	size_t num_edges = directed ? edges.size() : edges.size() * 2;
	detail::smp_vec_store::ptr sources = detail::smp_vec_store::create(
			num_edges, get_scalar_type<fg::vertex_id_t>());
	detail::smp_vec_store::ptr dests = detail::smp_vec_store::create(
			num_edges, get_scalar_type<fg::vertex_id_t>());
	detail::smp_vec_store::ptr vals = detail::smp_vec_store::create(
			num_edges, get_scalar_type<float>());
	size_t idx = 0;
	BOOST_FOREACH(edge_t e, edges) {
		if (directed) {
			sources->set<fg::vertex_id_t>(idx, e.first);
			dests->set<fg::vertex_id_t>(idx, e.second);
			vals->set<float>(idx, idx);
		}
		else {
			sources->set<fg::vertex_id_t>(idx * 2, e.first);
			dests->set<fg::vertex_id_t>(idx * 2, e.second);
			vals->set<float>(idx * 2, idx);

			sources->set<fg::vertex_id_t>(idx * 2 + 1, e.second);
			dests->set<fg::vertex_id_t>(idx * 2 + 1, e.first);
			vals->set<float>(idx * 2 + 1, idx);
		}
		idx++;
	}

	data_frame::ptr df = data_frame::create();
	df->add_vec("source", sources);
	df->add_vec("dest", dests);
	if (with_attr)
		df->add_vec("attr", vals);
	return fg::edge_list::create(df, directed);
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

	detail::mem_matrix_store::ptr out1
		= detail::NUMA_row_tall_matrix_store::create(num_rows, 10, num_nodes,
				get_scalar_type<float>());
	spm->multiply(in_mat, out1);
	print_cols(out1);

	detail::mem_matrix_store::ptr out2
		= detail::NUMA_col_tall_matrix_store::create(num_rows, 10, num_nodes,
				get_scalar_type<float>());
	spm->multiply(in_mat, out2);
	print_cols(out2);

	dense_matrix::ptr m1 = dense_matrix::create(out1);
	dense_matrix::ptr m2 = dense_matrix::create(out2);
	scalar_variable::ptr diff = m1->minus(*m2)->abs()->sum();
	printf("diff: %f\n", scalar_variable::get_val<float>(*diff));
}

void test_multiply_block(fg::edge_list::ptr el)
{
	printf("Multiply on 2D-partitioned matrix\n");
	const block_2d_size block_size(1024, 1024);

	auto oned_mat = fg::create_1d_matrix(el);
	vector_vector::ptr adj = oned_mat.first;
	size_t entry_size = el->get_attr_size();
	const scalar_type *entry_type = NULL;
	if (el->has_attr())
		entry_type = &el->get_attr_type();
	std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr> mat
		= fg::create_2d_matrix(adj, oned_mat.second, block_size, entry_type);
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
		spm = create_sparse_matrix(fg, &get_scalar_type<float>());
	}
	else
		spm = create_sparse_matrix(fg, NULL);

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
	spm->multiply(in_mat, out);
	print_cols(out);
}

void test_multiply_fg(fg::edge_list::ptr el)
{
	printf("Multiply on FlashGraph matrix\n");
	fg::FG_graph::ptr fg = fg::create_fg_graph("test", el);
	test_spmm_fg(fg);
}

fm::col_vec::ptr get_diag_direct(fg::edge_list::ptr el, size_t num_vertices)
{
	auto src = std::dynamic_pointer_cast<const fm::detail::mem_vec_store>(
			el->get_source());
	auto dst = std::dynamic_pointer_cast<const fm::detail::mem_vec_store>(
			el->get_dest());
	fm::detail::mem_vec_store::const_ptr attr;
	if (el->has_attr())
		attr = std::dynamic_pointer_cast<const fm::detail::mem_vec_store>(
				el->get_attr());
	fm::detail::mem_matrix_store::ptr res = fm::detail::mem_matrix_store::create(
			num_vertices, 1, fm::matrix_layout_t::L_COL,
			attr ? attr->get_type() : fm::get_scalar_type<bool>(), -1);
	res->reset_data();
	for (size_t i = 0; i < src->get_length(); i++) {
		fg::vertex_id_t id = src->get<fg::vertex_id_t>(i);
		if (id == fg::INVALID_VERTEX_ID)
			continue;

		assert(id < res->get_num_rows());
		if (id == dst->get<fg::vertex_id_t>(i)) {
			if (attr)
				memcpy(res->get(id, 0), attr->get(i), attr->get_type().get_size());
			else
				res->set<bool>(id, 0, 1);
		}
	}
	return fm::col_vec::create(res);
}

void test_diag(fg::edge_list::ptr el)
{
	fg::FG_graph::ptr fg = fg::create_fg_graph("test", el);
	fm::sparse_matrix::ptr spm;
	if (el->has_attr())
		spm = create_sparse_matrix(fg, &el->get_attr_type());
	else
		spm = create_sparse_matrix(fg, NULL);
	fm::dense_matrix::ptr vec = spm->get_diag();
	assert(vec->get_type() == spm->get_type());
	assert(vec->get_num_rows() == spm->get_num_rows());
	fm::dense_matrix::ptr true_diag = get_diag_direct(el, spm->get_num_rows());
	if (vec->get_type() == fm::get_scalar_type<bool>())
		vec = vec->cast_ele_type(fm::get_scalar_type<int>());
	if (true_diag->get_type() == fm::get_scalar_type<bool>())
		true_diag = true_diag->cast_ele_type(fm::get_scalar_type<int>());
	auto comp_sum = vec->minus(*true_diag)->abs()->sum();
	if (comp_sum->get_type() == fm::get_scalar_type<float>())
		assert(fm::scalar_variable::get_val<float>(*comp_sum) == 0);
	else
		assert(fm::scalar_variable::get_val<int>(*comp_sum) == 0);
}

void test_diag()
{
	fg::set_remove_self_edge(false);
	printf("test on a directed graph without attributes\n");
	fg::edge_list::ptr el = create_rand_el(false, true);
	test_diag(el);
	printf("test on an undirected graph without attributes\n");
	el = create_rand_el(false, false);
	test_diag(el);
	printf("test on a directed graph with attributes\n");
	el = create_rand_el(true, true);
	test_diag(el);
	printf("test on an undirected graph with attributes\n");
	el = create_rand_el(true, false);
	test_diag(el);
	fg::set_remove_self_edge(true);
}

int main()
{
	init_flash_matrix(NULL);
	test_diag();
	fg::edge_list::ptr el = create_rand_el(false, true);
	test_multiply_fg(el);
	test_multiply_block(el);

	el = create_rand_el(true, true);
	test_multiply_fg(el);
	test_multiply_block(el);
}
