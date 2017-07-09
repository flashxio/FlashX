#include <unordered_set>

#include "fm_utils.h"
#include "sparse_matrix.h"
#include "data_frame.h"

using namespace fm;

typedef std::pair<uint32_t, uint32_t> edge_t;

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
			edges.size(), get_scalar_type<uint32_t>());
	detail::smp_vec_store::ptr dests = detail::smp_vec_store::create(
			edges.size(), get_scalar_type<uint32_t>());
	detail::smp_vec_store::ptr vals = detail::smp_vec_store::create(
			edges.size(), get_scalar_type<float>());
	size_t idx = 0;
	BOOST_FOREACH(edge_t e, edges) {
		sources->set<uint32_t>(idx, e.first);
		dests->set<uint32_t>(idx, e.second);
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

void test_spmm_block(sparse_matrix::ptr spm)
{
	printf("test SpMM on 2D-partitioned matrix\n");
	size_t num_cols = spm->get_num_cols();
	size_t num_rows = spm->get_num_rows();
	detail::mem_matrix_store::ptr in_mat
		= detail::NUMA_row_tall_matrix_store::create(num_cols, 10, num_nodes,
				get_scalar_type<float>());
	int val = 0;
	for (size_t i = 0; i < in_mat->get_num_rows(); i++)
		for (size_t j = 0; j < in_mat->get_num_cols(); j++)
			in_mat->set<float>(i, j, val++);

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

	dense_matrix::ptr EM_in = dense_matrix::create(in_mat);
	EM_in = EM_in->conv_store(false, -1);
	dense_matrix::ptr EM_out1 = spm->multiply(EM_in);
	dense_matrix::ptr EM_out2 = spm->multiply(EM_in,
			EM_in->get_num_rows() * 3 * EM_in->get_entry_size());

	diff = EM_out1->minus(*m1)->abs()->sum();
	assert(scalar_variable::get_val<float>(*diff) == 0);
	diff = EM_out2->minus(*m1)->abs()->sum();
	assert(scalar_variable::get_val<float>(*diff) == 0);
}

void test_multiply_block(data_frame::ptr el)
{
	printf("Multiply on 2D-partitioned matrix\n");
	const block_2d_size block_size(1024, 1024);

	const scalar_type *entry_type = NULL;
	if (el->get_num_vecs() >= 3)
		entry_type = &el->get_vec(2)->get_type();
	sparse_matrix::ptr spm = create_2d_matrix(el, block_size,
			entry_type, false);

	test_spmm_block(spm);
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "test conf_file\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);
	int num_nodes = matrix_conf.get_num_nodes();

	data_frame::ptr el = create_rand_el(false);
	test_multiply_block(el);

	el = create_rand_el(true);
	test_multiply_block(el);

	destroy_flash_matrix();
}
