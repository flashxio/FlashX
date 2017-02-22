#include "sink_matrix.h"
#include "bulk_operate.h"
#include "agg_matrix_store.h"
#include "IPW_matrix_store.h"
#include "groupby_matrix_store.h"
#include "sparse_matrix.h"

using namespace fm;
using namespace fm::detail;

void test_transpose()
{
	matrix_store::const_ptr data = matrix_store::create(1000000, 10,
			matrix_layout_t::L_COL, get_scalar_type<double>(), -1, true);
	matrix_store::const_ptr data2 = matrix_store::create(1000000, 20,
			matrix_layout_t::L_COL, get_scalar_type<double>(), -1, true);
	matrix_store::const_ptr data3 = matrix_store::create(1000000, 4,
			matrix_layout_t::L_COL, get_scalar_type<double>(), -1, true);
	bulk_operate::const_ptr add = bulk_operate::conv2ptr(
			data->get_type().get_basic_ops().get_add());
	bulk_operate::const_ptr multiply = bulk_operate::conv2ptr(
			data->get_type().get_basic_ops().get_multiply());
	agg_operate::const_ptr agg_op = agg_operate::create(add, add);
	sink_store::const_ptr sink;
	matrix_store::const_ptr t;
	std::vector<matrix_store::const_ptr> sinks;

	// Agg matrix.
	sink = agg_matrix_store::create(data, matrix_margin::MAR_COL, agg_op);
	printf("agg sink: %ld, %ld\n", sink->get_num_rows(), sink->get_num_cols());
	assert(sink->get_num_rows() == 1);
	assert(sink->get_num_cols() == data->get_num_cols());
	t = sink->transpose();
	printf("t sink: %ld, %ld\n", t->get_num_rows(), t->get_num_cols());
	assert(t->get_num_rows() == data->get_num_cols());
	assert(t->get_num_cols() == 1);

	// Block agg matrix.
	sinks.resize(3);
	for (size_t i = 0; i < sinks.size(); i++)
		sinks[i] = sink;
	sink = block_sink_store::create(sinks, 1, sinks.size());
	printf("block sink: %ld, %ld\n", sink->get_num_rows(), sink->get_num_cols());
	assert(sink->get_num_rows() == 1);
	assert(sink->get_num_cols() == sinks.size() * data->get_num_cols());
	t = sink->transpose();
	printf("t sink: %ld, %ld\n", t->get_num_rows(), t->get_num_cols());
	assert(t->get_num_rows() == sinks.size() * data->get_num_cols());
	assert(t->get_num_cols() == 1);

	// IPW matrix.
	sink = IPW_matrix_store::create(data->transpose(), data2, multiply, add);
	printf("IPW sink: %ld, %ld\n", sink->get_num_rows(), sink->get_num_cols());
	assert(sink->get_num_rows() == data->get_num_cols());
	assert(sink->get_num_cols() == data2->get_num_cols());
	t = sink->transpose();
	printf("t sink: %ld, %ld\n", t->get_num_rows(), t->get_num_cols());
	assert(t->get_num_rows() == data2->get_num_cols());
	assert(t->get_num_cols() == data->get_num_cols());

	// Block IPW matrix.
	size_t num_block_rows = 3;
	size_t num_block_cols = 2;
	sinks.resize(num_block_rows * num_block_cols);
	for (size_t i = 0; i < num_block_rows; i++) {
		for (size_t j = 0; j < num_block_cols; j++) {
			matrix_store::const_ptr right = j == num_block_cols - 1 ? data3 : data2;
			sinks[num_block_cols * i + j] = IPW_matrix_store::create(
					data->transpose(), right, multiply, add);
		}
	}
	sink = block_sink_store::create(sinks, num_block_rows, num_block_cols);
	printf("block sink: %ld, %ld\n", sink->get_num_rows(), sink->get_num_cols());
	assert(sink->get_num_rows() == data->get_num_cols() * num_block_rows);
	assert(sink->get_num_cols() == (num_block_cols - 1) * data2->get_num_cols()
			+ data3->get_num_cols());
	t = sink->transpose();
	printf("t sink: %ld, %ld\n", t->get_num_rows(), t->get_num_cols());
	assert(t->get_num_rows() == (num_block_cols - 1) * data2->get_num_cols()
			+ data3->get_num_cols());
	assert(t->get_num_cols() == data->get_num_cols() * num_block_rows);

	// Groupby matrix
	factor_col_vector::const_ptr labels = factor_col_vector::create(factor(10),
			dense_matrix::create_randu<int>(0, 10, data->get_num_rows(), 1,
				matrix_layout_t::L_COL));
	sink = groupby_matrix_store::create(data, labels, matrix_margin::MAR_ROW,
			agg_op);
	printf("groupby sink: %ld, %ld\n", sink->get_num_rows(), sink->get_num_cols());
	assert(sink->get_num_rows() == 10);
	assert(sink->get_num_cols() == data->get_num_cols());
	t = sink->transpose();
	printf("t sink: %ld, %ld\n", t->get_num_rows(), t->get_num_cols());
	assert(t->get_num_rows() == data->get_num_cols());
	assert(t->get_num_cols() == 10);

	// Block groupby matrix.
	sinks.resize(3);
	for (size_t i = 0; i < sinks.size(); i++)
		sinks[i] = sink;
	sink = block_sink_store::create(sinks, 1, sinks.size());
	printf("block sink: %ld, %ld\n", sink->get_num_rows(), sink->get_num_cols());
	assert(sink->get_num_rows() == 10);
	assert(sink->get_num_cols() == sinks.size() * data->get_num_cols());
	t = sink->transpose();
	printf("t sink: %ld, %ld\n", t->get_num_rows(), t->get_num_cols());
	assert(t->get_num_rows() == sinks.size() * data->get_num_cols());
	assert(t->get_num_cols() == 10);
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

	test_transpose();

	destroy_flash_matrix();
}
