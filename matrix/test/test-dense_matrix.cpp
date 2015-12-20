#include "matrix_config.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"

using namespace fm;

class mat_init: public type_set_operate<int>
{
public:
	virtual void set(int *arr, size_t num_eles, off_t row_idx,
			            off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = i;
	}
};

dense_matrix::ptr create_mat(size_t ncol, matrix_layout_t layout)
{
	struct timeval start, end;
	gettimeofday(&start, NULL);
	dense_matrix::ptr mat = dense_matrix::create(4L * 1024 * 1024 * 1024,
			ncol, layout, get_scalar_type<int>(), mat_init(),
			matrix_conf.get_num_nodes());
	gettimeofday(&end, NULL);
	printf("creating a %s-major matrix (%ld,%ld) takes %.3f seconds\n",
			layout == matrix_layout_t::L_ROW ? "row" : "col",
			mat->get_num_rows(), mat->get_num_cols(), time_diff(start, end));
	return mat;
}

void test_conv_layout()
{
}

void test_inner_prod()
{
}

void test_agg()
{
}

void test_agg_mat()
{
}

void test_groupby()
{
}

void test_mapply_rows()
{
}

void test_mapply_cols()
{
}

void test_mapply2()
{
}

void test_sapply()
{
}

void test_apply_scalar()
{
}

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "test conf_file test_name\n");
		return -1;
	}
	std::string conf_file = argv[1];
	std::string test_name = argv[2];
	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	if (test_name == "conv_layout")
		test_conv_layout();
	else if (test_name == "inner_prod")
		test_inner_prod();
	else if (test_name == "agg")
		test_agg();
	else if (test_name == "agg_mat")
		test_agg_mat();
	else if (test_name == "groupby")
		test_groupby();
	else if (test_name == "mapply_rows")
		test_mapply_rows();
	else if (test_name == "mapply_cols")
		test_mapply_cols();
	else if (test_name == "mapply2")
		test_mapply2();
	else if (test_name == "sapply")
		test_sapply();
	else if (test_name == "apply_scalar")
		test_apply_scalar();
	else {
		fprintf(stderr, "unknow test\n");
		return -1;
	}
	destroy_flash_matrix();
}
