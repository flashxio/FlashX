#include "matrix_config.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"

using namespace fm;

size_t length = 4096L * 1024 * 1024;

class mat_init: public type_set_operate<int>
{
public:
	virtual void set(int *arr, size_t num_eles, off_t row_idx,
			            off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = i;
	}
};

class mat_init1: public type_set_operate<float>
{
	size_t num_vertices;
public:
	mat_init1(size_t num_vertices) {
		this->num_vertices = num_vertices;
	}

	virtual void set(float *arr, size_t num_eles, off_t row_idx,
			            off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = 1.0 / num_vertices;
	}
};

void test_mapply()
{
	struct timeval start, end;
	int num_nodes = matrix_conf.get_num_nodes();
	printf("mapply on %d nodes\n", num_nodes);
	gettimeofday(&start, NULL);
	detail::matrix_store::ptr tmp = detail::matrix_store::create(
			length, 1, matrix_layout_t::L_ROW,
			get_scalar_type<float>(), num_nodes, true);
	detail::matrix_store::ptr tmp1 = detail::matrix_store::create(
			length, 1, matrix_layout_t::L_ROW,
			get_scalar_type<float>(), num_nodes, true);
	tmp->set_data(mat_init());
	printf("init vec1 of %ld eles\n", tmp->get_num_rows());
	tmp1->set_data(mat_init1(tmp1->get_num_rows()));
	printf("init vec2 of %ld eles\n", tmp1->get_num_rows());
	gettimeofday(&end, NULL);
	printf("init takes %.3f seconds\n", time_diff(start, end));

	start = end;
	dense_matrix::ptr pr = dense_matrix::create(tmp);
	dense_matrix::ptr pr1 = dense_matrix::create(tmp1);
	pr = pr->multiply_scalar<float>(0.85);
	pr = pr->add_scalar<float>(0.15 / length);
	dense_matrix::ptr diff = pr1->minus(*pr);
	diff = diff->abs();
	gettimeofday(&end, NULL);
	printf("virtual computation takes %.3f seconds\n", time_diff(start, end));
	dense_matrix::ptr comp = diff->lt_scalar<float>(0.00001);
	scalar_variable::ptr convg = comp->sum();
	assert(convg->get_type() == get_scalar_type<size_t>());
	gettimeofday(&end, NULL);
	printf("test takes %.3f seconds\n", time_diff(start, end));
}

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
	else if (test_name == "mapply")
		test_mapply();
	else {
		fprintf(stderr, "unknow test\n");
		return -1;
	}
	destroy_flash_matrix();
}
