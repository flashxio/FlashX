#include <stdio.h>

#include <boost/format.hpp>

#include "common.h"

#include "dense_matrix.h"
#include "mem_matrix_store.h"
#include "mem_worker_thread.h"
#include "local_matrix_store.h"
#include "sparse_matrix.h"

using namespace fm;

int num_nodes = -1;

class set_col_operate: public type_set_operate<double>
{
	size_t num_cols;
public:
	set_col_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(double *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = (row_idx + i) * num_cols + col_idx;
		}
	}
};

class set_row_operate: public type_set_operate<double>
{
	size_t num_cols;
public:
	set_row_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(double *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = row_idx * num_cols + col_idx + i;
		}
	}
};

class init_vec: public type_set_vec_operate<double>
{
public:
	virtual void set(double *arr, size_t num_eles, off_t start_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = start_idx + i;
	}
};

void test_MM_blas_chain(size_t nrow, size_t ncol, bool in_mem = true)
{
	printf("multiply matrix: M(%ld x %ld) * M(%ld %ld) twice\n",
			nrow, ncol, ncol, ncol);
	dense_matrix::ptr m1
		= dense_matrix::create(nrow, ncol, matrix_layout_t::L_COL,
				get_scalar_type<double>(), set_col_operate(ncol), num_nodes,
				in_mem);
	dense_matrix::ptr m2
		= dense_matrix::create(ncol, ncol, matrix_layout_t::L_COL,
				get_scalar_type<double>(), set_col_operate(ncol), num_nodes,
				in_mem);

	struct timeval start, end;
	{
		gettimeofday(&start, NULL);
		dense_matrix::ptr res1 = m1->multiply(*m2, matrix_layout_t::L_NONE, true);
		res1->materialize_self();
		gettimeofday(&end, NULL);
		printf("It takes %.3f seconds to multiply column matrix in parallel\n",
				time_diff(start, end));
	}

	printf("multiply two small matrices in a row\n");
	{
		gettimeofday(&start, NULL);
		dense_matrix::ptr res1 = m1->multiply(*m2, matrix_layout_t::L_NONE,
				true);
		dense_matrix::ptr res2 = res1->multiply(*m2, matrix_layout_t::L_NONE,
				true);
		res2->materialize_self();
		gettimeofday(&end, NULL);
		printf("It takes %.3f seconds to multiply column matrix in parallel\n",
				time_diff(start, end));
	}

	printf("multiply three small matrices in a row\n");
	{
		gettimeofday(&start, NULL);
		dense_matrix::ptr res1 = m1->multiply(*m2, matrix_layout_t::L_NONE,
				true);
		dense_matrix::ptr res2 = res1->multiply(*m2, matrix_layout_t::L_NONE, true);
		dense_matrix::ptr res3 = res2->multiply(*m2, matrix_layout_t::L_NONE, true);
		res3->materialize_self();
		gettimeofday(&end, NULL);
		printf("It takes %.3f seconds to multiply column matrix in parallel\n",
				time_diff(start, end));
	}

	printf("multiply four small matrices in a row\n");
	{
		gettimeofday(&start, NULL);
		dense_matrix::ptr res1 = m1->multiply(*m2, matrix_layout_t::L_NONE,
				true);
		dense_matrix::ptr res2 = res1->multiply(*m2, matrix_layout_t::L_NONE, true);
		dense_matrix::ptr res3 = res2->multiply(*m2, matrix_layout_t::L_NONE, true);
		dense_matrix::ptr res4 = res3->multiply(*m2, matrix_layout_t::L_NONE, true);
		res4->materialize_self();
		gettimeofday(&end, NULL);
		printf("It takes %.3f seconds to multiply column matrix in parallel\n",
				time_diff(start, end));
	}

	printf("multiply five small matrices in a row\n");
	{
		gettimeofday(&start, NULL);
		dense_matrix::ptr res1 = m1->multiply(*m2, matrix_layout_t::L_NONE,
				true);
		dense_matrix::ptr res2 = res1->multiply(*m2, matrix_layout_t::L_NONE, true);
		dense_matrix::ptr res3 = res2->multiply(*m2, matrix_layout_t::L_NONE, true);
		dense_matrix::ptr res4 = res3->multiply(*m2, matrix_layout_t::L_NONE, true);
		dense_matrix::ptr res5 = res4->multiply(*m2, matrix_layout_t::L_NONE, true);
		res5->materialize_self();
		gettimeofday(&end, NULL);
		printf("It takes %.3f seconds to multiply column matrix in parallel\n",
				time_diff(start, end));
	}
}

dense_matrix::ptr test_MM_blas(size_t nrow, size_t ncol, size_t right_ncol,
		bool in_mem = true)
{
	printf("multiply matrix: M(%ld x %ld) * M(%ld %ld)\n",
			nrow, ncol, ncol, right_ncol);

	struct timeval start, end;
	gettimeofday(&start, NULL);
	dense_matrix::ptr m1
		= dense_matrix::create(nrow, ncol, matrix_layout_t::L_COL,
				get_scalar_type<double>(), set_col_operate(ncol), num_nodes,
				in_mem);
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to construct input column matrix\n",
			time_diff(start, end));
	dense_matrix::ptr m2
		= dense_matrix::create(ncol, right_ncol, matrix_layout_t::L_COL,
				get_scalar_type<double>(), set_col_operate(right_ncol), num_nodes,
				in_mem);

	gettimeofday(&start, NULL);
	dense_matrix::ptr res1 = m1->multiply(*m2, matrix_layout_t::L_NONE, true);
	res1->materialize_self();
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to multiply column matrix in parallel\n",
			time_diff(start, end));
	return res1;
}

dense_matrix::ptr test_scale_cols(size_t nrow, size_t ncol, bool in_mem = true)
{
	printf("scale cols: M(%ld x %ld)\n", nrow, ncol);

	struct timeval start, end;
	dense_matrix::ptr m1
		= dense_matrix::create(nrow, ncol, matrix_layout_t::L_COL,
				get_scalar_type<double>(), set_col_operate(ncol), num_nodes,
				in_mem);
	vector::ptr v = vector::create(ncol, get_scalar_type<double>(), -1, true,
			init_vec());

	gettimeofday(&start, NULL);
	dense_matrix::ptr res1 = m1->scale_cols(v);
	res1->materialize_self();
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to scale cols\n", time_diff(start, end));
	return res1;
}

dense_matrix::ptr test_add(size_t nrow, size_t ncol, bool in_mem = true)
{
	printf("add matrix: M(%ld x %ld)\n", nrow, ncol);

	struct timeval start, end;
	dense_matrix::ptr m1
		= dense_matrix::create(nrow, ncol, matrix_layout_t::L_COL,
				get_scalar_type<double>(), set_col_operate(ncol), num_nodes,
				in_mem);
	dense_matrix::ptr m2
		= dense_matrix::create(nrow, ncol, matrix_layout_t::L_COL,
				get_scalar_type<double>(), set_col_operate(ncol), num_nodes,
				in_mem);

	gettimeofday(&start, NULL);
	dense_matrix::ptr res1 = m1->add(*m2);
	res1->materialize_self();
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to add matrix\n", time_diff(start, end));
	return res1;
}

dense_matrix::ptr test_scale_mat(size_t nrow, size_t ncol, bool in_mem = true)
{
	printf("scale matrix: M(%ld x %ld)\n", nrow, ncol);

	struct timeval start, end;
	dense_matrix::ptr m1
		= dense_matrix::create(nrow, ncol, matrix_layout_t::L_COL,
				get_scalar_type<double>(), set_col_operate(ncol), num_nodes,
				in_mem);

	gettimeofday(&start, NULL);
	dense_matrix::ptr res1 = m1->multiply_scalar<double>(2);
	res1->materialize_self();
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to scale matrix\n", time_diff(start, end));
	return res1;
}

void test_norm(size_t nrow, size_t ncol, bool in_mem = true)
{
	printf("norm: M(%ld x %ld)\n", nrow, ncol);

	struct timeval start, end;
	dense_matrix::ptr m1
		= dense_matrix::create(nrow, ncol, matrix_layout_t::L_COL,
				get_scalar_type<double>(), set_col_operate(ncol), num_nodes,
				in_mem);

	gettimeofday(&start, NULL);
	m1->col_norm2();
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to norm cols\n", time_diff(start, end));
}

void matrix_mul_tests(size_t long_dim)
{
	for (size_t short_dim = 4; short_dim <= 64; short_dim *= 2)
		test_MM_blas(short_dim, long_dim, short_dim);
	for (size_t short_dim = 4; short_dim <= 64; short_dim *= 2)
		test_MM_blas(long_dim, short_dim, short_dim);
	for (size_t short_dim = 4; short_dim <= 64; short_dim *= 2)
		test_MM_blas(short_dim, long_dim, short_dim, false);
	for (size_t short_dim = 4; short_dim <= 64; short_dim *= 2)
		test_MM_blas(long_dim, short_dim, short_dim, false);

	size_t short_dim = 16;
	test_MM_blas_chain(long_dim, short_dim);
	test_MM_blas_chain(long_dim, short_dim, false);
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
	if (matrix_conf.get_num_nodes() > 1)
		num_nodes = matrix_conf.get_num_nodes();

	size_t long_dim = 1024 * 1024 * 200;
	size_t short_dim = 16;
	matrix_mul_tests(long_dim);
	test_norm(long_dim, short_dim);
	test_scale_mat(long_dim, short_dim);
	test_add(long_dim, short_dim);
	test_scale_cols(long_dim, short_dim);
	test_norm(long_dim, short_dim, false);
	test_scale_mat(long_dim, short_dim, false);
	test_add(long_dim, short_dim, false);
	test_scale_cols(long_dim, short_dim, false);

	destroy_flash_matrix();
}
