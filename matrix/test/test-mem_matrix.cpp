#include <stdio.h>

#include "common.h"

#include "mem_dense_matrix.h"
#include "mem_worker_thread.h"

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

/*
 * This multiplies a tall column-wise matrix with a small column-wise matrix.
 */
template<class Type>
mem_dense_matrix::ptr test_MM1(size_t nrow, size_t ncol, size_t right_ncol)
{
	struct timeval start, end;

	printf("inner product of tall col-wise matrix: M(%ld x %ld) * M(%ld %ld)\n",
			nrow, ncol, ncol, right_ncol);
	gettimeofday(&start, NULL);
	mem_dense_matrix::ptr m1
		= mem_dense_matrix::create(nrow, ncol, matrix_layout_t::L_COL,
				get_scalar_type<Type>(), set_col_operate(ncol), num_nodes);
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to construct input column matrix\n",
			time_diff(start, end));
	mem_dense_matrix::ptr m2
		= mem_dense_matrix::create(ncol, right_ncol, matrix_layout_t::L_COL,
				get_scalar_type<Type>(), set_col_operate(ncol), num_nodes);

	gettimeofday(&start, NULL);
	dense_matrix::ptr res1 = m1->multiply(*m2, m1->store_layout());
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to multiply column matrix in parallel\n",
			time_diff(start, end));
	assert(res1->get_num_rows() == m1->get_num_rows());
	assert(res1->get_num_cols() == m2->get_num_cols());
	printf("The result matrix has %ld rows and %ld columns\n",
			res1->get_num_rows(), res1->get_num_cols());

	return mem_dense_matrix::cast(res1);
}

/*
 * This multiplies a tall row-wise matrix with a small column-wise matrix.
 */
template<class Type>
mem_dense_matrix::ptr test_MM2(size_t nrow, size_t ncol, size_t right_ncol)
{
	struct timeval start, end;

	printf("inner product of tall row-wise matrix: M(%ld x %ld) * M(%ld %ld)\n",
			nrow, ncol, ncol, right_ncol);
	gettimeofday(&start, NULL);
	mem_dense_matrix::ptr m1 = mem_dense_matrix::create(nrow, ncol,
				matrix_layout_t::L_ROW, get_scalar_type<Type>(),
				set_row_operate(ncol), num_nodes);
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to construct input row matrix\n",
			time_diff(start, end));
	mem_dense_matrix::ptr m2 = mem_dense_matrix::create(ncol, right_ncol,
				matrix_layout_t::L_COL, get_scalar_type<Type>(),
				set_col_operate(ncol), num_nodes);

	gettimeofday(&start, NULL);
	dense_matrix::ptr res1 = m1->multiply(*m2, m1->store_layout());
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to multiply row matrix in parallel\n",
			time_diff(start, end));
	assert(res1->get_num_rows() == m1->get_num_rows());
	assert(res1->get_num_cols() == m2->get_num_cols());
	printf("The result matrix has %ld rows and %ld columns\n",
			res1->get_num_rows(), res1->get_num_cols());

	return mem_dense_matrix::cast(res1);
}

#if 0
/*
 * This multiplies a large column-wise matrix with a small column-wise matrix
 * directly. So this implementation potentially generates many CPU cache misses.
 */
template<class Type>
mem_dense_matrix::ptr test_MM3(size_t nrow, size_t ncol, size_t right_ncol)
{
	struct timeval start, end;

	printf("test tall col-wise matrix (bad impl): M(%ld x %ld) * M(%ld %ld)\n",
			nrow, ncol, ncol, right_ncol);
	gettimeofday(&start, NULL);
	mem_dense_matrix::ptr m1 = mem_dense_matrix::create(nrow, ncol,
				matrix_layout_t::L_COL, get_scalar_type<Type>(),
				set_col_operate(ncol), num_nodes);
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to construct input column matrix\n",
			time_diff(start, end));
	mem_dense_matrix::ptr m2 = mem_dense_matrix::create(ncol, right_ncol,
				matrix_layout_t::L_COL, get_scalar_type<Type>(),
				set_col_operate(ncol), num_nodes);
	
	mem_dense_matrix::ptr res_m = mem_dense_matrix::create(nrow, right_ncol,
				matrix_layout_t::L_COL, get_scalar_type<Type>(), num_nodes);

	gettimeofday(&start, NULL);
	const detail::mem_col_matrix_store &left_m
		= dynamic_cast<const detail::mem_col_matrix_store &>(m1->get_data());
	const detail::mem_col_matrix_store &right_m
		= dynamic_cast<const detail::mem_col_matrix_store &>(m2->get_data());
	const detail::mem_col_matrix_store &res_col_m
		= dynamic_cast<const detail::mem_col_matrix_store &>(res_m->get_data());
#pragma omp parallel for
	for (size_t i = 0; i < nrow; i++) {
		for (size_t j = 0; j < right_ncol; j++) {
			const Type *right_col = (const Type *) right_m.get_col(j);
			Type *out_col = (Type *) res_col_m.get_col(j);
			out_col[i] = 0;
			for (size_t k = 0; k < ncol; k++) {
				const Type *left_col = (const Type *) left_m.get_col(k);
				out_col[i] += left_col[i] * right_col[k];
			}
		}
	}
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to multiply column matrix in parallel\n",
			time_diff(start, end));
	return res_m;
}
#endif

template<class Type>
mem_dense_matrix::ptr test_MV1(size_t nrow, size_t ncol)
{
	printf("test a tall col-wise matrix: M(%ld x %ld) * v(%ld)\n",
			nrow, ncol, ncol);

	struct timeval start, end;
	gettimeofday(&start, NULL);
	mem_dense_matrix::ptr m1 = mem_dense_matrix::create(nrow, ncol,
				matrix_layout_t::L_COL, get_scalar_type<Type>(),
				set_col_operate(ncol), num_nodes);
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to construct input column matrix\n",
			time_diff(start, end));
	mem_dense_matrix::ptr m2 = mem_dense_matrix::create(ncol, 1,
				matrix_layout_t::L_COL, get_scalar_type<Type>(),
				set_col_operate(ncol), num_nodes);

	gettimeofday(&start, NULL);
	dense_matrix::ptr res1 = m1->multiply(*m2, m1->store_layout());
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to multiply column matrix in parallel\n",
			time_diff(start, end));
	assert(res1->get_num_rows() == m1->get_num_rows());
	assert(res1->get_num_cols() == m2->get_num_cols());
	printf("The result matrix has %ld rows and %ld columns\n",
			res1->get_num_rows(), res1->get_num_cols());
	return mem_dense_matrix::cast(res1);
}

template<class Type>
mem_dense_matrix::ptr test_MV2(size_t nrow, size_t ncol)
{
	printf("test a wide row-wise matrix: M(%ld x %ld) * v(%ld)\n",
			nrow, ncol, ncol);

	struct timeval start, end;
	gettimeofday(&start, NULL);
	mem_dense_matrix::ptr m1 = mem_dense_matrix::create(nrow, ncol,
				matrix_layout_t::L_ROW, get_scalar_type<Type>(),
				set_row_operate(ncol), num_nodes);
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to construct input row matrix\n",
			time_diff(start, end));
	mem_dense_matrix::ptr m2 = mem_dense_matrix::create(ncol, 1,
				matrix_layout_t::L_COL, get_scalar_type<Type>(),
				set_col_operate(ncol), num_nodes);

	gettimeofday(&start, NULL);
	dense_matrix::ptr res1 = m1->multiply(*m2, m1->store_layout());
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to multiply row matrix in parallel\n",
			time_diff(start, end));
	assert(res1->get_num_rows() == m1->get_num_rows());
	assert(res1->get_num_cols() == m2->get_num_cols());
	printf("The result matrix has %ld rows and %ld columns\n",
			res1->get_num_rows(), res1->get_num_cols());
	return mem_dense_matrix::cast(res1);
}

template<class Type>
void check_result(mem_dense_matrix::ptr m1, mem_dense_matrix::ptr m2)
{
	assert(m1->get_num_rows() == m2->get_num_rows());
	assert(m1->get_num_cols() == m2->get_num_cols());
#pragma omp parallel for
	for (size_t i = 0; i < m1->get_num_rows(); i++) {
		for (size_t j = 0; j < m1->get_num_cols(); j++) {
			assert(m1->get<Type>(i, j) == m2->get<Type>(i, j));
		}
	}
}

void matrix_mul_tests()
{
	size_t nrow = 1024 * 1024 * 124;
	size_t ncol = 20;
	printf("Multiplication of a large and tall matrix and a small square matrix\n");
	mem_dense_matrix::ptr res1;
	mem_dense_matrix::ptr res2;
	res1 = test_MM1<double>(nrow, ncol, ncol);
	res2 = test_MM2<double>(nrow, ncol, ncol);
	check_result<double>(res1, res2);
#if 0
	res2 = test_MM3<double>(nrow, ncol, ncol);
	check_result<double>(res1, res2);
#endif
}

void matrix_vec_mul_tests()
{
	size_t nrow = 1024 * 1024 * 124;
	size_t ncol = 120;
	printf("Multiplication of a large (tall/wide) matrix and a vector\n");
	test_MV1<double>(nrow, ncol);
	test_MV2<double>(ncol, nrow);
}

int main(int argc, char *argv[])
{
	int num_threads = 8;
	if (argc >= 3) {
		num_nodes = atoi(argv[1]);
		num_threads = atoi(argv[2]);
	}
	int real_num_nodes = num_nodes > 0 ? num_nodes : 1;
	detail::mem_thread_pool::init_global_mem_threads(real_num_nodes,
			num_threads / real_num_nodes);

	matrix_mul_tests();
	printf("\n\n");
	matrix_vec_mul_tests();
}
