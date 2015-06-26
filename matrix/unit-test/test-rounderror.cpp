#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <cblas.h>
#include <assert.h>
#include <string.h>

#include <vector>

#include "rand_gen.h"
#include "generic_type.h"
#include "mem_matrix_store.h"
#include "mem_vec_store.h"
#include "vector.h"
#include "dense_matrix.h"

using namespace fm;

void init_vec(std::vector<double> &vec)
{
	for (size_t i = 0; i < vec.size(); i++)
		vec[i] = ((double) random()) / ((double) random());
}

void mv(const std::vector<double> &A, const std::vector<double> &B,
		std::vector<double> &C)
{
	size_t num_rows = A.size() / B.size();
	size_t num_cols = B.size();
	for (size_t i = 0; i < num_rows; i++) {
		double sum = 0;
		for (size_t j = 0; j < num_cols; j++) {
			long double tmp1 = A[j * num_rows + i];
			long double tmp2 = B[j];
			sum += tmp2 * tmp1;
		}
		C[i] = sum;
	}
}

void multiply_col_row(const long double *A_col, int Alen, const long double *B_row,
		int Blen, long double *C_mat)
{
	for (int i = 0; i < Alen; i++)
		for (int j = 0; j < Blen; j++)
			C_mat[j * Alen + i] += A_col[i] * B_row[j];
}

/*
 * Multiplication on two column-wise matrices.
 */
void mm_cc(const std::vector<double> &A, int num_Arows, int num_Acols,
		const std::vector<double> &B, int num_Brows, int num_Bcols,
		std::vector<double> &C)
{
	std::vector<long double> C_tmp(C.size());

	assert(num_Acols == num_Brows);
	for (size_t i = 0; i < num_Acols; i++) {
		long double B_row[num_Bcols];
		for (size_t j = 0; j < num_Bcols; j++) {
			B_row[j] = B[j * num_Brows + i];
		}
		long double A_col[num_Arows];
		for (size_t j = 0; j < num_Arows; j++) {
			A_col[j] = A[i * num_Arows + j];
		}
		multiply_col_row(A_col, num_Arows, B_row, num_Bcols, C_tmp.data());
	}
	for (size_t i = 0; i < C.size(); i++)
		C[i] = C_tmp[i];
}

/*
 * Multiplication on a row-wise matrix and a column-wise matrix.
 */
void mm_rc(const std::vector<double> &A, int num_Arows, int num_Acols,
		const std::vector<double> &B, int num_Brows, int num_Bcols,
		std::vector<double> &C)
{
	assert(num_Acols == num_Brows);
	for (int i = 0; i < num_Arows; i++) {
		const double *Arow = &A[num_Acols * i];
		for (int j = 0; j < num_Bcols; j++) {
			const double *Bcol = &B[num_Brows * j];
			long double sum = 0;
			for (int k = 0; k < num_Acols; k++) {
				long double tmp1 = Arow[k];
				long double tmp2 = Bcol[k];
				sum += tmp1 * tmp2;
			}
			C[j * num_Arows + i] = sum;
		}
	}
}

void mm_cc1(const std::vector<double> &A, int num_Arows, int num_Acols,
		const std::vector<double> &B, int num_Brows, int num_Bcols,
		std::vector<double> &C)
{
	std::vector<double> B_buf(num_Brows);
	std::vector<double> C_buf(num_Arows);
	std::vector<double> C_buf1(num_Arows);
	for (int i = 0; i < num_Bcols; i++) {
		memcpy(B_buf.data(), B.data() + num_Brows * i,
				sizeof(double) * num_Brows);
		mv(A, B_buf, C_buf);
		memcpy(C.data() + num_Arows * i, C_buf.data(),
				sizeof(double) * num_Arows);

		cblas_dgemv(CblasColMajor, CblasNoTrans, num_Arows, num_Acols,
				1, A.data(), num_Arows, B_buf.data(), 1, 0, C_buf1.data(), 1);
		for (int i = 0; i < num_Arows; i++) {
			if (C_buf[i] != C_buf1[i])
				printf("%d: %g\n", i, C_buf[i] - C_buf1[i]);
			assert(C_buf[i] == C_buf1[i]);
		}
	}
}

void test_mv()
{
	printf("test mv\n");
	std::vector<double> A(100 * 10);
	std::vector<double> B(10);
	std::vector<double> C1(100);
	std::vector<double> C2(100);
	init_vec(A);
	init_vec(B);

	cblas_dgemv(CblasColMajor, CblasNoTrans, 100, 10,
			1, A.data(), 100, B.data(), 1, 0, C2.data(), 1);
	for (size_t i = 0; i < A.size(); i++)
		A[i] *= 1;
	for (size_t i = 0; i < B.size(); i++)
		B[i] *= 1;
	mv(A, B, C1);
	for (size_t i = 0; i < C1.size(); i++) {
		if (C1[i] - C2[i] != 0)
			printf("%ld: %g\n", i, C1[i] - C2[i]);
		assert(C1[i] == C2[i]);
	}
}

void copy_col2row(const std::vector<double> &mat1, std::vector<double> &mat2,
		int num_rows, int num_cols)
{
	for (int i = 0; i < num_rows; i++)
		for (int j = 0; j < num_cols; j++)
			mat2[i * num_cols + j] = mat1[num_rows * j + i];
}

void test_mm()
{
	size_t long_dim = 10000;
	printf("test mm\n");
	// Column-major
	std::vector<double> A(long_dim * 10);
	// Row-major
	std::vector<double> A_tmp(long_dim * 10);
	// Column-major
	std::vector<double> B(10 * 10);
	// Column-major
	std::vector<double> C1(long_dim * 10);
	// Column-major
	std::vector<double> C2(long_dim * 10);
	// Column-major
	std::vector<double> C3(long_dim * 10);
	rand_gen::ptr gen = get_scalar_type<double>().create_rand_gen(
			scalar_variable_impl<double>(0), scalar_variable_impl<double>(1));
	gen->gen(A.data(), A.size());
	gen->gen(B.data(), B.size());
	copy_col2row(A, A_tmp, long_dim, 10);

	detail::mem_col_matrix_store::ptr A_store
		= detail::mem_col_matrix_store::create(long_dim, 10,
				get_scalar_type<double>());
	for (size_t i = 0; i < A_store->get_num_cols(); i++)
		memcpy(A_store->get_col(i), A.data() + i * A_store->get_num_rows(),
				A_store->get_num_rows() * sizeof(double));
	dense_matrix::ptr A_mat = dense_matrix::create(A_store);

	detail::mem_row_matrix_store::ptr A_tmp_store
		= detail::mem_row_matrix_store::create(long_dim, 10,
				get_scalar_type<double>());
	for (size_t i = 0; i < A_tmp_store->get_num_rows(); i++)
		memcpy(A_tmp_store->get_row(i), A_tmp.data() + i * A_tmp_store->get_num_cols(),
				A_tmp_store->get_num_cols() * sizeof(double));
	dense_matrix::ptr A_tmp_mat = dense_matrix::create(A_tmp_store);

	detail::mem_col_matrix_store::ptr B_store
		= detail::mem_col_matrix_store::create(10, 10, get_scalar_type<double>());
	for (size_t i = 0; i < B_store->get_num_cols(); i++)
		memcpy(B_store->get_col(i), B.data() + i * B_store->get_num_rows(),
				B_store->get_num_rows() * sizeof(double));
	dense_matrix::ptr B_mat = dense_matrix::create(B_store);

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			long_dim, 10, 10, 1, A.data(), long_dim, B.data(), 10, 0, C2.data(),
			long_dim);
	for (size_t i = 0; i < A.size(); i++)
		A[i] *= 1;
	for (size_t i = 0; i < B.size(); i++)
		B[i] *= 1;
	for (size_t i = 0; i < A_tmp.size(); i++)
		A_tmp[i] *= 1;
	mm_cc(A, long_dim, 10, B, 10, 10, C1);
	mm_rc(A_tmp, long_dim, 10, B, 10, 10, C3);

	dense_matrix::ptr res = A_mat->multiply(*B_mat);
	assert(res->store_layout() == matrix_layout_t::L_COL);
	assert(res->get_type() == get_scalar_type<double>());
	res->materialize_self();
	const detail::mem_matrix_store &mem_res
		= dynamic_cast<const detail::mem_matrix_store &>(res->get_data());
	for (size_t i = 0; i < res->get_num_cols(); i++)
		for (size_t j = 0; j < res->get_num_rows(); j++)
			assert(mem_res.get<double>(j, i) == C1[i * res->get_num_rows() + j]);

	printf("There are %ld elements in C3\n", C3.size());
	for (size_t i = 0; i < C3.size(); i++) {
		if (C3[i] - C2[i] != 0)
			printf("%ld: %g (%g, %g)\n", i, C3[i] - C2[i], C3[i], C2[i]);
//		assert(C3[i] == C2[i]);
	}
	printf("There are %ld elements in C1\n", C1.size());
	for (size_t i = 0; i < C1.size(); i++) {
		if (C1[i] - C2[i] != 0)
			printf("%ld: %g (%g, %g)\n", i, C1[i] - C2[i], C1[i], C2[i]);
//		assert(C1[i] == C2[i]);
	}
}

void test_vec_scal()
{
	printf("test scale a vector\n");
	detail::mem_col_matrix_store::ptr mat_store
		= detail::mem_col_matrix_store::create(1000, 1,
				get_scalar_type<double>());
	rand_gen::ptr gen = mat_store->get_type().create_rand_gen(
			scalar_variable_impl<double>(0), scalar_variable_impl<double>(1));
	double *col = (double *) mat_store->get_col(0);
	gen->gen(mat_store->get_col(0), mat_store->get_num_rows());
	std::vector<double> copy1(mat_store->get_num_rows());
	memcpy(copy1.data(), col,
			mat_store->get_entry_size() * mat_store->get_num_rows());
	std::vector<double> copy2 = copy1;

	double scal = 0.001;
	// BLAS
	cblas_dscal(copy1.size(), scal, copy1.data(), 1);

	// Naive solution.
	for (size_t i = 0; i < copy2.size(); i++)
		copy2[i] = ((long double) copy2[i]) * ((long double) scal);
	for (size_t i = 0; i < copy2.size(); i++)
		assert(copy2[i] == copy1[i]);

	// FlashMatrix solution.
	dense_matrix::ptr mat = dense_matrix::create(mat_store);
	dense_matrix::ptr res = mat->multiply_scalar(scal);
	res->materialize_self();
	const detail::mem_matrix_store &mem_res
		= dynamic_cast<const detail::mem_matrix_store &>(res->get_data());
	for (size_t i = 0; i < res->get_num_rows(); i++)
		assert(mem_res.get<double>(i, 0) == copy1[i]);
}

void test_vec_dot()
{
	printf("test dot product on vectors\n");
	rand_gen::ptr gen = get_scalar_type<double>().create_rand_gen(
			scalar_variable_impl<double>(0), scalar_variable_impl<double>(1));
	std::vector<double> v1(1000);
	std::vector<double> v2(1000);
	gen->gen(v1.data(), v1.size());
	gen->gen(v2.data(), v2.size());

	// BLAS
	double res = cblas_ddot(v1.size(), v1.data(), 1, v2.data(), 1);

	// Naive solution.
	long double res1 = 0;
	for (size_t i = 0; i < v1.size(); i++)
		res1 += ((long double) v1[i]) * ((long double) v2[i]);
	printf("res diff: %g\n", ((double) res1) - res);
	assert(((double) res1) == res);
}

void test_norm2()
{
	printf("test norm2 on a vector\n");
	detail::mem_col_matrix_store::ptr mat_store
		= detail::mem_col_matrix_store::create(1000, 1,
				get_scalar_type<double>());
	rand_gen::ptr gen = mat_store->get_type().create_rand_gen(
			scalar_variable_impl<double>(0), scalar_variable_impl<double>(1));
	double *col = (double *) mat_store->get_col(0);
	gen->gen(col, mat_store->get_num_rows());

	// BLAS
	double res = cblas_dnrm2(mat_store->get_num_rows(), col, 1);

	// Naive solution.
	long double res1 = 0;
	for (size_t i = 0; i < mat_store->get_num_rows(); i++)
		res1 += ((long double) col[i]) * ((long double) col[i]);
	res1 = sqrtl(res1);
	printf("res diff: %g\n", ((double) res1) - res);
	assert(((double) res1) == res);

	// FlashMatrix solution.
	dense_matrix::ptr mat = dense_matrix::create(mat_store);
	double res2 = mat->norm2();
	assert(res2 == res);
}

int main()
{
	test_vec_scal();
	test_vec_dot();
	test_norm2();
	test_mv();
	test_mm();
}
