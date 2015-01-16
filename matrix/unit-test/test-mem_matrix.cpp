#include <stdio.h>

#include "mem_dense_matrix.h"

using namespace fm;

class set_col_operate: public type_set_operate<int>
{
	size_t num_cols;
public:
	set_col_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(int *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = (row_idx + i) * num_cols + col_idx;
		}
	}
};

class set_row_operate: public type_set_operate<int>
{
	size_t num_cols;
public:
	set_row_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(int *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = row_idx * num_cols + col_idx + i;
		}
	}
};

/*
 * This is a naive implementation of matrix multiplication.
 * It should be correct
 */
I_mem_dense_matrix::ptr naive_multiply(const I_mem_dense_matrix &m1,
		const I_mem_dense_matrix &m2)
{
	I_mem_dense_matrix::ptr res = I_mem_dense_matrix::create(
			m1.get_num_rows(), m2.get_num_cols(), matrix_layout_t::L_ROW);
	for (size_t i = 0; i < m1.get_num_rows(); i++) {
		for (size_t j = 0; j < m2.get_num_cols(); j++) {
			int sum = 0;
			for (size_t k = 0; k < m1.get_num_cols(); k++) {
				sum += m1.get(i, k) * m2.get(k, j);
			}
			res->set(i, j, sum);
		}
	}
	return res;
}

void verify_result(const I_mem_dense_matrix &m1, const I_mem_dense_matrix &m2)
{
	assert(m1.get_num_rows() == m2.get_num_rows());
	assert(m1.get_num_cols() == m2.get_num_cols());

	for (size_t i = 0; i < m1.get_num_rows(); i++)
		for (size_t j = 0; j < m1.get_num_cols(); j++)
			assert(m1.get(i, j) == m2.get(i, j));
}

void test1()
{
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	I_mem_dense_matrix::ptr m2 = I_mem_dense_matrix::create(10, 9,
			matrix_layout_t::L_COL, set_col_operate(9));
	I_mem_dense_matrix::ptr correct = naive_multiply(*m1, *m2);

	printf("Test multiply on col_matrix in on thread\n");
	I_mem_dense_matrix::ptr res1 = multiply<int, int, int>(*m1, *m2);
	verify_result(*res1, *correct);

	printf("Test multiply on col_matrix in parallel\n");
	res1 = par_multiply<int, int, int>(*m1, *m2);
	verify_result(*res1, *correct);
}

void test2()
{
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(10, 100,
			matrix_layout_t::L_ROW, set_row_operate(100));
	I_mem_dense_matrix::ptr m2 = I_mem_dense_matrix::create(100, 9,
			matrix_layout_t::L_COL, set_col_operate(9));
	I_mem_dense_matrix::ptr correct = naive_multiply(*m1, *m2);

	printf("Test multiply on row_matrix X col_matrix in on thread\n");
	I_mem_dense_matrix::ptr res1 = multiply<int, int, int>(*m1, *m2);
	verify_result(*res1, *correct);

	printf("Test multiply on row_matrix X col_matrix in parallel\n");
	res1 = par_multiply<int, int, int>(*m1, *m2);
	verify_result(*res1, *correct);
}

int main()
{
	test1();
	test2();
}
