#include <stdio.h>

#include "mem_dense_matrix.h"

using namespace fm;

int main()
{
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::COL);
	for (size_t i = 0; i < m1->get_num_rows(); i++) {
		for (size_t j = 0; j < m1->get_num_cols(); j++)
			m1->set(i, j, i * m1->get_num_cols() + j);
	}

	I_mem_dense_matrix::ptr m2 = I_mem_dense_matrix::create(10, 9,
			matrix_layout_t::COL);
	for (size_t i = 0; i < m2->get_num_rows(); i++) {
		for (size_t j = 0; j < m2->get_num_cols(); j++)
			m2->set(i, j, i * m2->get_num_cols() + j);
	}

	I_mem_dense_matrix::ptr m3 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::ROW);
	for (size_t i = 0; i < m3->get_num_rows(); i++) {
		for (size_t j = 0; j < m3->get_num_cols(); j++)
			m3->set(i, j, i * m3->get_num_cols() + j);
	}

	I_mem_dense_matrix::ptr res1 = multiply<int, int, int>(*m1, *m2);
	assert(res1->get_num_rows() == m1->get_num_rows());
	assert(res1->get_num_cols() == m2->get_num_cols());
	printf("The result matrix has %ld rows and %ld columns\n",
			res1->get_num_rows(), res1->get_num_cols());

	I_mem_dense_matrix::ptr res2 = multiply<int, int, int>(*m3, *m2);
	assert(res2->get_num_rows() == m3->get_num_rows());
	assert(res2->get_num_cols() == m2->get_num_cols());

	for (size_t i = 0; i < res1->get_num_rows(); i++) {
		for (size_t j = 0; j < res1->get_num_cols(); j++) {
			assert(res1->get(i, j) == res2->get(i, j));
		}
	}
}
