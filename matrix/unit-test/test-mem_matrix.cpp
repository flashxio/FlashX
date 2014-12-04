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

int main()
{
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	I_mem_dense_matrix::ptr m2 = I_mem_dense_matrix::create(10, 9,
			matrix_layout_t::L_COL, set_col_operate(9));
	I_mem_dense_matrix::ptr m3 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_ROW, set_row_operate(10));

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
