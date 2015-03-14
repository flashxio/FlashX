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
	printf("Test multiplication on tall matrix stored column wise\n");
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	I_mem_dense_matrix::ptr m2 = I_mem_dense_matrix::create(10, 9,
			matrix_layout_t::L_COL, set_col_operate(9));
	I_mem_dense_matrix::ptr correct = naive_multiply(*m1, *m2);

	printf("Test multiply on col_matrix in one thread\n");
	I_mem_dense_matrix::ptr res1 = multiply<int, int, int>(*m1, *m2);
	verify_result(*res1, *correct);

	printf("Test multiply on col_matrix in parallel\n");
	res1 = par_multiply<int, int, int>(*m1, *m2);
	verify_result(*res1, *correct);
}

void test2()
{
	printf("Test multiplication on wide matrix stored row wise\n");
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(10, 100,
			matrix_layout_t::L_ROW, set_row_operate(100));
	I_mem_dense_matrix::ptr m2 = I_mem_dense_matrix::create(100, 9,
			matrix_layout_t::L_COL, set_col_operate(9));
	I_mem_dense_matrix::ptr correct = naive_multiply(*m1, *m2);

	printf("Test multiply on row_matrix X col_matrix in one thread\n");
	I_mem_dense_matrix::ptr res1 = multiply<int, int, int>(*m1, *m2);
	verify_result(*res1, *correct);

	printf("Test multiply on row_matrix X col_matrix in parallel\n");
	res1 = par_multiply<int, int, int>(*m1, *m2);
	verify_result(*res1, *correct);
}

void test3()
{
	printf("Test multiplication on tall matrix stored row wise\n");
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_ROW, set_col_operate(10));
	I_mem_dense_matrix::ptr m2 = I_mem_dense_matrix::create(10, 9,
			matrix_layout_t::L_COL, set_col_operate(9));
	I_mem_dense_matrix::ptr correct = naive_multiply(*m1, *m2);

	printf("Test multiply on row_matrix in one thread\n");
	I_mem_dense_matrix::ptr res1 = multiply<int, int, int>(*m1, *m2);
	verify_result(*res1, *correct);

	printf("Test multiply on row_matrix in parallel\n");
	res1 = par_multiply<int, int, int>(*m1, *m2);
	verify_result(*res1, *correct);
}

void test4()
{
	printf("test submatrix of a column-wise matrix\n");
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	mem_col_dense_matrix::ptr col_m
		= mem_col_dense_matrix::cast(m1->get_matrix());
	std::vector<off_t> idxs(3);
	idxs[0] = 1;
	idxs[1] = 5;
	idxs[2] = 3;
	mem_col_dense_matrix::ptr sub_m = mem_col_dense_matrix::cast(col_m->get_cols(idxs));
	assert(sub_m != NULL);
	assert(sub_m->get_num_rows() == col_m->get_num_rows());
	assert(sub_m->get_num_cols() == idxs.size());
	assert(sub_m->store_layout() == col_m->store_layout());
	assert(sub_m->get_entry_size() == col_m->get_entry_size());
	assert(sub_m->get_type() == col_m->get_type());
	for (size_t i = 0; i < idxs.size(); i++)
		assert(memcmp(sub_m->get_col(i), col_m->get_col(idxs[i]),
				sub_m->get_entry_size() * sub_m->get_num_rows()) == 0);
}

void test5()
{
	printf("test conv col-wise to row-wise, row-wise to col-wise\n");
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(10000, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	mem_col_dense_matrix::ptr col_m
		= mem_col_dense_matrix::cast(m1->get_matrix());
	mem_row_dense_matrix::ptr row_m = col_m->get_row_store();
	I_mem_dense_matrix::ptr c1 = I_mem_dense_matrix::create(col_m);
	I_mem_dense_matrix::ptr c2 = I_mem_dense_matrix::create(row_m);
	for (size_t i = 0; i < col_m->get_num_rows(); i++) {
		for (size_t j = 0; j < col_m->get_num_cols(); j++)
			assert(c1->get(i, j) == c2->get(i, j));
	}
	col_m = row_m->get_col_store();
	c1 = I_mem_dense_matrix::create(col_m);
	c2 = I_mem_dense_matrix::create(row_m);
	for (size_t i = 0; i < col_m->get_num_rows(); i++)
		for (size_t j = 0; j < col_m->get_num_cols(); j++)
			assert(c1->get(i, j) == c2->get(i, j));
}

int main()
{
	test1();
	test2();
	test3();
	test4();
	test5();
}
