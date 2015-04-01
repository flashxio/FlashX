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

void test_multiply_scalar()
{
	printf("Test scalar multiplication\n");
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	mem_dense_matrix::ptr orig = mem_dense_matrix::cast(m1->get_matrix());
	mem_dense_matrix::ptr res = mem_dense_matrix::cast(orig->multiply_scalar(10));
	for (size_t i = 0; i < res->get_num_rows(); i++)
		for (size_t j = 0; j < res->get_num_cols(); j++)
			assert(res->get<int>(i, j) == orig->get<int>(i, j) * 10);

	m1 = I_mem_dense_matrix::create(100, 10, matrix_layout_t::L_ROW,
			set_row_operate(10));
	orig = mem_dense_matrix::cast(m1->get_matrix());
	res = mem_dense_matrix::cast(orig->multiply_scalar(10));
	for (size_t i = 0; i < res->get_num_rows(); i++)
		for (size_t j = 0; j < res->get_num_cols(); j++)
			assert(res->get<int>(i, j) == orig->get<int>(i, j) * 10);
}

void test_ele_wise()
{
	printf("Test element-wise operations\n");
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	I_mem_dense_matrix::ptr m2 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	mem_dense_matrix::ptr res = mem_dense_matrix::cast(
			m1->get_matrix()->add(*m2->get_matrix()));
	for (size_t i = 0; i < res->get_num_rows(); i++)
		for (size_t j = 0; j < res->get_num_cols(); j++)
			assert(res->get<int>(i, j) == m1->get_matrix()->get<int>(i, j) * 2);
}

void test_multiply_col()
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

void test_agg_col()
{
	printf("Test aggregation on tall matrix stored column wise\n");
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	const bulk_operate &op
		= m1->get_matrix()->get_type().get_basic_ops().get_add();
	scalar_variable::ptr res = m1->get_matrix()->aggregate(op);
	assert(res->get_type() == m1->get_matrix()->get_type());
	int sum = *(int *) res->get_raw();
	int num_eles = m1->get_num_rows() * m1->get_num_cols();
	assert(sum == (num_eles - 1) * num_eles / 2);
}

void test_multiply_wide_row()
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

void test_multiply_tall_row()
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

void test_copy()
{
	printf("Test copy on column-wise matrices\n");
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL);
	m1->get_matrix()->reset_data();
	I_mem_dense_matrix::ptr m2 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	bool ret = m1->get_matrix()->copy_from(*m2->get_matrix());
	assert(ret);

	mem_dense_matrix::ptr col_m1
		= mem_dense_matrix::cast(m1->get_matrix());
	mem_dense_matrix::ptr col_m2
		= mem_dense_matrix::cast(m2->get_matrix());
	for (size_t i = 0; i < col_m1->get_num_rows(); i++)
		for (size_t j = 0; j < col_m1->get_num_cols(); j++)
			assert(col_m1->get<int>(i, j) == col_m2->get<int>(i, j));
}

void test_agg_row()
{
	printf("Test aggregation on tall matrix stored row wise\n");
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_ROW, set_row_operate(10));
	const bulk_operate &op
		= m1->get_matrix()->get_type().get_basic_ops().get_add();
	scalar_variable::ptr res = m1->get_matrix()->aggregate(op);
	assert(res->get_type() == m1->get_matrix()->get_type());
	int sum = *(int *) res->get_raw();
	int num_eles = m1->get_num_rows() * m1->get_num_cols();
	assert(sum == (num_eles - 1) * num_eles / 2);
}

void test_submatrix()
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

	std::vector<off_t> idxs2(2);
	idxs2[0] = 0;
	idxs2[1] = 1;
	mem_col_dense_matrix::ptr subsub_m = mem_col_dense_matrix::cast(sub_m->get_cols(idxs2));
	assert(subsub_m != NULL);
	assert(subsub_m->get_num_cols() == idxs2.size());
	for (size_t i = 0; i < idxs2.size(); i++)
		assert(memcmp(subsub_m->get_col(i), col_m->get_col(idxs[i]),
					col_m->get_entry_size() * col_m->get_num_rows()) == 0);
}

void test_copy_sub()
{
	printf("Test copy on column-wise sub-matrices\n");
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 3,
			matrix_layout_t::L_COL);
	mem_dense_matrix::ptr mem_m1
		= mem_dense_matrix::cast(m1->get_matrix());
	mem_m1->reset_data();
	I_mem_dense_matrix::ptr m2 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL, set_col_operate(10));

	std::vector<off_t> idxs(3);
	idxs[0] = 1;
	idxs[1] = 5;
	idxs[2] = 3;
	mem_dense_matrix::ptr sub_m = mem_dense_matrix::cast(
			mem_col_dense_matrix::cast(m2->get_matrix())->get_cols(idxs));
	bool ret = mem_m1->copy_from(*sub_m);
	assert(ret);
	assert(sub_m->get_num_cols() == m1->get_num_cols());

	for (size_t i = 0; i < mem_m1->get_num_rows(); i++)
		for (size_t j = 0; j < mem_m1->get_num_cols(); j++)
			assert(mem_m1->get<int>(i, j) == sub_m->get<int>(i, j));

	mem_dense_matrix::ptr copy_m = mem_dense_matrix::cast(sub_m->deep_copy());
	assert(copy_m->get_num_cols() == idxs.size());
	for (size_t i = 0; i < copy_m->get_num_rows(); i++)
		for (size_t j = 0; j < copy_m->get_num_cols(); j++)
			assert(copy_m->get<int>(i, j) == sub_m->get<int>(i, j));
}

void test_set_data_sub()
{
	printf("Test set data on a column-wise sub-matrix\n");
	I_mem_dense_matrix::ptr m = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	std::vector<off_t> idxs(3);
	idxs[0] = 1;
	idxs[1] = 5;
	idxs[2] = 3;
	mem_dense_matrix::ptr sub_m = mem_dense_matrix::cast(
			mem_col_dense_matrix::cast(m->get_matrix())->get_cols(idxs));
	sub_m->reset_data();
	for (size_t i = 0; i < sub_m->get_num_rows(); i++)
		for (size_t j = 0; j < sub_m->get_num_cols(); i++)
			assert(sub_m->get<int>(i, j) == 0);
	sub_m->set_data(set_col_operate(sub_m->get_num_cols()));
	int val = 0;
	for (size_t i = 0; i < sub_m->get_num_rows(); i++)
		for (size_t j = 0; j < sub_m->get_num_cols(); i++)
			assert(sub_m->get<int>(i, j) == val++);
}

void test_agg_sub_col()
{
	printf("Test aggregation on a column-wise submatrix\n");
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	mem_col_dense_matrix::ptr col_m
		= mem_col_dense_matrix::cast(m1->get_matrix());
	std::vector<off_t> idxs(3);
	idxs[0] = 1;
	idxs[1] = 5;
	idxs[2] = 3;
	mem_col_dense_matrix::ptr sub_m
		= mem_col_dense_matrix::cast(col_m->get_cols(idxs));
	assert(sub_m != NULL);

	const bulk_operate &op = sub_m->get_type().get_basic_ops().get_add();
	scalar_variable::ptr res = sub_m->aggregate(op);
	assert(res->get_type() == sub_m->get_type());
	size_t sum = *(int *) res->get_raw();
	size_t ncol = m1->get_num_cols();
	size_t nrow = m1->get_num_rows();
	size_t sub_ncol = sub_m->get_num_cols();
	size_t expected = sub_ncol * ncol * (nrow - 1) * nrow / 2;
	for (size_t i = 0; i < idxs.size(); i++)
		expected += idxs[i] * nrow;
	assert(sum == expected);
}

void test_agg_sub_row()
{
	printf("Test aggregation on a row-wise submatrix\n");
	I_mem_dense_matrix::ptr m1 = I_mem_dense_matrix::create(100, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	mem_col_dense_matrix::ptr col_m
		= mem_col_dense_matrix::cast(m1->get_matrix());
	std::vector<off_t> idxs(3);
	idxs[0] = 1;
	idxs[1] = 5;
	idxs[2] = 3;
	mem_col_dense_matrix::ptr sub_col_m
		= mem_col_dense_matrix::cast(col_m->get_cols(idxs));
	mem_row_dense_matrix::ptr sub_row_m
		= mem_row_dense_matrix::cast(sub_col_m->transpose());

	const bulk_operate &op = sub_col_m->get_type().get_basic_ops().get_add();
	scalar_variable::ptr col_res = sub_col_m->aggregate(op);
	assert(col_res->get_type() == sub_col_m->get_type());
	scalar_variable::ptr row_res = sub_row_m->aggregate(op);
	assert(row_res->get_type() == sub_row_m->get_type());
	assert(*(int *) col_res->get_raw() == *(int *) row_res->get_raw());
}

void test_conv_row_col()
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
	test_multiply_scalar();
	test_ele_wise();
	test_multiply_col();
	test_agg_col();
	test_multiply_wide_row();
	test_multiply_tall_row();
	test_copy();
	test_agg_row();
	test_submatrix();
	test_copy_sub();
	test_agg_sub_col();
	test_agg_sub_row();
	test_conv_row_col();
}
