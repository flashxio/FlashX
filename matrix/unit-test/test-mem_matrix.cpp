#include <stdio.h>
#include <cblas.h>

#include "mem_dense_matrix.h"
#include "mem_vector.h"
#include "mem_worker_thread.h"

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

class set_col_long_operate: public type_set_operate<size_t>
{
	size_t num_cols;
public:
	set_col_long_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(size_t *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = (row_idx + i) * num_cols + col_idx;
		}
	}
};

class set_row_long_operate: public type_set_operate<size_t>
{
	size_t num_cols;
public:
	set_row_long_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(size_t *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = row_idx * num_cols + col_idx + i;
		}
	}
};

size_t long_dim = 10000000;

/*
 * This is a naive implementation of matrix multiplication.
 * It should be correct
 */
mem_dense_matrix::ptr naive_multiply(const mem_dense_matrix &m1,
		const mem_dense_matrix &m2)
{
	mem_dense_matrix::ptr res = mem_dense_matrix::create(
			m1.get_num_rows(), m2.get_num_cols(), matrix_layout_t::L_ROW,
			get_scalar_type<int>());
	detail::mem_matrix_store &res_store
		= (detail::mem_matrix_store &) res->get_data();
#pragma omp parallel for
	for (size_t i = 0; i < m1.get_num_rows(); i++) {
		for (size_t j = 0; j < m2.get_num_cols(); j++) {
			int sum = 0;
			for (size_t k = 0; k < m1.get_num_cols(); k++) {
				sum += m1.get<int>(i, k) * m2.get<int>(k, j);
			}
			res_store.set<int>(i, j, sum);
		}
	}
	return res;
}

void verify_result(const mem_dense_matrix &m1, const mem_dense_matrix &m2)
{
	assert(m1.get_num_rows() == m2.get_num_rows());
	assert(m1.get_num_cols() == m2.get_num_cols());

#pragma omp parallel for
	for (size_t i = 0; i < m1.get_num_rows(); i++)
		for (size_t j = 0; j < m1.get_num_cols(); j++)
			assert(m1.get<int>(i, j) == m2.get<int>(i, j));
}

void test_multiply_scalar(int num_nodes)
{
	printf("Test scalar multiplication\n");
	mem_dense_matrix::ptr orig = mem_dense_matrix::create(long_dim, 10,
			matrix_layout_t::L_COL, get_scalar_type<int>(), set_col_operate(10),
			num_nodes);
	mem_dense_matrix::ptr res = mem_dense_matrix::cast(orig->multiply_scalar(10));
	detail::mem_matrix_store &orig_store1
		= (detail::mem_matrix_store &) orig->get_data();
	detail::mem_matrix_store &res_store1
		= (detail::mem_matrix_store &) res->get_data();
#pragma omp parallel for
	for (size_t i = 0; i < res_store1.get_num_rows(); i++)
		for (size_t j = 0; j < res_store1.get_num_cols(); j++)
			assert(res_store1.get<int>(i, j) == orig_store1.get<int>(i, j) * 10);

	orig = mem_dense_matrix::create(long_dim, 10, matrix_layout_t::L_ROW,
			get_scalar_type<int>(), set_row_operate(10), num_nodes);
	res = mem_dense_matrix::cast(orig->multiply_scalar(10));
	detail::mem_matrix_store &orig_store2
		= (detail::mem_matrix_store &) orig->get_data();
	detail::mem_matrix_store &res_store2
		= (detail::mem_matrix_store &) res->get_data();
#pragma omp parallel for
	for (size_t i = 0; i < res_store2.get_num_rows(); i++)
		for (size_t j = 0; j < res_store2.get_num_cols(); j++)
			assert(res_store2.get<int>(i, j) == orig_store2.get<int>(i, j) * 10);
}

void test_ele_wise(int num_nodes)
{
	printf("Test element-wise operations\n");
	mem_dense_matrix::ptr m1 = mem_dense_matrix::create(long_dim, 10,
			matrix_layout_t::L_COL, get_scalar_type<int>(), set_col_operate(10),
			num_nodes);
	mem_dense_matrix::ptr m2 = mem_dense_matrix::create(long_dim, 10,
			matrix_layout_t::L_COL, get_scalar_type<int>(), set_col_operate(10),
			num_nodes);
	mem_dense_matrix::ptr res = mem_dense_matrix::cast(m1->add(*m2));
	detail::mem_matrix_store &res_store = (detail::mem_matrix_store &) res->get_data();
	detail::mem_matrix_store &m1_store
		= (detail::mem_matrix_store &) m1->get_data();
#pragma omp parallel for
	for (size_t i = 0; i < res_store.get_num_rows(); i++)
		for (size_t j = 0; j < res_store.get_num_cols(); j++)
			assert(res_store.get<int>(i, j) == m1_store.get<int>(i, j) * 2);
}

void test_multiply_col(int num_nodes)
{
	printf("Test multiplication on tall matrix stored column wise\n");
	mem_dense_matrix::ptr m1 = mem_dense_matrix::create(long_dim, 10,
			matrix_layout_t::L_COL, get_scalar_type<int>(), set_col_operate(10),
			num_nodes);
	mem_dense_matrix::ptr m2 = mem_dense_matrix::create(10, 9,
			matrix_layout_t::L_COL, get_scalar_type<int>(), set_col_operate(9),
			num_nodes);
	mem_dense_matrix::ptr correct = naive_multiply(*m1, *m2);

	printf("Test multiply on col_matrix\n");
	mem_dense_matrix::ptr res1 = mem_dense_matrix::cast(m1->multiply(*m2,
				m1->store_layout()));
	verify_result(*res1, *correct);
}

void test_agg_col(int num_nodes)
{
	printf("Test aggregation on tall matrix stored column wise\n");
	mem_dense_matrix::ptr m1 = mem_dense_matrix::create(long_dim, 10,
			matrix_layout_t::L_COL, get_scalar_type<size_t>(),
			set_col_long_operate(10), num_nodes);
	const bulk_operate &op
		= m1->get_type().get_basic_ops().get_add();
	scalar_variable::ptr res = m1->aggregate(op);
	assert(res->get_type() == m1->get_type());
	size_t sum = *(size_t *) res->get_raw();
	size_t num_eles = m1->get_num_rows() * m1->get_num_cols();
	assert(sum == (num_eles - 1) * num_eles / 2);
}

void test_multiply_wide_row(int num_nodes)
{
	printf("Test multiplication on wide matrix stored row wise\n");
	mem_dense_matrix::ptr m1 = mem_dense_matrix::create(10, long_dim,
			matrix_layout_t::L_ROW, get_scalar_type<int>(),
			set_row_operate(long_dim), num_nodes);
	mem_dense_matrix::ptr m2 = mem_dense_matrix::create(long_dim, 9,
			matrix_layout_t::L_COL, get_scalar_type<int>(), set_col_operate(9),
			num_nodes);
	mem_dense_matrix::ptr correct = naive_multiply(*m1, *m2);

	printf("Test multiply on row_matrix X col_matrix\n");
	mem_dense_matrix::ptr res1 = mem_dense_matrix::cast(m1->multiply(*m2,
				m1->store_layout()));
	verify_result(*res1, *correct);
}

void test_multiply_tall_row(int num_nodes)
{
	printf("Test multiplication on tall matrix stored row wise\n");
	mem_dense_matrix::ptr m1 = mem_dense_matrix::create(long_dim, 10,
			matrix_layout_t::L_ROW, get_scalar_type<int>(), set_col_operate(10),
			num_nodes);
	mem_dense_matrix::ptr m2 = mem_dense_matrix::create(10, 9,
			matrix_layout_t::L_COL, get_scalar_type<int>(), set_col_operate(9),
			num_nodes);
	mem_dense_matrix::ptr correct = naive_multiply(*m1, *m2);

	printf("Test multiply on row_matrix\n");
	mem_dense_matrix::ptr res1 = mem_dense_matrix::cast(m1->multiply(*m2,
				m1->store_layout()));
	verify_result(*res1, *correct);
}

void test_agg_row(int num_nodes)
{
	printf("Test aggregation on tall matrix stored row wise\n");
	mem_dense_matrix::ptr m1 = mem_dense_matrix::create(long_dim, 10,
			matrix_layout_t::L_ROW, get_scalar_type<size_t>(),
			set_row_long_operate(10), num_nodes);
	const bulk_operate &op
		= m1->get_type().get_basic_ops().get_add();
	scalar_variable::ptr res = m1->aggregate(op);
	assert(res->get_type() == m1->get_type());
	size_t sum = *(size_t *) res->get_raw();
	size_t num_eles = m1->get_num_rows() * m1->get_num_cols();
	assert(sum == (num_eles - 1) * num_eles / 2);
}

void test_agg_sub_col(int num_nodes)
{
	printf("Test aggregation on a column-wise submatrix\n");
	mem_dense_matrix::ptr col_m = mem_dense_matrix::create(long_dim, 10,
			matrix_layout_t::L_COL, get_scalar_type<size_t>(),
			set_col_long_operate(10), num_nodes);
	std::vector<off_t> idxs(3);
	idxs[0] = 1;
	idxs[1] = 5;
	idxs[2] = 3;
	mem_dense_matrix::ptr sub_m = mem_dense_matrix::cast(col_m->get_cols(idxs));
	assert(sub_m != NULL);

	const bulk_operate &op = sub_m->get_type().get_basic_ops().get_add();
	scalar_variable::ptr res = sub_m->aggregate(op);
	assert(res->get_type() == sub_m->get_type());
	size_t sum = *(size_t *) res->get_raw();
	size_t ncol = col_m->get_num_cols();
	size_t nrow = col_m->get_num_rows();
	size_t sub_ncol = sub_m->get_num_cols();
	size_t expected = sub_ncol * ncol * (nrow - 1) * nrow / 2;
	for (size_t i = 0; i < idxs.size(); i++)
		expected += idxs[i] * nrow;
	assert(sum == expected);
}

void test_agg_sub_row(int num_nodes)
{
	printf("Test aggregation on a row-wise submatrix\n");
	mem_dense_matrix::ptr col_m = mem_dense_matrix::create(long_dim, 10,
			matrix_layout_t::L_COL, get_scalar_type<size_t>(),
			set_col_long_operate(10), num_nodes);
	std::vector<off_t> idxs(3);
	idxs[0] = 1;
	idxs[1] = 5;
	idxs[2] = 3;
	mem_dense_matrix::ptr sub_col_m
		= mem_dense_matrix::cast(col_m->get_cols(idxs));
	mem_dense_matrix::ptr sub_row_m
		= mem_dense_matrix::cast(sub_col_m->transpose());

	const bulk_operate &op = sub_col_m->get_type().get_basic_ops().get_add();
	scalar_variable::ptr col_res = sub_col_m->aggregate(op);
	assert(col_res->get_type() == sub_col_m->get_type());
	scalar_variable::ptr row_res = sub_row_m->aggregate(op);
	assert(row_res->get_type() == sub_row_m->get_type());
	assert(*(size_t *) col_res->get_raw() == *(size_t *) row_res->get_raw());
}

#if 0

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

#endif

void test_rand_init()
{
	printf("test rand init\n");
	mem_dense_matrix::ptr m = mem_dense_matrix::create_rand<double>(-1.0, 1.0,
			long_dim / 100, 10, matrix_layout_t::L_COL);
	double sum = 0;
	for (size_t i = 0; i < m->get_num_rows(); i++)
		for (size_t j = 0; j < m->get_num_cols(); j++) {
			double v = m->get<double>(i, j);
			assert(v >= -1.0 && v <= 1.0);
			sum += v;
		}
	printf("sum: %f\n", sum);
}

#if 0

void test_flatten()
{
	printf("test flatten a matrix to a vector\n");
	I_mem_dense_matrix::ptr m = I_mem_dense_matrix::create(10000, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	mem_vector::ptr vec = m->get_matrix()->flatten(true);
	for (size_t i = 0; i < vec->get_length(); i++)
		assert((size_t) vec->get<int>(i) == i);

	m = I_mem_dense_matrix::create(10000, 10, matrix_layout_t::L_ROW,
			set_row_operate(10));
	vec = m->get_matrix()->flatten(true);
	for (size_t i = 0; i < vec->get_length(); i++)
		assert((size_t) vec->get<int>(i) == i);

	m = I_mem_dense_matrix::create(10000, 10, matrix_layout_t::L_COL,
			set_col_operate(10));
	vec = m->get_matrix()->flatten(false);
	for (size_t i = 0; i < vec->get_length(); i++) {
		size_t row_idx = i % m->get_num_rows();
		size_t col_idx = i / m->get_num_rows();
		assert((size_t) vec->get<int>(i) == row_idx * m->get_num_cols() + col_idx);
	}

	m = I_mem_dense_matrix::create(10000, 10, matrix_layout_t::L_ROW,
			set_row_operate(10));
	vec = m->get_matrix()->flatten(false);
	for (size_t i = 0; i < vec->get_length(); i++) {
		size_t row_idx = i % m->get_num_rows();
		size_t col_idx = i / m->get_num_rows();
		assert((size_t) vec->get<int>(i) == row_idx * m->get_num_cols() + col_idx);
	}
}

class time2_apply_operate: public arr_apply_operate
{
public:
	time2_apply_operate(size_t num_out_eles): arr_apply_operate(num_out_eles) {
	}
	virtual void run(const mem_vector &in, mem_vector &out) const {
		assert(out.get_length() == get_num_out_eles());
		for (size_t i = 0; i < get_num_out_eles(); i++)
			out.set<int>(i, in.get<int>(i) * 2);
	}
	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<int>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<int>();
	}
};

void test_apply()
{
	printf("test applying to a matrix\n");
	I_mem_dense_matrix::ptr m = I_mem_dense_matrix::create(long_dim, 10,
			matrix_layout_t::L_COL, set_col_operate(10));
	size_t out_num_cols = m->get_num_cols() / 2;
	mem_dense_matrix::ptr res = mem_dense_matrix::cast(
			m->get_matrix()->apply(apply_margin::MAR_ROW,
			time2_apply_operate(out_num_cols)));
	assert(res->get_num_rows() == m->get_num_rows());
	assert(res->get_num_cols() == out_num_cols);
	for (size_t i = 0; i < res->get_num_rows(); i++) {
		for (size_t j = 0; j < res->get_num_cols(); j++)
			assert(m->get(i, j) * 2 == res->get<int>(i, j));
	}

	m = I_mem_dense_matrix::create(long_dim, 10, matrix_layout_t::L_ROW,
			set_row_operate(10));
	res = mem_dense_matrix::cast(m->get_matrix()->apply(apply_margin::MAR_ROW,
				time2_apply_operate(out_num_cols)));
	printf("output matrix: %ld, %ld\n", res->get_num_rows(), res->get_num_cols());
	assert(res->get_num_rows() == m->get_num_rows());
	assert(res->get_num_cols() == out_num_cols);
	for (size_t i = 0; i < res->get_num_rows(); i++) {
		for (size_t j = 0; j < res->get_num_cols(); j++)
			assert(m->get(i, j) * 2 == res->get<int>(i, j));
	}

	m = I_mem_dense_matrix::create(long_dim, 10, matrix_layout_t::L_COL,
			set_col_operate(10));
	size_t out_num_rows = m->get_num_rows() / 2;
	res = mem_dense_matrix::cast(m->get_matrix()->apply(apply_margin::MAR_COL,
			time2_apply_operate(out_num_rows)));
	assert(res->get_num_rows() == out_num_rows);
	assert(res->get_num_cols() == m->get_num_cols());
	for (size_t i = 0; i < res->get_num_rows(); i++) {
		for (size_t j = 0; j < res->get_num_cols(); j++)
			assert(m->get(i, j) * 2 == res->get<int>(i, j));
	}

	m = I_mem_dense_matrix::create(long_dim, 10, matrix_layout_t::L_ROW,
			set_row_operate(10));
	res = mem_dense_matrix::cast(m->get_matrix()->apply(apply_margin::MAR_COL,
				time2_apply_operate(out_num_rows)));
	assert(res->get_num_rows() == out_num_rows);
	assert(res->get_num_cols() == m->get_num_cols());
	for (size_t i = 0; i < res->get_num_rows(); i++) {
		for (size_t j = 0; j < res->get_num_cols(); j++)
			assert(m->get(i, j) * 2 == res->get<int>(i, j));
	}
}

#endif

void test_scale_cols(int num_nodes)
{
	printf("Test scale cols\n");
	mem_dense_matrix::ptr orig = mem_dense_matrix::create(long_dim, 10,
			matrix_layout_t::L_COL, get_scalar_type<int>(), set_col_operate(10),
			num_nodes);
	mem_vector::ptr vals = mem_vector::create(orig->get_num_cols(),
			orig->get_type());
	for (size_t i = 0; i < vals->get_length(); i++)
		vals->set<int>(i, i);
	mem_dense_matrix::ptr res = mem_dense_matrix::cast(orig->scale_cols(*vals));
	detail::mem_matrix_store &orig_store1
		= (detail::mem_matrix_store &) orig->get_data();
	detail::mem_matrix_store &res_store1
		= (detail::mem_matrix_store &) res->get_data();
#pragma omp parallel for
	for (size_t i = 0; i < res_store1.get_num_rows(); i++)
		for (size_t j = 0; j < res_store1.get_num_cols(); j++)
			assert(res_store1.get<int>(i, j)
					== orig_store1.get<int>(i, j) * vals->get<int>(j));

	orig = mem_dense_matrix::create(long_dim, 10, matrix_layout_t::L_ROW,
			get_scalar_type<int>(), set_row_operate(10), num_nodes);
	res = mem_dense_matrix::cast(orig->scale_cols(*vals));
	detail::mem_matrix_store &orig_store2
		= (detail::mem_matrix_store &) orig->get_data();
	detail::mem_matrix_store &res_store2
		= (detail::mem_matrix_store &) res->get_data();
#pragma omp parallel for
	for (size_t i = 0; i < res_store2.get_num_rows(); i++)
		for (size_t j = 0; j < res_store2.get_num_cols(); j++)
			assert(res_store2.get<int>(i, j)
					== orig_store2.get<int>(i, j) * vals->get<int>(j));
}

int main(int argc, char *argv[])
{
	int num_nodes = 1;
	int num_threads = 8;
	if (argc >= 3) {
		num_nodes = atoi(argv[1]);
		num_threads = atoi(argv[2]);
	}
	detail::mem_thread_pool::init_global_mem_threads(num_nodes,
			num_threads / num_nodes);

	test_scale_cols(-1);
	test_scale_cols(num_nodes);
	test_multiply_scalar(-1);
	test_multiply_scalar(num_nodes);
	test_ele_wise(-1);
	test_ele_wise(num_nodes);
	test_multiply_col(-1);
	test_multiply_col(num_nodes);
	test_agg_col(-1);
	test_agg_col(num_nodes);
	test_multiply_wide_row(-1);
	test_multiply_wide_row(num_nodes);
	test_multiply_tall_row(-1);
	test_multiply_tall_row(num_nodes);
	test_agg_row(-1);
	test_agg_row(num_nodes);
	test_agg_sub_col(-1);
	test_agg_sub_row(-1);
	test_rand_init();
#if 0
	test_conv_row_col();
	test_flatten();
	test_apply();
#endif
}
