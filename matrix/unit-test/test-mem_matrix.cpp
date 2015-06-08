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
	m1.materialize_self();
	m2.materialize_self();
	detail::mem_matrix_store::ptr res_store = detail::mem_matrix_store::create(
			m1.get_num_rows(), m2.get_num_cols(), matrix_layout_t::L_ROW,
			get_scalar_type<int>(), -1);
#pragma omp parallel for
	for (size_t i = 0; i < m1.get_num_rows(); i++) {
		for (size_t j = 0; j < m2.get_num_cols(); j++) {
			int sum = 0;
			for (size_t k = 0; k < m1.get_num_cols(); k++) {
				sum += m1.get<int>(i, k) * m2.get<int>(k, j);
			}
			res_store->set<int>(i, j, sum);
		}
	}
	return mem_dense_matrix::create(res_store);
}

void verify_result(const mem_dense_matrix &m1, const mem_dense_matrix &m2)
{
	assert(m1.get_num_rows() == m2.get_num_rows());
	assert(m1.get_num_cols() == m2.get_num_cols());

	m1.materialize_self();
	m2.materialize_self();
#pragma omp parallel for
	for (size_t i = 0; i < m1.get_num_rows(); i++)
		for (size_t j = 0; j < m1.get_num_cols(); j++)
			assert(m1.get<int>(i, j) == m2.get<int>(i, j));
}

enum matrix_val_t
{
	DEFAULT,
	SEQ,
	NUM_TYPES,
} matrix_val;

mem_dense_matrix::ptr create_seq_matrix(size_t nrow, size_t ncol,
		matrix_layout_t layout, int num_nodes, const scalar_type &type)
{
	if (type == get_scalar_type<int>()) {
		if (layout == matrix_layout_t::L_COL)
			return mem_dense_matrix::create(nrow, ncol, layout,
					type, set_col_operate(ncol), num_nodes);
		else
			return mem_dense_matrix::create(nrow, ncol, layout,
					type, set_row_operate(ncol), num_nodes);
	}
	else if (type == get_scalar_type<size_t>()) {
		if (layout == matrix_layout_t::L_COL)
			return mem_dense_matrix::create(nrow, ncol, layout,
					type, set_col_long_operate(ncol), num_nodes);
		else
			return mem_dense_matrix::create(nrow, ncol, layout,
					type, set_row_long_operate(ncol), num_nodes);
	}
	else
		return mem_dense_matrix::ptr();
}

mem_dense_matrix::ptr create_matrix(size_t nrow, size_t ncol,
		matrix_layout_t layout, int num_nodes,
		const scalar_type &type = get_scalar_type<int>())
{
	switch (matrix_val) {
		case matrix_val_t::DEFAULT:
			if (layout == matrix_layout_t::L_COL)
				return mem_dense_matrix::create(nrow, ncol, layout,
						type, num_nodes);
			else
				return mem_dense_matrix::create(nrow, ncol, layout,
						type, num_nodes);
		case matrix_val_t::SEQ:
			return create_seq_matrix(nrow, ncol, layout, num_nodes, type);
		default:
			assert(0);
			return mem_dense_matrix::ptr();
	}
}

void test_multiply_scalar(int num_nodes)
{
	printf("Test scalar multiplication\n");
	mem_dense_matrix::ptr orig = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes);
	mem_dense_matrix::ptr res = mem_dense_matrix::cast(orig->multiply_scalar(10));
	detail::mem_matrix_store &orig_store1
		= (detail::mem_matrix_store &) orig->get_data();
	detail::mem_matrix_store &res_store1
		= (detail::mem_matrix_store &) res->get_data();
	res_store1.materialize_self();
	orig_store1.materialize_self();
#pragma omp parallel for
	for (size_t i = 0; i < res_store1.get_num_rows(); i++)
		for (size_t j = 0; j < res_store1.get_num_cols(); j++)
			assert(res_store1.get<int>(i, j) == orig_store1.get<int>(i, j) * 10);

	orig = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes);
	res = mem_dense_matrix::cast(orig->multiply_scalar(10));
	detail::mem_matrix_store &orig_store2
		= (detail::mem_matrix_store &) orig->get_data();
	detail::mem_matrix_store &res_store2
		= (detail::mem_matrix_store &) res->get_data();
	res_store2.materialize_self();
	orig_store2.materialize_self();
#pragma omp parallel for
	for (size_t i = 0; i < res_store2.get_num_rows(); i++)
		for (size_t j = 0; j < res_store2.get_num_cols(); j++)
			assert(res_store2.get<int>(i, j) == orig_store2.get<int>(i, j) * 10);
}

void test_ele_wise(int num_nodes)
{
	printf("Test element-wise operations\n");
	mem_dense_matrix::ptr m1 = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes);
	mem_dense_matrix::ptr m2 = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes);
	mem_dense_matrix::ptr res = mem_dense_matrix::cast(m1->add(*m2));
	detail::mem_matrix_store &res_store = (detail::mem_matrix_store &) res->get_data();
	detail::mem_matrix_store &m1_store
		= (detail::mem_matrix_store &) m1->get_data();
	res_store.materialize_self();
	m1_store.materialize_self();
#pragma omp parallel for
	for (size_t i = 0; i < res_store.get_num_rows(); i++)
		for (size_t j = 0; j < res_store.get_num_cols(); j++)
			assert(res_store.get<int>(i, j) == m1_store.get<int>(i, j) * 2);
}

void test_multiply_col(int num_nodes)
{
	printf("Test multiplication on tall matrix stored column wise\n");
	mem_dense_matrix::ptr m1 = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes);
	mem_dense_matrix::ptr m2 = create_matrix(10, 9,
			matrix_layout_t::L_COL, num_nodes);
	mem_dense_matrix::ptr correct = naive_multiply(*m1, *m2);

	printf("Test multiply on col_matrix\n");
	mem_dense_matrix::ptr res1 = mem_dense_matrix::cast(m1->multiply(*m2));
	verify_result(*res1, *correct);
}

void test_agg_col(int num_nodes)
{
	printf("Test aggregation on tall matrix stored column wise\n");
	mem_dense_matrix::ptr m1 = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes, get_scalar_type<size_t>());
	const bulk_operate &op
		= m1->get_type().get_basic_ops().get_add();
	scalar_variable::ptr res = m1->aggregate(op);
	assert(res->get_type() == m1->get_type());
	assert(res->get_type() == get_scalar_type<size_t>());
	size_t sum = *(size_t *) res->get_raw();
	size_t num_eles = m1->get_num_rows() * m1->get_num_cols();
	if (matrix_val == matrix_val_t::DEFAULT)
		assert(sum == 0);
	else
		assert(sum == (num_eles - 1) * num_eles / 2);
}

void test_multiply_matrix(int num_nodes)
{
	mem_dense_matrix::ptr m1, m2, correct, res;

	printf("Test multiplication on wide row matrix X tall column matrix\n");
	m1 = create_matrix(10, long_dim, matrix_layout_t::L_ROW, num_nodes);
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_COL, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = mem_dense_matrix::cast(m1->multiply(*m2));
	verify_result(*res, *correct);

	printf("Test multiplication on wide row matrix X tall row matrix\n");
	m1 = create_matrix(10, long_dim, matrix_layout_t::L_ROW, num_nodes);
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_ROW, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = mem_dense_matrix::cast(m1->multiply(*m2));
	verify_result(*res, *correct);

	printf("Test multiplication on wide column matrix X tall column matrix\n");
	m1 = create_matrix(10, long_dim, matrix_layout_t::L_COL, num_nodes);
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_COL, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = mem_dense_matrix::cast(m1->multiply(*m2));
	verify_result(*res, *correct);

	printf("Test multiplication on wide column matrix X tall row matrix\n");
	m1 = create_matrix(10, long_dim, matrix_layout_t::L_COL, num_nodes);
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_ROW, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = mem_dense_matrix::cast(m1->multiply(*m2));
	verify_result(*res, *correct);

	printf("Test multiplication on tall row matrix X small row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes);
	m2 = create_matrix(10, 9, matrix_layout_t::L_ROW, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = mem_dense_matrix::cast(m1->multiply(*m2));
	verify_result(*res, *correct);

	printf("Test multiplication on tall row matrix X small column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes);
	m2 = create_matrix(10, 9, matrix_layout_t::L_COL, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = mem_dense_matrix::cast(m1->multiply(*m2));
	verify_result(*res, *correct);

	printf("Test multiplication on tall column matrix X small row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes);
	m2 = create_matrix(10, 9, matrix_layout_t::L_ROW, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = mem_dense_matrix::cast(m1->multiply(*m2));
	verify_result(*res, *correct);

	printf("Test multiplication on tall column matrix X small column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes);
	m2 = create_matrix(10, 9, matrix_layout_t::L_COL, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = mem_dense_matrix::cast(m1->multiply(*m2));
	verify_result(*res, *correct);
}

void test_agg_row(int num_nodes)
{
	printf("Test aggregation on tall matrix stored row wise\n");
	mem_dense_matrix::ptr m1 = create_matrix(long_dim, 10,
			matrix_layout_t::L_ROW, num_nodes, get_scalar_type<size_t>());
	const bulk_operate &op
		= m1->get_type().get_basic_ops().get_add();
	scalar_variable::ptr res = m1->aggregate(op);
	assert(res->get_type() == m1->get_type());
	size_t sum = *(size_t *) res->get_raw();
	size_t num_eles = m1->get_num_rows() * m1->get_num_cols();
	if (matrix_val == matrix_val_t::DEFAULT)
		assert(sum == 0);
	else
		assert(sum == (num_eles - 1) * num_eles / 2);
}

void test_agg_sub_col(int num_nodes)
{
	printf("Test aggregation on a column-wise submatrix\n");
	mem_dense_matrix::ptr col_m = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes, get_scalar_type<size_t>());
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
	if (matrix_val == matrix_val_t::DEFAULT)
		assert(sum == 0);
	else
		assert(sum == expected);
}

void test_agg_sub_row(int num_nodes)
{
	printf("Test aggregation on a row-wise submatrix\n");
	mem_dense_matrix::ptr col_m = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes, get_scalar_type<size_t>());
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

#endif

void test_scale_cols1(mem_dense_matrix::ptr orig)
{
	mem_vector::ptr vals = mem_vector::cast(
			create_vector<int>(0, orig->get_num_cols() - 1, 1));
	mem_dense_matrix::ptr res = mem_dense_matrix::cast(orig->scale_cols(vals));
	detail::mem_matrix_store &orig_store1
		= (detail::mem_matrix_store &) orig->get_data();
	detail::mem_matrix_store &res_store1
		= (detail::mem_matrix_store &) res->get_data();
	res_store1.materialize_self();
	orig_store1.materialize_self();
#pragma omp parallel for
	for (size_t i = 0; i < res_store1.get_num_rows(); i++)
		for (size_t j = 0; j < res_store1.get_num_cols(); j++)
			assert(res_store1.get<int>(i, j)
					== orig_store1.get<int>(i, j) * vals->get<int>(j));
}

void test_scale_cols(int num_nodes)
{
	printf("Test scale cols of tall column matrix\n");
	mem_dense_matrix::ptr orig = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes);
	test_scale_cols1(orig);

	printf("Test scale cols of tall row matrix\n");
	orig = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes);
	test_scale_cols1(orig);

	printf("Test scale cols of wide column matrix\n");
	orig = create_matrix(10, long_dim, matrix_layout_t::L_COL, num_nodes);
	test_scale_cols1(orig);

	printf("Test scale cols of wide row matrix\n");
	orig = create_matrix(10, long_dim, matrix_layout_t::L_ROW, num_nodes);
	test_scale_cols1(orig);
}

void test_scale_rows1(mem_dense_matrix::ptr orig)
{
	mem_vector::ptr vals = mem_vector::cast(
			create_vector<int>(0, orig->get_num_rows() - 1, 1));
	mem_dense_matrix::ptr res = mem_dense_matrix::cast(orig->scale_rows(vals));
	detail::mem_matrix_store &orig_store1
		= (detail::mem_matrix_store &) orig->get_data();
	detail::mem_matrix_store &res_store1
		= (detail::mem_matrix_store &) res->get_data();
	res_store1.materialize_self();
	orig_store1.materialize_self();
#pragma omp parallel for
	for (size_t i = 0; i < res_store1.get_num_rows(); i++)
		for (size_t j = 0; j < res_store1.get_num_cols(); j++)
			assert(res_store1.get<int>(i, j)
					== orig_store1.get<int>(i, j) * vals->get<int>(i));
}

void test_scale_rows(int num_nodes)
{
	printf("Test scale rows of wide row matrix\n");
	mem_dense_matrix::ptr orig = create_matrix(10, long_dim,
			matrix_layout_t::L_ROW, num_nodes);
	test_scale_rows1(orig);

	printf("Test scale rows of wide column matrix\n");
	orig = create_matrix(10, long_dim, matrix_layout_t::L_COL, num_nodes);
	test_scale_rows1(orig);

	printf("Test scale rows of tall row matrix\n");
	orig = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes);
	test_scale_rows1(orig);

	printf("Test scale rows of tall column matrix\n");
	orig = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes);
	test_scale_rows1(orig);
}

void test_create_const()
{
	printf("test create const matrix\n");
	mem_dense_matrix::ptr mat = mem_dense_matrix::create(10000, 10,
			matrix_layout_t::L_COL, get_scalar_type<int>());
	for (size_t i = 0; i < mat->get_num_rows(); i++)
		for (size_t j = 0; j < mat->get_num_cols(); j++)
			assert(mat->get<int>(i, j) == 0);

	mat = mem_dense_matrix::create_const<int>(1, 10000, 10,
			matrix_layout_t::L_COL);
	for (size_t i = 0; i < mat->get_num_rows(); i++)
		for (size_t j = 0; j < mat->get_num_cols(); j++)
			assert(mat->get<int>(i, j) == 1);
}

class sum_apply_op: public arr_apply_operate
{
public:
	sum_apply_op(): arr_apply_operate(1) {
	}

	void run(const local_vec_store &in, local_vec_store &out) const {
		assert(in.get_type() == get_scalar_type<int>());
		assert(out.get_type() == get_scalar_type<long>());
		long res = 0;
		for (size_t i = 0; i < in.get_length(); i++)
			res += in.get<int>(i);
		out.set<long>(0, res);
	}

	const scalar_type &get_input_type() const {
		return get_scalar_type<int>();
	}

	const scalar_type &get_output_type() const {
		return get_scalar_type<long>();
	}
};

void test_apply1(mem_dense_matrix::ptr mat)
{
	size_t num_rows = mat->get_num_rows();
	size_t num_cols = mat->get_num_cols();
	dense_matrix::ptr res;
	mem_vector::ptr res_vec;

	printf("Test apply on rows of a %s %s-wise matrix\n",
			mat->is_wide() ? "wide" : "tall",
			mat->store_layout() == matrix_layout_t::L_ROW ? "row" : "column");
	res = mat->apply(apply_margin::MAR_ROW,
			arr_apply_operate::const_ptr(new sum_apply_op()));
	assert(res->get_num_cols() == 1 && res->get_num_rows() == mat->get_num_rows());
	assert(res->is_type<long>());
	res_vec = mem_vector::cast(res->get_col(0));
	for (size_t i = 0; i < res_vec->get_length(); i++)
		assert(res_vec->get<long>(i)
				== i * num_cols * num_cols + (num_cols - 1) * num_cols / 2);

	printf("Test apply on columns of a %s %s-wise matrix\n",
			mat->is_wide() ? "wide" : "tall",
			mat->store_layout() == matrix_layout_t::L_ROW ? "row" : "column");
	res = mat->apply(apply_margin::MAR_COL,
			arr_apply_operate::const_ptr(new sum_apply_op()));
	assert(res->get_num_rows() == 1 && res->get_num_cols() == mat->get_num_cols());
	assert(res->is_type<long>());
	res_vec = mem_vector::cast(res->get_row(0));
	for (size_t i = 0; i < res_vec->get_length(); i++)
		assert(res_vec->get<long>(i)
				== (num_rows - 1) * num_rows / 2 * num_cols + num_rows * i);
}

void test_apply()
{
	detail::mem_matrix_store::ptr store;
	mem_dense_matrix::ptr mat;

	// Tall row-wise matrix
	store = detail::mem_matrix_store::create(long_dim, 10,
			matrix_layout_t::L_ROW, get_scalar_type<int>(), -1);
	store->set_data(set_row_operate(store->get_num_cols()));
	mat = mem_dense_matrix::create(store);
	test_apply1(mat);

	// Tall col-wise matrix
	store = detail::mem_matrix_store::create(long_dim, 10,
			matrix_layout_t::L_COL, get_scalar_type<int>(), -1);
	store->set_data(set_col_operate(store->get_num_cols()));
	mat = mem_dense_matrix::create(store);
	test_apply1(mat);

	// Wide row-wise matrix
	store = detail::mem_matrix_store::create(10, long_dim,
			matrix_layout_t::L_ROW, get_scalar_type<int>(), -1);
	store->set_data(set_row_operate(store->get_num_cols()));
	mat = mem_dense_matrix::create(store);
	test_apply1(mat);

	// wide col-wise matrix
	store = detail::mem_matrix_store::create(10, long_dim,
			matrix_layout_t::L_COL, get_scalar_type<int>(), -1);
	store->set_data(set_col_operate(store->get_num_cols()));
	mat = mem_dense_matrix::create(store);
	test_apply1(mat);
}

void test_conv_vec2mat()
{
	printf("convert a vector to a matrix\n");
	size_t len = 10000;
	size_t num_rows = 1000;
	vector::ptr vec = create_vector<int>(0, len, 1);
	dense_matrix::ptr mat = vec->conv2mat(num_rows, len / num_rows, true);
	assert(mat);
	assert(mat->store_layout() == matrix_layout_t::L_ROW);

	mat = vec->conv2mat(num_rows, len / num_rows, false);
	assert(mat);
	assert(mat->store_layout() == matrix_layout_t::L_COL);
}

void test_write2file1(detail::mem_matrix_store::ptr mat)
{
	char *tmp_file_name = tempnam(".", "tmp.mat");
	if (mat->store_layout() == matrix_layout_t::L_ROW)
		mat->set_data(set_row_operate(mat->get_num_cols()));
	else
		mat->set_data(set_col_operate(mat->get_num_cols()));
	bool ret = mat->write2file(tmp_file_name);
	assert(ret);

	detail::mem_matrix_store::ptr read_mat = detail::mem_matrix_store::load(tmp_file_name);
	assert(read_mat);
	assert(read_mat->get_num_rows() == mat->get_num_rows());
	assert(read_mat->get_num_cols() == mat->get_num_cols());
	assert(read_mat->get_type() == mat->get_type());
	assert(read_mat->store_layout() == mat->store_layout());
	for (size_t i = 0; i < mat->get_num_rows(); i++) {
		for (size_t j = 0; j < mat->get_num_cols(); j++)
			assert(mat->get<long>(i, j) == read_mat->get<long>(i, j));
	}

	unlink(tmp_file_name);
}

void test_write2file()
{
	detail::mem_matrix_store::ptr mat;
	printf("write a tall row matrix\n");
	mat = detail::mem_matrix_store::create(1000000, 10,
			matrix_layout_t::L_ROW, get_scalar_type<int>(), -1);
	test_write2file1(mat);
	printf("write a tall column matrix\n");
	mat = detail::mem_matrix_store::create(1000000, 10,
			matrix_layout_t::L_COL, get_scalar_type<int>(), -1);
	test_write2file1(mat);
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

	test_write2file();
	test_create_const();
	test_apply();
	test_conv_vec2mat();

	for (int i = 0; i < matrix_val_t::NUM_TYPES; i++) {
		matrix_val = (matrix_val_t) i;
		printf("matrix val type: %d\n", i);

		test_scale_cols(-1);
		test_scale_cols(num_nodes);
		test_scale_rows(-1);
		test_scale_rows(num_nodes);
		test_multiply_scalar(-1);
		test_multiply_scalar(num_nodes);
		test_ele_wise(-1);
		test_ele_wise(num_nodes);
		test_multiply_col(-1);
		test_multiply_col(num_nodes);
		test_agg_col(-1);
		test_agg_col(num_nodes);
		test_multiply_matrix(-1);
		test_multiply_matrix(num_nodes);
		test_agg_row(-1);
		test_agg_row(num_nodes);
		test_agg_sub_col(-1);
		test_agg_sub_row(-1);
		test_rand_init();
#if 0
		test_conv_row_col();
		test_flatten();
#endif
	}
}
