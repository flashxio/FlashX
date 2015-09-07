#include <stdio.h>
#include <cblas.h>

#include "vector.h"
#include "mem_worker_thread.h"
#include "dense_matrix.h"
#include "mem_matrix_store.h"
#include "local_matrix_store.h"
#include "sparse_matrix.h"
#include "EM_dense_matrix.h"
#include "matrix_stats.h"

#include "eigensolver/block_dense_matrix.h"
#include "eigensolver/collected_col_matrix_store.h"

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

size_t long_dim = 9999999;

/*
 * This is a naive implementation of matrix multiplication.
 * It should be correct
 */
dense_matrix::ptr naive_multiply(const dense_matrix &m1, const dense_matrix &m2)
{
	m1.materialize_self();
	m2.materialize_self();
	detail::mem_matrix_store::ptr res_store = detail::mem_matrix_store::create(
			m1.get_num_rows(), m2.get_num_cols(), matrix_layout_t::L_ROW,
			get_scalar_type<int>(), -1);
	detail::mem_matrix_store::const_ptr mem_m1;
	detail::mem_matrix_store::const_ptr mem_m2;
	if (m1.is_in_mem())
		mem_m1 = detail::mem_matrix_store::cast(m1.get_raw_store());
	else {
		dense_matrix::ptr mem_mat = m1.conv_store(true, -1);
		mem_m1 = detail::mem_matrix_store::cast(mem_mat->get_raw_store());
	}
	if (m2.is_in_mem())
		mem_m2 = detail::mem_matrix_store::cast(m2.get_raw_store());
	else {
		dense_matrix::ptr mem_mat = m2.conv_store(true, -1);
		mem_m2 = detail::mem_matrix_store::cast(mem_mat->get_raw_store());
	}
#pragma omp parallel for
	for (size_t i = 0; i < m1.get_num_rows(); i++) {
		for (size_t j = 0; j < m2.get_num_cols(); j++) {
			int sum = 0;
			for (size_t k = 0; k < m1.get_num_cols(); k++) {
				sum += mem_m1->get<int>(i, k) * mem_m2->get<int>(k, j);
			}
			res_store->set<int>(i, j, sum);
		}
	}
	return dense_matrix::create(res_store);
}

dense_matrix::ptr blas_multiply(const dense_matrix &m1, const dense_matrix &m2)
{
	dense_matrix::ptr tmp1 = m1.conv2(matrix_layout_t::L_COL);
	dense_matrix::ptr tmp2 = m2.conv2(matrix_layout_t::L_COL);
	tmp1->materialize_self();
	tmp2->materialize_self();
	detail::mem_col_matrix_store::ptr col_res = detail::mem_col_matrix_store::create(
			tmp1->get_num_rows(), tmp2->get_num_cols(), get_scalar_type<double>());
	detail::mem_matrix_store::const_ptr mem_m1;
	detail::mem_matrix_store::const_ptr mem_m2;
	if (tmp1->is_in_mem())
		mem_m1 = detail::mem_matrix_store::cast(tmp1->get_raw_store());
	else {
		dense_matrix::ptr mem_tmp = tmp1->conv_store(true, -1);
		mem_m1 = detail::mem_matrix_store::cast(mem_tmp->get_raw_store());
	}
	if (tmp2->is_in_mem())
		mem_m2 = detail::mem_matrix_store::cast(tmp2->get_raw_store());
	else {
		dense_matrix::ptr mem_tmp = tmp2->conv_store(true, -1);
		mem_m2 = detail::mem_matrix_store::cast(mem_tmp->get_raw_store());
	}
	assert(mem_m1->get_num_nodes() < 0);
	assert(mem_m2->get_num_nodes() < 0);

	detail::mem_col_matrix_store::const_ptr col_m1
		= detail::mem_col_matrix_store::cast(mem_m1);
	detail::mem_col_matrix_store::const_ptr col_m2
		= detail::mem_col_matrix_store::cast(mem_m2);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			mem_m1->get_num_rows(), mem_m2->get_num_cols(),
			mem_m1->get_num_cols(), 1,
			(const double *) col_m1->get_data().get_raw(),
			mem_m1->get_num_rows(),
			(const double *) col_m2->get_data().get_raw(),
			mem_m2->get_num_rows(), 0,
			(double *) col_res->get_data().get_raw(),
			col_res->get_num_rows());
	return dense_matrix::create(col_res);
}

struct approx_equal_func
{
public:
	bool operator()(const char *raw1, const char *raw2) const {
		double v1 = *(const double *) raw1;
		double v2 = *(const double *) raw2;
		double diff = v1 - v2;
		if (diff == 0)
			return true;

		if (v1 < 0)
			v1 = -v1;
		if (v2 < 0)
			v2 = -v2;
		if (diff < 0)
			diff = -diff;
		diff /= std::min(v1, v2);
		return diff < 1e-13;
	}
};

template<class T>
struct equal_func
{
public:
	bool operator()(const char *v1, const char *v2) const {
		return *(const T *) v1 == *(const T *) v2;
	}
};

template<class T1, class T2>
struct equal_func2
{
public:
	bool operator()(const char *v1, const char *v2) const {
		return *(const T1 *) v1 == *(const T2 *) v2;
	}
};

template<class T>
struct scale_equal_func
{
	T scale1;
	T scale2;
public:
	scale_equal_func(T scale1, T scale2) {
		this->scale1 = scale1;
		this->scale2 = scale2;
	}

	bool operator()(const char *v1, const char *v2) const {
		return (*(const T *) v1) * scale1 == (*(const T *) v2) * scale2;
	}
};

template<class Func>
void verify_result(const dense_matrix &m1, const dense_matrix &m2,
		const Func &func)
{
	assert(m1.get_num_rows() == m2.get_num_rows());
	assert(m1.get_num_cols() == m2.get_num_cols());

	m1.materialize_self();
	m2.materialize_self();

	detail::mem_matrix_store::const_ptr mem_m1;
	detail::mem_matrix_store::const_ptr mem_m2;
	if (m1.is_in_mem())
		mem_m1 = detail::mem_matrix_store::cast(m1.get_raw_store());
	else {
		dense_matrix::ptr mem_mat = m1.conv_store(true, -1);
		mem_m1 = detail::mem_matrix_store::cast(mem_mat->get_raw_store());
	}
	if (m2.is_in_mem())
		mem_m2 = detail::mem_matrix_store::cast(m2.get_raw_store());
	else {
		dense_matrix::ptr mem_mat = m2.conv_store(true, -1);
		mem_m2 = detail::mem_matrix_store::cast(mem_mat->get_raw_store());
	}

#pragma omp parallel for
	for (size_t i = 0; i < m1.get_num_rows(); i++)
		for (size_t j = 0; j < m1.get_num_cols(); j++)
			assert(func(mem_m1->get(i, j), mem_m2->get(i, j)));
}

enum matrix_val_t
{
	SEQ,
	DEFAULT,
	NUM_TYPES,
} matrix_val = matrix_val_t::SEQ;

dense_matrix::ptr create_seq_matrix(size_t nrow, size_t ncol,
		matrix_layout_t layout, int num_nodes, const scalar_type &type,
		bool in_mem)
{
	if (type == get_scalar_type<int>()) {
		if (layout == matrix_layout_t::L_COL)
			return dense_matrix::create(nrow, ncol, layout,
					type, set_col_operate(ncol), num_nodes, in_mem);
		else
			return dense_matrix::create(nrow, ncol, layout,
					type, set_row_operate(ncol), num_nodes, in_mem);
	}
	else if (type == get_scalar_type<size_t>()) {
		if (layout == matrix_layout_t::L_COL)
			return dense_matrix::create(nrow, ncol, layout,
					type, set_col_long_operate(ncol), num_nodes, in_mem);
		else
			return dense_matrix::create(nrow, ncol, layout,
					type, set_row_long_operate(ncol), num_nodes, in_mem);
	}
	else if (type == get_scalar_type<double>()) {
		if (layout == matrix_layout_t::L_COL)
			return dense_matrix::create(nrow, ncol, layout,
					get_scalar_type<size_t>(), set_col_long_operate(ncol),
					num_nodes, in_mem)->cast_ele_type(get_scalar_type<double>());
		else
			return dense_matrix::create(nrow, ncol, layout,
					get_scalar_type<size_t>(), set_row_long_operate(ncol),
					num_nodes, in_mem)->cast_ele_type(get_scalar_type<double>());
	}
	else
		return dense_matrix::ptr();
}

bool in_mem = true;

dense_matrix::ptr create_matrix(size_t nrow, size_t ncol,
		matrix_layout_t layout, int num_nodes,
		const scalar_type &type = get_scalar_type<int>())
{
	switch (matrix_val) {
		case matrix_val_t::DEFAULT:
			if (layout == matrix_layout_t::L_COL)
				return dense_matrix::create(nrow, ncol, layout,
						type, num_nodes, in_mem);
			else
				return dense_matrix::create(nrow, ncol, layout,
						type, num_nodes, in_mem);
		case matrix_val_t::SEQ:
			return create_seq_matrix(nrow, ncol, layout, num_nodes, type,
					in_mem);
		default:
			assert(0);
			return dense_matrix::ptr();
	}
}

void test_multiply_scalar(int num_nodes)
{
	printf("Test scalar multiplication\n");
	dense_matrix::ptr orig = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes);
	dense_matrix::ptr res = orig->multiply_scalar(10);
	verify_result(*res, *orig, scale_equal_func<int>(1, 10));

	orig = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes);
	res = orig->multiply_scalar(10);
	assert(res->is_virtual());
	assert(res->is_in_mem() == orig->is_in_mem());
	verify_result(*res, *orig, scale_equal_func<int>(1, 10));
}

void test_ele_wise(int num_nodes)
{
	printf("Test element-wise operations\n");
	dense_matrix::ptr m1 = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes);
	dense_matrix::ptr m2 = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes);
	dense_matrix::ptr res = m1->add(*m2);
	assert(res->is_virtual());
	assert(res->is_in_mem() == m1->is_in_mem());
	verify_result(*res, *m1, scale_equal_func<int>(1, 2));
}

void test_multiply_col(int num_nodes)
{
	printf("Test multiplication on tall matrix stored column wise\n");
	dense_matrix::ptr m1 = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes);
	dense_matrix::ptr m2 = create_matrix(10, 9,
			matrix_layout_t::L_COL, num_nodes);
	dense_matrix::ptr correct = naive_multiply(*m1, *m2);

	printf("Test multiply on col_matrix\n");
	dense_matrix::ptr res1 = m1->multiply(*m2);
	assert(res1->is_virtual());
	assert(res1->is_in_mem() == m1->is_in_mem());
	verify_result(*res1, *correct, equal_func<int>());
}

void test_agg_col(int num_nodes)
{
	printf("Test aggregation on tall matrix stored column wise\n");
	dense_matrix::ptr m1 = create_matrix(long_dim, 10,
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

void test_multiply_double(int num_nodes)
{
	dense_matrix::ptr m1, m2, correct, res;

	printf("Test multiplication on tall row matrix X small row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<double>());
	m2 = create_matrix(10, 9, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<double>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE, true);
	assert(res->store_layout() == matrix_layout_t::L_ROW);
	assert(res->is_virtual());
	assert(res->is_in_mem() == m1->is_in_mem());
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func());

	printf("Test multiplication on tall row matrix X small row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<double>());
	m2 = create_matrix(10, 9, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<double>());
	res = m1->multiply(*m2, matrix_layout_t::L_COL, true);
	assert(res->store_layout() == matrix_layout_t::L_COL);
	assert(res->is_virtual());
	assert(res->is_in_mem() == m1->is_in_mem());
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func());

	printf("Test multiplication on tall row matrix X small column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<double>());
	m2 = create_matrix(10, 9, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<double>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE, true);
	assert(res->store_layout() == matrix_layout_t::L_ROW);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func());

	printf("Test multiplication on tall column matrix X small row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<double>());
	m2 = create_matrix(10, 9, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<double>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE, true);
	assert(res->store_layout() == matrix_layout_t::L_COL);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func());

	printf("Test multiplication on tall column matrix X small row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<double>());
	m2 = create_matrix(10, 9, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<double>());
	res = m1->multiply(*m2, matrix_layout_t::L_ROW, true);
	assert(res->store_layout() == matrix_layout_t::L_ROW);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func());

	printf("Test multiplication on tall column matrix X small column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<double>());
	m2 = create_matrix(10, 9, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<double>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE, true);
	assert(res->store_layout() == matrix_layout_t::L_COL);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func());

	printf("Test multiplication on wide row matrix X tall column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<double>());
	m1 = m1->transpose();
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<double>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE, true);
	assert(res->store_layout() == matrix_layout_t::L_ROW);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func());

	printf("Test multiplication on wide row matrix X tall column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<double>());
	m1 = m1->transpose();
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<double>());
	res = m1->multiply(*m2, matrix_layout_t::L_COL, true);
	assert(res->store_layout() == matrix_layout_t::L_COL);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func());

	printf("Test multiplication on wide row matrix X tall row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<double>());
	m1 = m1->transpose();
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<double>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE, true);
	assert(res->store_layout() == matrix_layout_t::L_ROW);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func());

	printf("Test multiplication on wide column matrix X tall column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<double>());
	m1 = m1->transpose();
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<double>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE, true);
	assert(res->store_layout() == matrix_layout_t::L_COL);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func());

	printf("Test multiplication on wide column matrix X tall column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<double>());
	m1 = m1->transpose();
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<double>());
	res = m1->multiply(*m2, matrix_layout_t::L_ROW, true);
	assert(res->store_layout() == matrix_layout_t::L_ROW);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func());

	printf("Test multiplication on wide column matrix X tall row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<double>());
	m1 = m1->transpose();
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<double>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE, true);
	assert(res->store_layout() == matrix_layout_t::L_COL);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func());
}

void test_multiply_matrix(int num_nodes)
{
	dense_matrix::ptr m1, m2, correct, res;

	printf("Test multiplication on wide row matrix X tall column matrix\n");
	m1 = create_matrix(10, long_dim, matrix_layout_t::L_ROW, num_nodes);
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_COL, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = m1->multiply(*m2);
	verify_result(*res, *correct, equal_func<int>());

	printf("Test multiplication on wide row matrix X tall row matrix\n");
	m1 = create_matrix(10, long_dim, matrix_layout_t::L_ROW, num_nodes);
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_ROW, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = m1->multiply(*m2);
	verify_result(*res, *correct, equal_func<int>());

	printf("Test multiplication on wide column matrix X tall column matrix\n");
	m1 = create_matrix(10, long_dim, matrix_layout_t::L_COL, num_nodes);
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_COL, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = m1->multiply(*m2);
	verify_result(*res, *correct, equal_func<int>());

	printf("Test multiplication on wide column matrix X tall row matrix\n");
	m1 = create_matrix(10, long_dim, matrix_layout_t::L_COL, num_nodes);
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_ROW, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = m1->multiply(*m2);
	verify_result(*res, *correct, equal_func<int>());

	printf("Test multiplication on tall row matrix X small row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes);
	m2 = create_matrix(10, 9, matrix_layout_t::L_ROW, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = m1->multiply(*m2);
	assert(res->is_virtual());
	assert(res->is_in_mem() == m1->is_in_mem());
	verify_result(*res, *correct, equal_func<int>());

	printf("Test multiplication on tall row matrix X small column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes);
	m2 = create_matrix(10, 9, matrix_layout_t::L_COL, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = m1->multiply(*m2);
	verify_result(*res, *correct, equal_func<int>());

	printf("Test multiplication on tall column matrix X small row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes);
	m2 = create_matrix(10, 9, matrix_layout_t::L_ROW, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = m1->multiply(*m2);
	verify_result(*res, *correct, equal_func<int>());

	printf("Test multiplication on tall column matrix X small column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes);
	m2 = create_matrix(10, 9, matrix_layout_t::L_COL, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = m1->multiply(*m2);
	verify_result(*res, *correct, equal_func<int>());
}

void test_agg_row(int num_nodes)
{
	printf("Test aggregation on tall matrix stored row wise\n");
	dense_matrix::ptr m1 = create_matrix(long_dim, 10,
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
	dense_matrix::ptr col_m = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes, get_scalar_type<size_t>());
	std::vector<off_t> idxs(4);
	idxs[0] = 1;
	idxs[1] = 5;
	idxs[2] = 3;
	idxs[3] = 4;
	dense_matrix::ptr sub_m = col_m->get_cols(idxs);
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
	dense_matrix::ptr col_m = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes, get_scalar_type<size_t>());
	std::vector<off_t> idxs(4);
	idxs[0] = 1;
	idxs[1] = 5;
	idxs[2] = 3;
	idxs[3] = 4;
	dense_matrix::ptr sub_col_m = col_m->get_cols(idxs);
	dense_matrix::ptr sub_row_m = sub_col_m->transpose();

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

void test_rand_init()
{
	printf("test rand init\n");
	dense_matrix::ptr m = dense_matrix::create_rand<double>(-1.0, 1.0,
			long_dim / 100, 10, matrix_layout_t::L_COL);
	double sum = 0;
	const detail::mem_matrix_store &mem_m
		= dynamic_cast<const detail::mem_matrix_store &>(*m);
	for (size_t i = 0; i < m->get_num_rows(); i++)
		for (size_t j = 0; j < m->get_num_cols(); j++) {
			double v = mem_m.get<double>(i, j);
			assert(v >= -1.0 && v <= 1.0);
			sum += v;
		}
	printf("sum: %f\n", sum);
}

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

void test_scale_cols1(dense_matrix::ptr orig)
{
	vector::ptr vals = create_vector<int>(0, orig->get_num_cols() - 1, 1);
	dense_matrix::ptr res = orig->scale_cols(vals);
	assert(res->is_virtual());
	assert(res->is_in_mem() == orig->is_in_mem());
	res->materialize_self();
	orig->materialize_self();
	if (res->is_in_mem()) {
		assert(orig->is_in_mem());
		const detail::mem_matrix_store &orig_store1
			= dynamic_cast<const detail::mem_matrix_store &>(orig->get_data());
		const detail::mem_matrix_store &res_store1
			= dynamic_cast<const detail::mem_matrix_store &>(res->get_data());
		const detail::smp_vec_store &val_store
			= dynamic_cast<const detail::smp_vec_store &>(vals->get_data());
#pragma omp parallel for
		for (size_t i = 0; i < res_store1.get_num_rows(); i++)
			for (size_t j = 0; j < res_store1.get_num_cols(); j++)
				assert(res_store1.get<int>(i, j)
						== orig_store1.get<int>(i, j) * val_store.get<int>(j));
	}
}

void test_scale_cols(int num_nodes)
{
	printf("Test scale cols of tall column matrix\n");
	dense_matrix::ptr orig = create_matrix(long_dim, 10,
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

void test_scale_rows1(dense_matrix::ptr orig)
{
	vector::ptr vals = create_vector<int>(0, orig->get_num_rows() - 1, 1);
	dense_matrix::ptr res = orig->scale_rows(vals);
	assert(res->is_virtual());
	assert(res->is_in_mem() == orig->is_in_mem());
	res->materialize_self();
	orig->materialize_self();
	if (res->is_in_mem()) {
		assert(orig->is_in_mem());
		const detail::mem_matrix_store &orig_store1
			= dynamic_cast<const detail::mem_matrix_store &>(orig->get_data());
		const detail::mem_matrix_store &res_store1
			= dynamic_cast<const detail::mem_matrix_store &>(res->get_data());
		const detail::smp_vec_store &val_store
			= dynamic_cast<const detail::smp_vec_store &>(vals->get_data());
#pragma omp parallel for
		for (size_t i = 0; i < res_store1.get_num_rows(); i++)
			for (size_t j = 0; j < res_store1.get_num_cols(); j++)
				assert(res_store1.get<int>(i, j)
						== orig_store1.get<int>(i, j) * val_store.get<int>(i));
	}
}

void test_scale_rows(int num_nodes)
{
	printf("Test scale rows of wide row matrix\n");
	dense_matrix::ptr orig = create_matrix(10, long_dim,
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

#if 0
void test_create_const()
{
	printf("test create const matrix\n");
	dense_matrix::ptr mat = dense_matrix::create(10000, 10,
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
#endif

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

void test_apply1(dense_matrix::ptr mat)
{
	size_t num_rows = mat->get_num_rows();
	size_t num_cols = mat->get_num_cols();
	dense_matrix::ptr res;
	vector::ptr res_vec;

	printf("Test apply on rows of a %s %s-wise matrix\n",
			mat->is_wide() ? "wide" : "tall",
			mat->store_layout() == matrix_layout_t::L_ROW ? "row" : "column");
	res = mat->apply(apply_margin::MAR_ROW,
			arr_apply_operate::const_ptr(new sum_apply_op()));
	if (res) {
		assert(res->get_num_cols() == 1 && res->get_num_rows() == mat->get_num_rows());
		assert(res->is_type<long>());
		res->materialize_self();
		res_vec = res->get_col(0);
		const detail::smp_vec_store &vstore1
			= dynamic_cast<const detail::smp_vec_store &>(res_vec->get_data());
		for (size_t i = 0; i < res_vec->get_length(); i++)
			assert(vstore1.get<long>(i)
					== i * num_cols * num_cols + (num_cols - 1) * num_cols / 2);
	}

	printf("Test apply on columns of a %s %s-wise matrix\n",
			mat->is_wide() ? "wide" : "tall",
			mat->store_layout() == matrix_layout_t::L_ROW ? "row" : "column");
	res = mat->apply(apply_margin::MAR_COL,
			arr_apply_operate::const_ptr(new sum_apply_op()));
	if (res) {
		assert(res->get_num_rows() == 1 && res->get_num_cols() == mat->get_num_cols());
		assert(res->is_type<long>());
		res->materialize_self();
		res_vec = res->get_row(0);
		const detail::smp_vec_store &vstore2
			= dynamic_cast<const detail::smp_vec_store &>(res_vec->get_data());
		for (size_t i = 0; i < res_vec->get_length(); i++)
			assert(vstore2.get<long>(i)
					== (num_rows - 1) * num_rows / 2 * num_cols + num_rows * i);
	}
}

void test_apply()
{
	dense_matrix::ptr mat;

	// Tall row-wise matrix
	mat = create_seq_matrix(long_dim, 10,
			matrix_layout_t::L_ROW, -1, get_scalar_type<int>(), in_mem);
	test_apply1(mat);

	// Tall col-wise matrix
	mat = create_seq_matrix(long_dim, 10,
			matrix_layout_t::L_COL, -1, get_scalar_type<int>(), in_mem);
	test_apply1(mat);

	// Wide row-wise matrix
	mat = create_seq_matrix(10, long_dim,
			matrix_layout_t::L_ROW, -1, get_scalar_type<int>(), in_mem);
	test_apply1(mat);

	// wide col-wise matrix
	mat = create_seq_matrix(10, long_dim,
			matrix_layout_t::L_COL, -1, get_scalar_type<int>(), in_mem);
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
			assert(mat->get<int>(i, j) == read_mat->get<int>(i, j));
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

void test_cast()
{
	printf("test cast type\n");
	dense_matrix::ptr mat, mat1;

	mat = dense_matrix::create_randu<int>(
			0, 1000, long_dim, 10, matrix_layout_t::L_ROW, -1, in_mem);
	mat1 = mat->cast_ele_type(get_scalar_type<long>());
	assert(mat1->is_virtual());
	assert(mat1->is_in_mem() == mat->is_in_mem());
	verify_result(*mat, *mat1, equal_func2<int, long>());

	mat = dense_matrix::create_randu<float>(
			0, 1000, long_dim, 10, matrix_layout_t::L_ROW, -1, in_mem);
	mat1 = mat->cast_ele_type(get_scalar_type<double>());
	verify_result(*mat, *mat1, equal_func2<float, double>());
}

void test_conv2(int num_nodes)
{
	// Test conv2
	printf("conv2 layout\n");
	dense_matrix::ptr mat = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes);
	dense_matrix::ptr mat1 = mat->conv2(matrix_layout_t::L_ROW);
	assert(mat1->is_virtual());
	assert(mat1->is_in_mem() == mat->is_in_mem());
	assert(mat->get_num_rows() == mat1->get_num_rows());
	assert(mat->get_num_cols() == mat1->get_num_cols());
	assert(mat1->store_layout() == matrix_layout_t::L_ROW);
	verify_result(*mat, *mat1, equal_func<int>());

	// Test transpose of conv2
	dense_matrix::ptr t_mat = mat->transpose();
	dense_matrix::ptr t_mat1 = mat1->transpose();
	assert(t_mat->get_num_rows() == t_mat1->get_num_rows());
	assert(t_mat->get_num_cols() == t_mat1->get_num_cols());
	assert(t_mat1->store_layout() == matrix_layout_t::L_COL);
	verify_result(*t_mat, *t_mat1, equal_func<int>());

	mat = mat1->conv2(matrix_layout_t::L_COL);
	assert(mat->get_num_rows() == mat1->get_num_rows());
	assert(mat->get_num_cols() == mat1->get_num_cols());
	verify_result(*mat, *mat1, equal_func<int>());
}

void test_sum_row_col1(dense_matrix::ptr mat)
{
	printf("test row sum on %s %s matrix\n", mat->is_wide() ? "wide" : "tall",
			mat->store_layout() == matrix_layout_t::L_COL ? "column" : "row");
	vector::ptr vec = mat->row_sum();
	assert(vec->is_in_mem());
	assert(vec->get_length() == mat->get_num_rows());

	detail::mem_vec_store::const_ptr mem_vec
		= detail::mem_vec_store::cast(vec->get_raw_store());
	detail::mem_matrix_store::const_ptr mem_m;
	if (mat->is_in_mem())
		mem_m = detail::mem_matrix_store::cast(mat->get_raw_store());
	else {
		dense_matrix::ptr mem_mat = mat->conv_store(true, -1);
		mem_m = detail::mem_matrix_store::cast(mem_mat->get_raw_store());
	}
	for (size_t i = 0; i < vec->get_length(); i++) {
		int sum = 0;
		for (size_t j = 0; j < mat->get_num_cols(); j++)
			sum += mem_m->get<int>(i, j);
		assert(sum == *(int *) mem_vec->get_sub_arr(i, i + 1));
	}

	printf("test col sum on %s %s matrix\n", mat->is_wide() ? "wide" : "tall",
			mat->store_layout() == matrix_layout_t::L_COL ? "column" : "row");
	vec = mat->col_sum();
	assert(vec->is_in_mem());
	assert(vec->get_length() == mat->get_num_cols());
	mem_vec = detail::mem_vec_store::cast(vec->get_raw_store());
	for (size_t i = 0; i < vec->get_length(); i++) {
		int sum = 0;
		for (size_t j = 0; j < mat->get_num_rows(); j++)
			sum += mem_m->get<int>(j, i);
		assert(sum == *(int *) mem_vec->get_sub_arr(i, i + 1));
	}
}

void test_sum_row_col(int num_nodes)
{
	dense_matrix::ptr mat = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes);
	test_sum_row_col1(mat);

	mat = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes);
	test_sum_row_col1(mat);

	mat = create_matrix(10, long_dim, matrix_layout_t::L_COL, num_nodes);
	test_sum_row_col1(mat);

	mat = create_matrix(10, long_dim, matrix_layout_t::L_ROW, num_nodes);
	test_sum_row_col1(mat);
}

void test_mapply_chain(int num_nodes, const scalar_type &type)
{
	printf("test mapply chain for %s\n",
			type == get_scalar_type<int>() ? "int" : "double");
	dense_matrix::ptr orig_mat1 = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes, type);
	dense_matrix::ptr orig_mat2 = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes, type);
	dense_matrix::ptr smat1 = create_seq_matrix(10, 10,
			matrix_layout_t::L_COL, num_nodes, type, true);
	dense_matrix::ptr smat2 = create_seq_matrix(10, 10,
			matrix_layout_t::L_COL, num_nodes, type, true);
	dense_matrix::ptr smat3 = create_seq_matrix(10, 10,
			matrix_layout_t::L_COL, num_nodes, type, true);
	dense_matrix::ptr smat4 = create_seq_matrix(10, 10,
			matrix_layout_t::L_COL, num_nodes, type, true);
	dense_matrix::ptr smat5 = create_seq_matrix(10, 10,
			matrix_layout_t::L_COL, num_nodes, type, true);

	printf("test a chain of mapply virtual matrices\n");
	detail::matrix_stats_t orig_stats = detail::matrix_stats;
	dense_matrix::ptr vmat1 = orig_mat1->multiply(*smat1,
			matrix_layout_t::L_NONE, true);
	dense_matrix::ptr vmat3 = orig_mat1->multiply(*smat5,
			matrix_layout_t::L_NONE, true);
	dense_matrix::ptr vmat4 = vmat3->multiply(*smat2,
			matrix_layout_t::L_NONE, true);
	dense_matrix::ptr vmat5 = orig_mat2->multiply(*smat3,
			matrix_layout_t::L_NONE, true);
	dense_matrix::ptr vmat6 = vmat4->multiply(*smat4,
			matrix_layout_t::L_NONE, true)->add(*vmat5);
	dense_matrix::ptr res = vmat6->transpose()->multiply(*vmat1,
			matrix_layout_t::L_NONE, true);
	detail::matrix_stats.print_diff(orig_stats);

	printf("materialize every matrix operations individually\n");
	detail::matrix_stats_t orig_stats1 = detail::matrix_stats;
	dense_matrix::ptr mat1 = orig_mat1->multiply(*smat1,
			matrix_layout_t::L_NONE, true);
	mat1->materialize_self();
	dense_matrix::ptr mat3 = orig_mat1->multiply(*smat5,
			matrix_layout_t::L_NONE, true);
	mat3->materialize_self();
	dense_matrix::ptr mat4 = vmat3->multiply(*smat2,
			matrix_layout_t::L_NONE, true);
	mat4->materialize_self();
	dense_matrix::ptr mat5 = orig_mat2->multiply(*smat3,
			matrix_layout_t::L_NONE, true);
	mat5->materialize_self();
	dense_matrix::ptr mat6 = vmat4->multiply(*smat4,
			matrix_layout_t::L_NONE, true)->add(*vmat5);
	mat6->materialize_self();
	dense_matrix::ptr res1 = mat6->transpose()->multiply(*mat1,
			matrix_layout_t::L_NONE, true);
	detail::matrix_stats.print_diff(orig_stats1);

	if (type == get_scalar_type<int>())
		verify_result(*res, *res1, equal_func<int>());
	else
		verify_result(*res, *res1, approx_equal_func());
}

class split_op: public detail::portion_mapply_op
{
public:
	split_op(): portion_mapply_op(0, 0, get_scalar_type<int>()) {
	}
	void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			const std::vector<detail::local_matrix_store::ptr> &outs) const;

	virtual portion_mapply_op::const_ptr transpose() const {
		return portion_mapply_op::const_ptr();
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return std::string();
	}
};

void split_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		const std::vector<detail::local_matrix_store::ptr> &outs) const
{
	assert(ins.size() == 1);
	assert(ins[0]->store_layout() == matrix_layout_t::L_COL);
	for (size_t i = 0; i < outs.size(); i++)
		assert(outs[i]->store_layout() == matrix_layout_t::L_COL);
	const detail::local_col_matrix_store &col_in
		= dynamic_cast<const detail::local_col_matrix_store &>(*ins[0]);
	size_t col_idx = 0;
	for (size_t j = 0; j < outs.size(); j++) {
		detail::local_col_matrix_store &col_out
			= dynamic_cast<detail::local_col_matrix_store &>(*outs[j]);
		for (size_t i = 0; i < col_out.get_num_cols(); i++) {
			assert(col_idx < col_in.get_num_cols());
			memcpy(col_out.get_col(i), col_in.get_col(col_idx),
					col_in.get_num_rows() * col_in.get_entry_size());
			col_idx++;
		}
	}
	assert(col_idx == ins[0]->get_num_cols());
}

void test_mul_output(int num_nodes)
{
	printf("test multiple output\n");
	dense_matrix::ptr mat = create_matrix(long_dim, 10, matrix_layout_t::L_COL,
			num_nodes, get_scalar_type<double>());
	std::vector<detail::matrix_store::const_ptr> in(1);
	in[0] = mat->get_raw_store();
	std::vector<detail::matrix_store::ptr> out(2);
	out[0] = detail::matrix_store::create(long_dim, 5, matrix_layout_t::L_COL,
			get_scalar_type<double>(), mat->get_data().get_num_nodes(),
			mat->is_in_mem());
	out[1] = detail::matrix_store::create(long_dim, 5, matrix_layout_t::L_COL,
			get_scalar_type<double>(), mat->get_data().get_num_nodes(),
			mat->is_in_mem());
	bool ret = __mapply_portion(in,
			detail::portion_mapply_op::const_ptr(new split_op()), out);
	assert(ret);

	dense_matrix::ptr out_mat0 = dense_matrix::create(out[0]);
	dense_matrix::ptr out_mat1 = dense_matrix::create(out[1]);
	dense_matrix::ptr in_mat = dense_matrix::create(in[0]);
	scalar_variable::ptr agg0 = out_mat0->aggregate(
			out[0]->get_type().get_basic_ops().get_add());
	scalar_variable::ptr agg1 = out_mat1->aggregate(
			out[1]->get_type().get_basic_ops().get_add());
	scalar_variable::ptr agg = in_mat->aggregate(
			out[0]->get_type().get_basic_ops().get_add());
	assert(*(double *) agg1->get_raw() - *(double *) agg0->get_raw()
			== out[0]->get_num_cols() * out[0]->get_num_cols() * out[0]->get_num_rows());
	assert(*(double *) agg0->get_raw() + *(double *) agg1->get_raw()
			== *(double *) agg->get_raw());
}

class set_col_operate1: public type_set_operate<int>
{
	size_t num_cols;
public:
	set_col_operate1(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(int *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = (row_idx + i) * num_cols + col_idx + 1;
		}
	}
};

void test_min()
{
	printf("test min\n");
	const bulk_operate &op = *get_scalar_type<int>().get_basic_ops().get_op(
			basic_ops::op_idx::MIN);

	dense_matrix::ptr mat = dense_matrix::create(2000, 1,
			matrix_layout_t::L_COL, get_scalar_type<int>(), set_col_operate1(4));
	scalar_variable::ptr res = mat->aggregate(op);
	assert(*(int *) res->get_raw() == 1);

	mat = dense_matrix::create(long_dim, 1, matrix_layout_t::L_COL,
			get_scalar_type<int>(), set_col_operate1(4));
	res = mat->aggregate(op);
	assert(*(int *) res->get_raw() == 1);
}

void test_apply_scalar()
{
	printf("test apply scalar\n");
	dense_matrix::ptr mat1 = dense_matrix::create_randu<double>(0, 1,
			long_dim, 10, matrix_layout_t::L_COL);
	dense_matrix::ptr mat2 = dense_matrix::create_randu<double>(0, 1,
			long_dim, 10, matrix_layout_t::L_COL);
	dense_matrix::ptr tmp = mat1->minus(*mat2);
	assert(tmp->get_type() == get_scalar_type<double>());
	tmp = tmp->abs();
	assert(tmp->get_type() == get_scalar_type<double>());
	tmp = tmp->lt_scalar<double>(1e-5);
	assert(tmp->get_type() == get_scalar_type<bool>());
	scalar_variable::ptr res = tmp->sum();
	assert(res->get_type() == get_scalar_type<size_t>());

	detail::mem_matrix_store::const_ptr store1
		= detail::mem_matrix_store::cast(mat1->get_raw_store());
	detail::mem_matrix_store::const_ptr store2
		= detail::mem_matrix_store::cast(mat2->get_raw_store());
	size_t num = 0;
	for (size_t i = 0; i < store1->get_num_rows(); i++)
		for (size_t j = 0; j < store2->get_num_cols(); j++) {
			if (std::abs(store1->get<double>(i, j) - store2->get<double>(i, j)) < 1e-5)
				num++;
		}
	assert(num == *(size_t *) res->get_raw());
	printf("two matrices have %ld elements with similar values\n", num);
}

class add_portion_op: public detail::portion_mapply_op
{
public:
	add_portion_op(size_t num_rows, size_t num_cols): detail::portion_mapply_op(
			num_rows, num_cols, get_scalar_type<double>()) {
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		out.reset_data();
		for (size_t i = 0; i < ins.size(); i++) {
			assert(ins[i]->store_layout() == matrix_layout_t::L_COL);
			detail::local_matrix_store::const_ptr in = ins[i]->get_portion(0, 0,
					out.get_num_rows(), out.get_num_cols());
			assert(in);
			detail::mapply2(out, *in,
					out.get_type().get_basic_ops().get_add(), out);
		}
	}

	virtual portion_mapply_op::const_ptr transpose() const {
		return portion_mapply_op::const_ptr();
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return std::string();
	}
};

void test_mapply_mixed(int num_nodes)
{
	printf("test serial and parallel mapply\n");
	std::vector<dense_matrix::ptr> mats(5);
	mats[0] = dense_matrix::create_randu<size_t>(0, 1000,
			long_dim, 10, matrix_layout_t::L_COL, -1, true);
	mats[1] = dense_matrix::create_randu<size_t>(0, 1000,
			long_dim, 11, matrix_layout_t::L_COL, -1, false);
	mats[2] = dense_matrix::create_randu<size_t>(0, 1000,
			long_dim, 12, matrix_layout_t::L_COL, -1, false);
	mats[3] = dense_matrix::create_randu<size_t>(0, 1000,
			long_dim, 13, matrix_layout_t::L_COL, num_nodes, true);
	mats[4] = dense_matrix::create_randu<size_t>(0, 1000,
			long_dim, 14, matrix_layout_t::L_COL, -1, false);

	detail::portion_mapply_op::const_ptr op(new add_portion_op(
				long_dim, 10));
	std::vector<detail::matrix_store::const_ptr> stores(5);
	for (size_t i = 0; i < stores.size(); i++)
		stores[i] = mats[i]->get_raw_store();
	dense_matrix::ptr par_res = dense_matrix::create(detail::__mapply_portion(
				stores, op, matrix_layout_t::L_COL, true));
	dense_matrix::ptr serial_res = dense_matrix::create(detail::__mapply_portion(
				stores, op, matrix_layout_t::L_COL, false));
	scalar_variable::ptr max_diff = par_res->minus(*serial_res)->abs()->max();
	printf("max diff: %ld\n", *(size_t *) max_diff->get_raw());
	assert(*(size_t *) max_diff->get_raw() == 0);
}

void test_sub_matrix()
{
	printf("test sub tall col-matrix\n");
	size_t long_dim = 16 * 1024 * 1024;
	dense_matrix::ptr mat = dense_matrix::create_randu<int>(0, 1000,
			long_dim, 10, matrix_layout_t::L_COL, -1, false);
	mat = mat->cast_ele_type(get_scalar_type<double>());
	mat->materialize_self();

	// Get the first sub-matrix.
	std::vector<off_t> sub_offs1(2);
	sub_offs1[0] = 1;
	sub_offs1[1] = 2;
	dense_matrix::ptr sub_mat1 = mat->get_cols(sub_offs1);
	// Get the in-mem copy of the first sub-matrix.
	detail::matrix_stats_t orig_stats1 = detail::matrix_stats;
	dense_matrix::ptr copy1 = sub_mat1->conv_store(true, -1);
	printf("copy the first sub matrix\n");
	detail::matrix_stats.print_diff(orig_stats1);

	{
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		dense_matrix::ptr diff = sub_mat1->minus(*copy1);
		scalar_variable::ptr max_var = diff->abs()->max();
		printf("max diff: %g\n", *(const double *) max_var->get_raw());
		detail::matrix_stats.print_diff(orig_stats);
		assert(*(const double *) max_var->get_raw() == 0);
	}

	// Get the second sub-matrix and its in-mem copy.
	std::vector<off_t> sub_offs2(2);
	sub_offs2[0] = 4;
	sub_offs2[1] = 5;
	dense_matrix::ptr sub_mat2 = mat->get_cols(sub_offs2);
	dense_matrix::ptr copy2 = sub_mat2->conv_store(true, -1);

	// Get the third matrix (in-mem).
	dense_matrix::ptr mat3 = dense_matrix::create_randu<int>(0, 1000,
			long_dim, 2, matrix_layout_t::L_COL, -1, true);
	mat3 = mat3->cast_ele_type(get_scalar_type<double>());
	mat3->materialize_self();

	// Create the block MV from the sub-matrices.
	size_t num_cols = sub_mat1->get_num_cols() + sub_mat2->get_num_cols()
		+ mat3->get_num_cols();
	eigen::block_multi_vector::ptr mv1 = eigen::block_multi_vector::create(
			long_dim, num_cols, 2, get_scalar_type<double>(), in_mem);
	mv1->set_block(0, sub_mat1);
	mv1->set_block(1, sub_mat2);
	mv1->set_block(2, mat3);

	// Create the block MV from the in-mem matrices.
	eigen::block_multi_vector::ptr mv2 = eigen::block_multi_vector::create(
			long_dim, num_cols, 2, get_scalar_type<double>(), in_mem);
	mv2->set_block(0, copy1);
	mv2->set_block(1, copy2);
	mv2->set_block(2, mat3);

	// Create the right matrix for GEMM.
	dense_matrix::ptr right_mat = dense_matrix::create_randu<int>(0, 1000,
			num_cols, 2, matrix_layout_t::L_COL, -1, true);
	detail::mem_col_matrix_store::const_ptr B
		= detail::mem_col_matrix_store::cast(right_mat->get_raw_store());
	scalar_variable_impl<double> alpha(2);
	scalar_variable_impl<double> beta(3);

	// Perform GEMM on the first block MV.
	eigen::block_multi_vector::ptr res1;
	dense_matrix::ptr res_mat1;
	{
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		res1 = eigen::block_multi_vector::create(
				long_dim, mv1->get_block_size(), mv1->get_block_size(),
				get_scalar_type<double>(), false);
		res1->set_block(0, create_seq_matrix(long_dim, mv1->get_block_size(),
					matrix_layout_t::L_COL, -1, get_scalar_type<double>(), false));
		res1 = res1->gemm(*mv1, B, alpha, beta);
		res_mat1 = res1->get_block(0);
		res_mat1->materialize_self();
		printf("The first GEMM on the submatrix\n");
		detail::matrix_stats.print_diff(orig_stats);
	}

	// Perform GEMM on the second block MV.
	eigen::block_multi_vector::ptr res2;
	dense_matrix::ptr res_mat2;
	{
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		res2 = eigen::block_multi_vector::create(
				long_dim, mv1->get_block_size(), mv1->get_block_size(),
				get_scalar_type<double>(), true);
		res2->set_block(0, create_seq_matrix(long_dim, mv1->get_block_size(),
					matrix_layout_t::L_COL, -1, get_scalar_type<double>(), true));
		res2 = res2->gemm(*mv2, B, alpha, beta);
		res_mat2 = res2->get_block(0);
		res_mat2->materialize_self();
		printf("The second GEMM on the in-mem matrix\n");
		detail::matrix_stats.print_diff(orig_stats);
	}

	// Compare the result of GEMM.
	dense_matrix::ptr diff = res_mat1->minus(*res_mat2);
	scalar_variable::ptr max_var = diff->abs()->max();
	printf("max diff: %g\n", *(const double *) max_var->get_raw());
	assert(*(const double *) max_var->get_raw() == 0);

	assert(!res_mat1->is_in_mem());
	assert(res_mat2->is_in_mem());

	// Perform MvTransMv.
	dense_matrix::ptr res3;
	{
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		res3 = mv1->MvTransMv(*res1);
		printf("The first MvTransMv\n");
		detail::matrix_stats.print_diff(orig_stats);
	}
	dense_matrix::ptr res4;
	{
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		res4 = mv2->MvTransMv(*res2);
		printf("The second MvTransMv\n");
		detail::matrix_stats.print_diff(orig_stats);
	}

	// Compare the result of MvTransMv.
	diff = res_mat1->minus(*res_mat2);
	max_var = diff->abs()->max();
	printf("max diff: %g\n", *(const double *) max_var->get_raw());
	assert(*(const double *) max_var->get_raw() == 0);
}

void test_copy(int num_nodes, bool in_mem)
{
	printf("test deep copy a dense matrix\n");
	dense_matrix::ptr mat = dense_matrix::create_randu<int>(0, 1000,
			long_dim, 10, matrix_layout_t::L_COL, num_nodes, in_mem);
	dense_matrix::ptr copy = mat->deep_copy();
	dense_matrix::ptr diff = mat->minus(*copy);
	scalar_variable::ptr max_var = diff->abs()->max();
	assert(*(const int *) max_var->get_raw() == 0);
}

void test_EM_persistent()
{
	std::string mat_name = "test.mat";
	printf("test creating a matrix from an existing matrix file\n");
	{
		// Remove the test matrix file.
		detail::EM_matrix_store::ptr store
			= detail::EM_matrix_store::create(mat_name);
		if (store)
			store->unset_persistent();
	}
	{
		// Create the test matrix file.
		dense_matrix::ptr mat = dense_matrix::create_randu<int>(0, 1000,
				long_dim, 10, matrix_layout_t::L_COL, -1, false);
		detail::EM_matrix_store::const_ptr store1
			= detail::EM_matrix_store::cast(mat->get_raw_store());
		bool ret = store1->set_persistent(mat_name);
		assert(ret);
		dense_matrix::ptr mat2 = dense_matrix::create(
				detail::EM_matrix_store::create(mat_name));
		assert(mat2);
		dense_matrix::ptr diff = mat->minus(*mat2);
		scalar_variable::ptr max_var = diff->abs()->max();
		assert(*(const int *) max_var->get_raw() == 0);
	}
	{
		// Access an existing matrix file.
		dense_matrix::ptr mat = dense_matrix::create(
				detail::EM_matrix_store::create(mat_name));
		assert(mat);
		dense_matrix::ptr test = mat->lt_scalar<int>(1001);
		scalar_variable::ptr sum = test->sum();
		assert(*(const size_t *) sum->get_raw()
				== mat->get_num_rows() * mat->get_num_cols());
		detail::EM_matrix_store::const_ptr store
			= detail::EM_matrix_store::cast(mat->get_raw_store());
		// Delete the matrix file.
		store->unset_persistent();
	}
	{
		// Test if the matrix file still exists.
		detail::EM_matrix_store::ptr store
			= detail::EM_matrix_store::create(mat_name);
		assert(store == NULL);
	}
}

void test_EM_matrix(int num_nodes)
{
	printf("test EM matrix\n");
	in_mem = false;

	matrix_val = matrix_val_t::SEQ;
	test_EM_persistent();
	test_sub_matrix();
	test_mapply_chain(-1, get_scalar_type<double>());
	test_mapply_chain(-1, get_scalar_type<int>());
	test_multiply_double(-1);
	test_cast();
	test_write2file();
	test_apply();
	test_conv_vec2mat();

	test_conv2(-1);
	test_scale_cols(-1);
	test_scale_rows(-1);
	test_multiply_scalar(-1);
	test_ele_wise(-1);
	test_multiply_col(-1);
	test_agg_col(-1);
	test_multiply_matrix(-1);
	test_agg_row(-1);
	test_agg_sub_col(-1);
	test_agg_sub_row(-1);
	test_sum_row_col(-1);
#if 0
	test_rand_init();
	test_conv_row_col();
	test_flatten();
#endif
}

void test_mem_matrix(int num_nodes)
{
	printf("test mem matrix\n");
	in_mem = true;

	matrix_val = matrix_val_t::SEQ;
	test_copy(-1, true);
	test_apply_scalar();
	test_min();
	test_mul_output(-1);
	test_mul_output(num_nodes);
	test_mapply_chain(-1, get_scalar_type<double>());
	test_mapply_chain(-1, get_scalar_type<int>());
	test_mapply_chain(num_nodes, get_scalar_type<int>());
	test_multiply_double(-1);
	test_cast();
	test_write2file();
	test_apply();
	test_conv_vec2mat();
	test_sum_row_col(-1);
	test_sum_row_col(num_nodes);

	matrix_val = matrix_val_t::SEQ;
	test_conv2(-1);
	test_conv2(num_nodes);
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
#if 0
	test_rand_init();
	test_conv_row_col();
	test_flatten();
#endif
}

void test_conv_store()
{
	in_mem = true;
	dense_matrix::ptr mat0 = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, -1, get_scalar_type<int>());
	assert(mat0->is_in_mem());
	printf("conv in-mem to EM\n");
	dense_matrix::ptr EM_mat = mat0->conv_store(false, -1);
	assert(!EM_mat->is_in_mem());
	printf("conv EM to in-mem\n");
	dense_matrix::ptr smp_mat = EM_mat->conv_store(true, -1);
	assert(smp_mat->is_in_mem());
	verify_result(*mat0, *smp_mat, equal_func<int>());

	if (matrix_conf.get_num_nodes() > 1) {
		printf("conv SMP to NUMA\n");
		dense_matrix::ptr numa_mat = mat0->conv_store(true,
				matrix_conf.get_num_nodes());
		verify_result(*mat0, *numa_mat, equal_func<int>());

		printf("conv NUMA to EM\n");
		EM_mat = numa_mat->conv_store(false, -1);
		assert(!EM_mat->is_in_mem());

		printf("conv EM to NUMA\n");
		numa_mat = EM_mat->conv_store(true, matrix_conf.get_num_nodes());
		assert(numa_mat->is_in_mem());
		verify_result(*mat0, *numa_mat, equal_func<int>());
	}
}

void test_bmv_multiply_tall()
{
	bool in_mem = false;
	printf("gemm tall on block multi-vector\n");
	eigen::block_multi_vector::ptr mv = eigen::block_multi_vector::create(
			long_dim, 14, 2, get_scalar_type<double>(), in_mem);
	for (size_t i = 0; i < mv->get_num_blocks(); i++)
		mv->set_block(i, create_seq_matrix(long_dim, mv->get_block_size(),
					matrix_layout_t::L_COL, -1, get_scalar_type<double>(),
					in_mem));
	dense_matrix::ptr mat = create_seq_matrix(mv->get_num_cols(),
			mv->get_block_size(), matrix_layout_t::L_COL, -1,
			get_scalar_type<double>(), true);
	mat->materialize_self();
	detail::mem_col_matrix_store::const_ptr B
		= detail::mem_col_matrix_store::cast(mat->get_raw_store());
	scalar_variable_impl<double> alpha(2);
	scalar_variable_impl<double> beta(3);

	{
		printf("gemm1\n");
		eigen::block_multi_vector::ptr res1 = eigen::block_multi_vector::create(
				long_dim, mv->get_block_size(), mv->get_block_size(),
				get_scalar_type<double>(), true);
		res1->set_block(0, create_seq_matrix(long_dim, mv->get_block_size(),
					matrix_layout_t::L_COL, -1, get_scalar_type<double>(), in_mem));
		res1->set_multiply_blocks(2);
		res1 = res1->gemm(*mv, B, alpha, beta);
		assert(res1->get_num_blocks() == 1);
		dense_matrix::ptr res_mat1 = res1->get_block(0);
		res_mat1->materialize_self();

		printf("gemm2\n");
		eigen::block_multi_vector::ptr res2 = eigen::block_multi_vector::create(
				long_dim, mv->get_block_size(), mv->get_block_size(),
				get_scalar_type<double>(), true);
		res2->set_block(0, create_seq_matrix(long_dim, mv->get_block_size(),
					matrix_layout_t::L_COL, -1, get_scalar_type<double>(), in_mem));
		res2->set_multiply_blocks(mv->get_num_blocks());
		res2 = res2->gemm(*mv, B, alpha, beta);
		assert(res2->get_num_blocks() == 1);
		dense_matrix::ptr res_mat2 = res2->get_block(0);
		res_mat2->materialize_self();

		dense_matrix::ptr diff = res_mat1->minus(*res_mat2);
		scalar_variable::ptr max_diff = diff->abs()->max();
		scalar_variable::ptr max1 = res_mat1->max();
		scalar_variable::ptr max2 = res_mat2->max();
		printf("max diff: %g, max mat1: %g, max mat2: %g\n",
				*(double *) max_diff->get_raw(), *(double *) max1->get_raw(),
				*(double *) max2->get_raw());
		assert(*(double *) max_diff->get_raw() == 0);
	}

	mat = create_seq_matrix(mv->get_num_cols(),
			mv->get_block_size() * 2, matrix_layout_t::L_COL, -1,
			get_scalar_type<double>(), true);
	mat->materialize_self();
	B = detail::mem_col_matrix_store::cast(mat->get_raw_store());

	{
		printf("gemm1\n");
		eigen::block_multi_vector::ptr res1 = eigen::block_multi_vector::create(
				long_dim, B->get_num_cols(), mv->get_block_size(),
				get_scalar_type<double>(), true);
		std::vector<detail::matrix_store::const_ptr> orig_res_stores(res1->get_num_blocks());
		for (size_t i = 0; i < res1->get_num_blocks(); i++) {
			res1->set_block(i, create_seq_matrix(long_dim, mv->get_block_size(),
						matrix_layout_t::L_COL, -1, get_scalar_type<double>(),
						in_mem));
			orig_res_stores[i] = res1->get_block(i)->get_raw_store();
		}
		res1->set_multiply_blocks(2);
		res1 = res1->gemm(*mv, B, alpha, beta);

		printf("gemm2\n");
		eigen::block_multi_vector::ptr res2 = eigen::block_multi_vector::create(
				long_dim, B->get_num_cols(), B->get_num_cols(),
				get_scalar_type<double>(), true);
		res2->set_block(0, dense_matrix::create(
					eigen::collected_matrix_store::create(orig_res_stores,
						res1->get_num_cols())));
		res2->set_multiply_blocks(mv->get_num_blocks());
		res2 = res2->gemm(*mv, B, alpha, beta);
		assert(res2->get_num_blocks() == 1);
		res2->get_block(0)->materialize_self();

		off_t col_off = 0;
		for (size_t i = 0; i < res1->get_num_blocks(); i++) {
			std::vector<off_t> col_idxs(res1->get_block_size());
			for (size_t j = 0; j < res1->get_block_size(); j++)
				col_idxs[j] = col_off++;
			dense_matrix::ptr res_mat1 = res1->get_block(i);
			dense_matrix::ptr res_mat2 = res2->get_block(0)->get_cols(col_idxs);

			dense_matrix::ptr diff = res_mat1->minus(*res_mat2);
			scalar_variable::ptr max_diff = diff->abs()->max();
			scalar_variable::ptr max1 = res_mat1->max();
			scalar_variable::ptr max2 = res_mat2->max();
			printf("max diff: %g, max mat1: %g, max mat2: %g\n",
					*(double *) max_diff->get_raw(), *(double *) max1->get_raw(),
					*(double *) max2->get_raw());
			assert(*(double *) max_diff->get_raw() == 0);
		}
	}
}

void test_bmv_multiply_wide()
{
	bool in_mem = false;
	printf("gemm wide on block multi-vector\n");
	eigen::block_multi_vector::ptr mv1 = eigen::block_multi_vector::create(
			long_dim, 14, 2, get_scalar_type<double>(), in_mem);
	for (size_t i = 0; i < mv1->get_num_blocks(); i++) {
		dense_matrix::ptr mat = dense_matrix::create_randu<int>(0, 10, long_dim,
				mv1->get_block_size(), matrix_layout_t::L_COL, -1, in_mem);
		mv1->set_block(i, mat->cast_ele_type(get_scalar_type<double>()));
	}

	eigen::block_multi_vector::ptr mv2 = eigen::block_multi_vector::create(
			long_dim, 2, 2, get_scalar_type<double>(), in_mem);
	for (size_t i = 0; i < mv2->get_num_blocks(); i++) {
		dense_matrix::ptr mat = dense_matrix::create_randu<int>(0, 10, long_dim,
				mv2->get_block_size(), matrix_layout_t::L_COL, -1, in_mem);
		mv2->set_block(i, mat->cast_ele_type(get_scalar_type<double>()));
	}

	printf("gemm0\n");
	mv1->set_multiply_blocks(mv1->get_num_blocks());
	dense_matrix::ptr res0 = mv1->MvTransMv(*mv2);

	printf("gemm1\n");
	mv1->set_multiply_blocks(2);
	dense_matrix::ptr res1 = mv1->MvTransMv(*mv2);

	printf("gemm2\n");
	mv2->set_multiply_blocks(2);
	dense_matrix::ptr res2 = mv2->MvTransMv(*mv1);
	res2 = res2->transpose();

	assert(res0->get_num_rows() == res1->get_num_rows());
	assert(res0->get_num_cols() == res1->get_num_cols());
	assert(res0->get_num_rows() == res2->get_num_rows());
	assert(res0->get_num_cols() == res2->get_num_cols());

	dense_matrix::ptr diff1 = res0->minus(*res1);
	scalar_variable::ptr max_diff1 = diff1->abs()->max();
	dense_matrix::ptr diff2 = res0->minus(*res2);
	scalar_variable::ptr max_diff2 = diff2->abs()->max();
	printf("max diff1: %g, max diff2: %g\n",
			*(double *) max_diff1->get_raw(),
			*(double *) max_diff2->get_raw());
	assert(*(double *) max_diff1->get_raw() == 0);
	assert(*(double *) max_diff2->get_raw() == 0);
}

void test_block_mv()
{
	test_bmv_multiply_wide();
	test_bmv_multiply_tall();
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
	int num_nodes = matrix_conf.get_num_nodes();

	test_block_mv();
	test_conv_store();
	test_mapply_mixed(num_nodes);
	test_mem_matrix(num_nodes);
	test_EM_matrix(num_nodes);

	destroy_flash_matrix();
}
