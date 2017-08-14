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
#include "mapply_matrix_store.h"
#include "factor.h"
#include "block_matrix.h"
#include "project_matrix_store.h"
#include "data_frame.h"
#include "mem_vec_store.h"

#include "eigensolver/block_dense_matrix.h"
#include "eigensolver/collected_col_matrix_store.h"

using namespace fm;

size_t long_dim = 99999;

/*
 * This is a naive implementation of matrix multiplication.
 * It should be correct
 */
dense_matrix::ptr naive_multiply(const dense_matrix &_m1, const dense_matrix &_m2)
{
	dense_matrix::ptr m1 = dense_matrix::create(_m1.get_raw_store());
	dense_matrix::ptr m2 = dense_matrix::create(_m2.get_raw_store());
	m1->materialize_self();
	m2->materialize_self();
	detail::mem_matrix_store::ptr res_store = detail::mem_matrix_store::create(
			m1->get_num_rows(), m2->get_num_cols(), matrix_layout_t::L_ROW,
			get_scalar_type<int>(), -1);
	detail::mem_matrix_store::const_ptr mem_m1;
	detail::mem_matrix_store::const_ptr mem_m2;
	if (m1->is_in_mem())
		mem_m1 = detail::mem_matrix_store::cast(m1->get_raw_store());
	else {
		dense_matrix::ptr mem_mat = m1->conv_store(true, -1);
		mem_m1 = detail::mem_matrix_store::cast(mem_mat->get_raw_store());
	}
	if (m2->is_in_mem())
		mem_m2 = detail::mem_matrix_store::cast(m2->get_raw_store());
	else {
		dense_matrix::ptr mem_mat = m2->conv_store(true, -1);
		mem_m2 = detail::mem_matrix_store::cast(mem_mat->get_raw_store());
	}
#pragma omp parallel for
	for (size_t i = 0; i < m1->get_num_rows(); i++) {
		for (size_t j = 0; j < m2->get_num_cols(); j++) {
			int sum = 0;
			for (size_t k = 0; k < m1->get_num_cols(); k++) {
				sum += mem_m1->get<int>(i, k) * mem_m2->get<int>(k, j);
			}
			res_store->set<int>(i, j, sum);
		}
	}
	return dense_matrix::create(res_store);
}

dense_matrix::ptr blas_multiply(const dense_matrix &m1, const dense_matrix &m2)
{
	assert(m1.get_type() == m2.get_type());
	dense_matrix::ptr tmp1 = dense_matrix::create(m1.get_raw_store());
	dense_matrix::ptr tmp2 = dense_matrix::create(m2.get_raw_store());
	tmp1 = tmp1->conv2(matrix_layout_t::L_COL);
	tmp2 = tmp2->conv2(matrix_layout_t::L_COL);
	tmp1->materialize_self();
	tmp2->materialize_self();
	detail::mem_col_matrix_store::ptr col_res = detail::mem_col_matrix_store::create(
			tmp1->get_num_rows(), tmp2->get_num_cols(), m1.get_type());
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
	if (m1.get_type() == get_scalar_type<double>()) {
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
				mem_m1->get_num_rows(), mem_m2->get_num_cols(),
				mem_m1->get_num_cols(), 1,
				(const double *) col_m1->get_data().get_raw(),
				mem_m1->get_num_rows(),
				(const double *) col_m2->get_data().get_raw(),
				mem_m2->get_num_rows(), 0,
				(double *) col_res->get_data().get_raw(),
				col_res->get_num_rows());
	}
	else {
		assert(m1.get_type() == get_scalar_type<float>());
		cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
				mem_m1->get_num_rows(), mem_m2->get_num_cols(),
				mem_m1->get_num_cols(), 1,
				(const float *) col_m1->get_data().get_raw(),
				mem_m1->get_num_rows(),
				(const float *) col_m2->get_data().get_raw(),
				mem_m2->get_num_rows(), 0,
				(float *) col_res->get_data().get_raw(),
				col_res->get_num_rows());
	}
	return dense_matrix::create(col_res);
}

template<class T>
T get_precision()
{
	return 0;
}

template<>
float get_precision<float>()
{
	return 1e-5;
}

template<>
double get_precision<double>()
{
	return 1e-13;
}

template<class T>
struct approx_equal_func
{
public:
	bool operator()(const char *raw1, const char *raw2) const {
		T v1 = *(const T *) raw1;
		T v2 = *(const T *) raw2;
		T diff = v1 - v2;
		if (diff == 0)
			return true;

		if (v1 < 0)
			v1 = -v1;
		if (v2 < 0)
			v2 = -v2;
		if (diff < 0)
			diff = -diff;
		diff /= std::min(v1, v2);
		bool ret = diff < get_precision<T>();
		if (!ret) {
			printf("v1: %f, v2: %f, diff: %g\n", v1, v2, diff);
		}
		return ret;
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
void verify_result(const dense_matrix &_m1, const dense_matrix &_m2,
		const Func &func)
{
	// It's possible the input matrices are block matrices.
	dense_matrix::ptr m1 = dense_matrix::create(_m1.get_raw_store());
	dense_matrix::ptr m2 = dense_matrix::create(_m2.get_raw_store());

	assert(m1->get_num_rows() == m2->get_num_rows());
	assert(m1->get_num_cols() == m2->get_num_cols());

	m1->materialize_self();
	m2->materialize_self();

	detail::mem_matrix_store::const_ptr mem_m1;
	detail::mem_matrix_store::const_ptr mem_m2;
	if (m1->is_in_mem() && !m1->is_virtual())
		mem_m1 = detail::mem_matrix_store::cast(m1->get_raw_store());
	else {
		dense_matrix::ptr mem_mat = m1->conv_store(true, -1);
		mem_m1 = detail::mem_matrix_store::cast(mem_mat->get_raw_store());
	}
	if (m2->is_in_mem() && !m2->is_virtual())
		mem_m2 = detail::mem_matrix_store::cast(m2->get_raw_store());
	else {
		dense_matrix::ptr mem_mat = m2->conv_store(true, -1);
		mem_m2 = detail::mem_matrix_store::cast(mem_mat->get_raw_store());
	}

	bool success = true;
#pragma omp parallel for
	for (size_t i = 0; i < m1->get_num_rows(); i++)
		for (size_t j = 0; j < m1->get_num_cols(); j++)
			if (!func(mem_m1->get(i, j), mem_m2->get(i, j))) {
				success = false;
				break;
			}
	if (!success) {
		m1->print();
		m2->print();
	}
}

enum matrix_val_t
{
	SEQ,
	SEQ_MATER,
	DEFAULT,
	SPARSE,
	NUM_TYPES,
} matrix_val = matrix_val_t::SEQ;

dense_matrix::ptr create_seq_matrix(size_t nrow, size_t ncol,
		matrix_layout_t layout, int num_nodes, const scalar_type &type,
		bool in_mem)
{
	if (type == get_scalar_type<int>())
		return dense_matrix::create_seq<int>(0, 1, nrow, ncol, layout, true,
				num_nodes, in_mem);
	else if (type == get_scalar_type<size_t>())
		return dense_matrix::create_seq<size_t>(0, 1, nrow, ncol, layout, true,
				num_nodes, in_mem);
	else if (type == get_scalar_type<double>())
		return dense_matrix::create_seq<size_t>(0, 1, nrow, ncol, layout, true,
				num_nodes, in_mem)->cast_ele_type(get_scalar_type<double>());
	else if (type == get_scalar_type<float>())
		return dense_matrix::create_seq<size_t>(0, 1, nrow, ncol, layout, true,
				num_nodes, in_mem)->cast_ele_type(get_scalar_type<float>());
	else
		return dense_matrix::ptr();
}

template<class T>
dense_matrix::ptr _create_seq_block_matrix(size_t nrow, size_t ncol,
		size_t block_size, matrix_layout_t layout, int num_nodes, bool in_mem)
{
	scalar_variable::ptr start_val(new scalar_variable_impl<T>(0));
	scalar_variable::ptr stride_val(new scalar_variable_impl<T>(1));
	return block_matrix::create_seq_layout(start_val, stride_val,
			nrow, ncol, layout, block_size, true, num_nodes, in_mem);
}

dense_matrix::ptr create_seq_block_matrix(size_t nrow, size_t ncol,
		size_t block_size, matrix_layout_t layout, int num_nodes,
		const scalar_type &type, bool in_mem)
{
	if (type == get_scalar_type<int>())
		return _create_seq_block_matrix<int>(nrow, ncol, block_size, layout,
				num_nodes, in_mem);
	else if (type == get_scalar_type<size_t>())
		return _create_seq_block_matrix<size_t>(nrow, ncol, block_size, layout,
				num_nodes, in_mem);
	else if (type == get_scalar_type<double>())
		return _create_seq_block_matrix<size_t>(nrow, ncol, block_size, layout,
				num_nodes, in_mem)->cast_ele_type(get_scalar_type<double>());
	else if (type == get_scalar_type<float>())
		return _create_seq_block_matrix<size_t>(nrow, ncol, block_size, layout,
				num_nodes, in_mem)->cast_ele_type(get_scalar_type<float>());
	else
		return dense_matrix::ptr();
}

dense_matrix::ptr create_const_matrix(size_t nrow, size_t ncol,
		matrix_layout_t layout, int num_nodes, const scalar_type &type,
		bool in_mem)
{
	if (type == get_scalar_type<int>())
		return dense_matrix::create_const<int>(1, nrow, ncol, layout,
				num_nodes, in_mem);
	else if (type == get_scalar_type<size_t>())
		return dense_matrix::create_const<size_t>(1, nrow, ncol, layout,
				num_nodes, in_mem);
	else if (type == get_scalar_type<double>())
		return dense_matrix::create_const<double>(1, nrow, ncol, layout,
				num_nodes, in_mem);
	else if (type == get_scalar_type<float>())
		return dense_matrix::create_const<float>(1, nrow, ncol, layout,
				num_nodes, in_mem);
	else
		return dense_matrix::ptr();
}

dense_matrix::ptr create_sparse_matrix(size_t nrow, size_t ncol,
		matrix_layout_t layout, int num_nodes, const scalar_type &type,
		bool in_mem)
{
	if (type == get_scalar_type<int>())
		return dense_matrix::create(detail::sparse_project_matrix_store::create_sparse_rand(
				nrow, ncol, layout, get_scalar_type<int>(), 0.001));
	else if (type == get_scalar_type<size_t>())
		return dense_matrix::create(detail::sparse_project_matrix_store::create_sparse_rand(
					nrow, ncol, layout, get_scalar_type<int>(),
					0.001))->cast_ele_type(get_scalar_type<size_t>());
	else if (type == get_scalar_type<double>())
		return dense_matrix::create(detail::sparse_project_matrix_store::create_sparse_rand(
				nrow, ncol, layout, get_scalar_type<double>(), 0.001));
	else if (type == get_scalar_type<float>())
		return dense_matrix::create(detail::sparse_project_matrix_store::create_sparse_rand(
					nrow, ncol, layout, get_scalar_type<int>(),
					0.001))->cast_ele_type(get_scalar_type<float>());
	else
		return dense_matrix::ptr();
}


bool in_mem = true;
size_t block_size = 0;

dense_matrix::ptr create_matrix(size_t nrow, size_t ncol,
		matrix_layout_t layout, int num_nodes,
		const scalar_type &type = get_scalar_type<int>())
{
	dense_matrix::ptr mat;
	switch (matrix_val) {
		case matrix_val_t::DEFAULT:
			return create_const_matrix(nrow, ncol, layout, num_nodes, type,
					in_mem);
		case matrix_val_t::SEQ:
			if (block_size == 0)
				mat = create_seq_matrix(nrow, ncol, layout, num_nodes, type,
						in_mem);
			else
				mat = create_seq_block_matrix(nrow, ncol, block_size, layout,
						num_nodes, type, in_mem);
			if (!in_mem)
				mat = mat->conv_store(false, -1);
			return mat;
		case matrix_val_t::SEQ_MATER:
			if (block_size == 0)
				mat = create_seq_matrix(nrow, ncol, layout, num_nodes, type,
						in_mem);
			else
				mat = create_seq_block_matrix(nrow, ncol, block_size, layout,
						num_nodes, type, in_mem);
			mat = mat->conv_store(in_mem, num_nodes);
			return mat;
		case matrix_val_t::SPARSE:
			return create_sparse_matrix(nrow, ncol, layout, num_nodes, type,
					in_mem);
		default:
			assert(0);
			return dense_matrix::ptr();
	}
}

void test_mapply_sink(dense_matrix::ptr sink)
{
	// Mapply on a sink matrix.
	printf("Test mapply on a sink matrix\n");
	assert(sink->get_raw_store()->is_sink());
	sink = sink->cast_ele_type(get_scalar_type<double>());
	dense_matrix tmp = *sink + *sink;
	tmp = tmp * 0.5;
	tmp.materialize_self();
	sink->materialize_self();
	verify_result(*sink, tmp, equal_func<double>());
}

void test_mapply_sink_t(dense_matrix::ptr sink)
{
	printf("Test mapply on a sink matrix and transpose\n");
	assert(sink->get_raw_store()->is_sink());
	sink = sink->cast_ele_type(get_scalar_type<double>());
	dense_matrix tmp = *sink + *sink;
	tmp = tmp * 0.5;
	sink = sink->transpose();
	dense_matrix::ptr tmp1 = tmp.transpose();
	tmp1->materialize_self();
	sink->materialize_self();
	verify_result(*sink, *tmp1, equal_func<double>());
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
	dense_matrix::ptr res = m1->add(*m1);
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
	assert(res1->is_in_mem() == m1->is_in_mem());
	verify_result(*res1, *correct, equal_func<int>());
}

template<class T>
void test_multiply(int num_nodes)
{
	dense_matrix::ptr m1, m2, correct, res;

	printf("Test multiplication on wide block matrix X tall block matrix\n");
	size_t orig_multiply_block_size = matrix_conf.get_max_multiply_block_size();
	block_size = 3;
	matrix_conf.set_max_multiply_block_size(6);
	m1 = create_matrix(20, long_dim, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	m2 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	correct = blas_multiply(*m1, *m2);
	res = m1->multiply(*m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test self cross prod\n");
	m2 = m1->transpose();
	correct = blas_multiply(*m1, *m2);
	res = m1->multiply(*m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	res = m1->multiply(*m2);
	test_mapply_sink(res);
	res = m1->multiply(*m2);
	test_mapply_sink_t(res);

	printf("Test transpose of self cross prod\n");
	res = m1->multiply(*m2);
	res = res->transpose();
	correct = blas_multiply(*m2->transpose(), *m1->transpose());
	verify_result(*res, *correct, approx_equal_func<T>());
	matrix_conf.set_max_multiply_block_size(orig_multiply_block_size);
	block_size = 0;

	std::vector<off_t> idxs(3);
	idxs[0] = 1;
	idxs[1] = 3;
	idxs[2] = 5;

	std::vector<dense_matrix::ptr> tmp_vec(1);
	printf("Test self cross product on wide row matrix\n");
	m1 = create_matrix(10, long_dim, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<T>());
	m2 = m1->transpose();
	res = m1->multiply(*m2);
	tmp_vec[0] = res;
	materialize(tmp_vec);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test self cross product on wide row submatrix\n");
	m1->materialize_self();
	m1 = m1->get_rows(idxs);
	m2 = m1->transpose();
	if (m1->is_in_mem())
		m1 = m1->deep_copy();
	if (m2->is_in_mem())
		m2 = m2->deep_copy();
	res = m1->multiply(*m2);
	tmp_vec[0] = res;
	materialize(tmp_vec);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test self cross product on wide col matrix\n");
	m1 = create_matrix(10, long_dim, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	m2 = m1->transpose();
	res = m1->multiply(*m2);
	res->materialize_self();
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test multiplication on tall row matrix X small row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<T>());
	m2 = create_matrix(10, 9, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<T>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE);
	assert(res->is_in_mem() == m1->is_in_mem());
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test multiplication on tall row matrix X small row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<T>());
	m2 = create_matrix(10, 9, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<T>());
	res = m1->multiply(*m2, matrix_layout_t::L_COL);
	assert(res->is_in_mem() == m1->is_in_mem());
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test multiplication on tall row matrix X small column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<T>());
	m2 = create_matrix(10, 9, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test multiplication on tall column matrix X small row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	m2 = create_matrix(10, 9, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<T>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test multiplication on tall column matrix X small row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	m2 = create_matrix(10, 9, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<T>());
	res = m1->multiply(*m2, matrix_layout_t::L_ROW);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test multiplication on tall column matrix X small column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	m2 = create_matrix(10, 9, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE);
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test multiplication on wide row matrix X tall column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	m1 = m1->transpose();
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE);
	res->materialize_self();
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test multiplication on wide row matrix X tall column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	m1 = m1->transpose();
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	res = m1->multiply(*m2, matrix_layout_t::L_COL);
	res->materialize_self();
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test multiplication on wide row matrix X tall row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	m1 = m1->transpose();
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<T>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE);
	res->materialize_self();
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test multiplication on wide column matrix X tall column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<T>());
	m1 = m1->transpose();
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE);
	res->materialize_self();
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test multiplication on wide column matrix X tall column matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<T>());
	m1 = m1->transpose();
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<T>());
	res = m1->multiply(*m2, matrix_layout_t::L_ROW);
	res->materialize_self();
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());

	printf("Test multiplication on wide column matrix X tall row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<T>());
	m1 = m1->transpose();
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<T>());
	res = m1->multiply(*m2, matrix_layout_t::L_NONE);
	res->materialize_self();
	correct = blas_multiply(*m1, *m2);
	verify_result(*res, *correct, approx_equal_func<T>());
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
	res = m1->multiply(*m2);
	verify_result(*res->transpose(), *correct->transpose(), equal_func<int>());

	printf("Test multiplication on wide row matrix X tall row matrix\n");
	m1 = create_matrix(10, long_dim, matrix_layout_t::L_ROW, num_nodes);
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_ROW, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = m1->multiply(*m2);
	verify_result(*res, *correct, equal_func<int>());
	res = m1->multiply(*m2);
	verify_result(*res->transpose(), *correct->transpose(), equal_func<int>());

	printf("Test multiplication on wide column matrix X tall column matrix\n");
	m1 = create_matrix(10, long_dim, matrix_layout_t::L_COL, num_nodes);
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_COL, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = m1->multiply(*m2);
	verify_result(*res, *correct, equal_func<int>());
	res = m1->multiply(*m2);
	verify_result(*res->transpose(), *correct->transpose(), equal_func<int>());

	printf("Test multiplication on wide column matrix X tall row matrix\n");
	m1 = create_matrix(10, long_dim, matrix_layout_t::L_COL, num_nodes);
	m2 = create_matrix(long_dim, 9, matrix_layout_t::L_ROW, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = m1->multiply(*m2);
	verify_result(*res, *correct, equal_func<int>());
	res = m1->multiply(*m2);
	verify_result(*res->transpose(), *correct->transpose(), equal_func<int>());

	printf("Test multiplication on tall row matrix X small row matrix\n");
	m1 = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes);
	m2 = create_matrix(10, 9, matrix_layout_t::L_ROW, num_nodes);
	correct = naive_multiply(*m1, *m2);
	res = m1->multiply(*m2);
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

void test_agg(int num_nodes, matrix_layout_t layout)
{
	printf("Test aggregation on tall %s-major matrix\n",
			layout == matrix_layout_t::L_COL ? "column" : "row");
	// Aggregate on the entire matrix and output a single value.
	dense_matrix::ptr m1 = create_matrix(long_dim, 10,
			layout, num_nodes, get_scalar_type<size_t>());
	bulk_operate::const_ptr op = bulk_operate::conv2ptr(
			m1->get_type().get_basic_ops().get_add());
	scalar_variable::ptr res = m1->aggregate(op);
	assert(res->get_type() == m1->get_type());
	assert(res->get_type() == get_scalar_type<size_t>());
	size_t sum = *(size_t *) res->get_raw();
	size_t num_eles = m1->get_num_rows() * m1->get_num_cols();
	if (matrix_val == matrix_val_t::DEFAULT)
		assert(sum == m1->get_num_rows() * m1->get_num_cols());
	else if (matrix_val_t::SEQ)
		assert(sum == (num_eles - 1) * num_eles / 2);

	// Aggregate on rows of a matrix and output a column
	dense_matrix::ptr sum_col = m1->aggregate(matrix_margin::MAR_ROW,
			agg_operate::create(op));
	if (sum_col->is_in_mem())
		sum_col->set_materialize_level(materialize_level::MATER_FULL);
	else {
		detail::matrix_store::ptr buf = detail::mem_matrix_store::create(
				sum_col->get_num_rows(), sum_col->get_num_cols(),
				sum_col->store_layout(), sum_col->get_type(), num_nodes);
		sum_col->set_materialize_level(materialize_level::MATER_FULL, buf);
	}
	res = sum_col->aggregate(op);
	sum_col->materialize_self();
	assert(sum_col->get_num_rows() == m1->get_num_rows());
	assert(sum_col->get_num_cols() == 1);
	detail::mem_matrix_store::const_ptr tmp
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				sum_col->get_raw_store());
	assert(tmp);
	for (size_t i = 0; i < sum_col->get_num_rows(); i++) {
		size_t ncol = m1->get_num_cols();
		if (matrix_val == matrix_val_t::SEQ)
			assert(tmp->get<size_t>(i, 0)
					== i * ncol * ncol + (ncol - 1) * ncol / 2);
		else if (matrix_val == matrix_val_t::DEFAULT)
			assert(tmp->get<size_t>(i, 0) == m1->get_num_cols());
	}

	// Aggregate on cols of a matrix.
	dense_matrix::ptr sum_row = m1->aggregate(matrix_margin::MAR_COL,
			agg_operate::create(op));
	sum_col = sum_row->transpose();
	sum_row->materialize_self();
	sum_col->materialize_self();
	detail::mem_matrix_store::const_ptr tmp1
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				sum_row->get_raw_store());
	detail::mem_matrix_store::const_ptr tmp2
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				sum_col->get_raw_store());
	assert(tmp1);
	assert(tmp2);
	assert(tmp1->get_raw_arr());
	assert(tmp2->get_raw_arr());
	num_eles = tmp1->get_num_rows() * tmp1->get_num_cols();
	assert(memcmp(tmp1->get_raw_arr(), tmp2->get_raw_arr(),
				num_eles * tmp1->get_entry_size()) == 0);
	dense_matrix::ptr tmp_mat = dense_matrix::create(tmp1);
	res = tmp_mat->aggregate(op);
	sum = *(size_t *) res->get_raw();
	num_eles = m1->get_num_rows() * m1->get_num_cols();
	if (matrix_val == matrix_val_t::DEFAULT)
		assert(sum == m1->get_num_rows() * m1->get_num_cols());
	else if (matrix_val_t::SEQ)
		assert(sum == (num_eles - 1) * num_eles / 2);

	// Mapply on a sink matrix.
	printf("Test mapply on an agg matrix\n");
	dense_matrix::ptr sink = m1->aggregate(matrix_margin::MAR_COL, agg_operate::create(op));
	test_mapply_sink(sink);
	sink = m1->aggregate(matrix_margin::MAR_COL, agg_operate::create(op));
	test_mapply_sink_t(sink);
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

	bulk_operate::const_ptr op = bulk_operate::conv2ptr(
			sub_m->get_type().get_basic_ops().get_add());
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
		assert(sum == sub_m->get_num_rows() * sub_m->get_num_cols());
	else if (matrix_val_t::SEQ)
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

	bulk_operate::const_ptr op = bulk_operate::conv2ptr(
			sub_col_m->get_type().get_basic_ops().get_add());
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
	vector::ptr vals = create_seq_vector<int>(0, orig->get_num_cols() - 1, 1);
	dense_matrix::ptr res = orig->scale_cols(col_vec::create(vals));
	assert(res->is_virtual());
	assert(res->is_in_mem() == orig->is_in_mem());
	res = dense_matrix::create(res->get_raw_store());
	res->materialize_self();
	orig = dense_matrix::create(orig->get_raw_store());
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
	vector::ptr vals = create_seq_vector<int>(0, orig->get_num_rows() - 1, 1);
	dense_matrix::ptr res = orig->scale_rows(col_vec::create(vals));
	assert(res->is_virtual());
	assert(res->is_in_mem() == orig->is_in_mem());
	res = dense_matrix::create(res->get_raw_store());
	res->materialize_self();
	orig = dense_matrix::create(orig->get_raw_store());
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
	virtual size_t get_num_out_eles(size_t num_input) const {
		return 1;
	}
};

void test_apply1(dense_matrix::ptr mat)
{
	size_t num_rows = mat->get_num_rows();
	size_t num_cols = mat->get_num_cols();
	dense_matrix::ptr res;

	printf("Test apply on rows of a %s %s-wise matrix\n",
			mat->is_wide() ? "wide" : "tall",
			mat->store_layout() == matrix_layout_t::L_ROW ? "row" : "column");
	res = mat->apply(matrix_margin::MAR_ROW,
			arr_apply_operate::const_ptr(new sum_apply_op()));
	if (res) {
		assert(res->get_num_cols() == 1 && res->get_num_rows() == mat->get_num_rows());
		assert(res->is_type<long>());
		res->materialize_self();
		col_vec::ptr res_vec = col_vec::create(res);
		std::vector<long> stdvec = res_vec->conv2std<long>();
		for (size_t i = 0; i < stdvec.size(); i++)
			assert(stdvec[i]
					== i * num_cols * num_cols + (num_cols - 1) * num_cols / 2);
	}

	printf("Test apply on columns of a %s %s-wise matrix\n",
			mat->is_wide() ? "wide" : "tall",
			mat->store_layout() == matrix_layout_t::L_ROW ? "row" : "column");
	res = mat->apply(matrix_margin::MAR_COL,
			arr_apply_operate::const_ptr(new sum_apply_op()));
	if (res) {
		assert(res->get_num_rows() == 1 && res->get_num_cols() == mat->get_num_cols());
		assert(res->is_type<long>());
		res->materialize_self();
		col_vec::ptr res_vec = col_vec::create(res);
		std::vector<long> stdvec = res_vec->conv2std<long>();
		for (size_t i = 0; i < stdvec.size(); i++)
			assert(stdvec[i]
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
	vector::ptr vec = create_seq_vector<int>(0, len, 1);
	dense_matrix::ptr mat = vec->conv2mat(num_rows, len / num_rows, true);
	assert(mat);
	assert(mat->store_layout() == matrix_layout_t::L_ROW);

	mat = vec->conv2mat(num_rows, len / num_rows, false);
	assert(mat);
	assert(mat->store_layout() == matrix_layout_t::L_COL);
}

void test_write2file1(detail::mem_matrix_store::const_ptr mat)
{
	char *tmp_file_name = tempnam(".", "tmp.mat");
	bool ret = mat->write2file(tmp_file_name);
	assert(ret);

	detail::mem_matrix_store::const_ptr read_mat = detail::mem_matrix_store::load(
			tmp_file_name, mat->get_num_nodes());
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
	dense_matrix::ptr mat;
	detail::mem_matrix_store::const_ptr store;

	printf("write a tall row matrix\n");
	mat = dense_matrix::create_seq<int>(0, 1, 1000000, 10,
			matrix_layout_t::L_ROW, true, -1, true);
	mat->materialize_self();
	store = detail::mem_matrix_store::cast(mat->get_raw_store());
	test_write2file1(store);
	printf("write a tall column matrix\n");
	mat = dense_matrix::create_seq<int>(0, 1, 1000000, 10,
			matrix_layout_t::L_COL, true, -1, true);
	mat->materialize_self();
	store = detail::mem_matrix_store::cast(mat->get_raw_store());
	test_write2file1(store);
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
	dense_matrix::ptr mat = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes);

	printf("conv2 layout of mem matrix of %ld, %ld\n", mat->get_num_rows(),
			mat->get_num_cols());
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

	mat = create_matrix(long_dim, 1, matrix_layout_t::L_COL, num_nodes);
	printf("conv2 layout of mem matrix of %ld, %ld\n", mat->get_num_rows(),
			mat->get_num_cols());
	mat1 = mat->conv2(matrix_layout_t::L_ROW);
	assert(mat->get_num_rows() == mat1->get_num_rows());
	assert(mat->get_num_cols() == mat1->get_num_cols());
	assert(mat1->store_layout() == matrix_layout_t::L_ROW);
	verify_result(*mat, *mat1, equal_func<int>());

	mat = mat->add_scalar<int>(1);
	assert(mat->is_virtual());
	printf("conv2 layout of virtual matrix of %ld, %ld\n", mat->get_num_rows(),
			mat->get_num_cols());
	mat1 = mat->conv2(matrix_layout_t::L_ROW);
	assert(mat1->is_in_mem() == mat->is_in_mem());
	assert(mat->get_num_rows() == mat1->get_num_rows());
	assert(mat->get_num_cols() == mat1->get_num_cols());
	assert(mat1->store_layout() == matrix_layout_t::L_ROW);
	verify_result(*mat, *mat1, equal_func<int>());
}

void test_sum_row_col1(dense_matrix::ptr mat)
{
	printf("test row sum on %s %s matrix\n", mat->is_wide() ? "wide" : "tall",
			mat->store_layout() == matrix_layout_t::L_COL ? "column" : "row");
	dense_matrix::ptr sum = mat->row_sum();
	sum->materialize_self();
	col_vec::ptr vec = col_vec::create(sum);
	assert(vec->get_length() == mat->get_num_rows());

	std::vector<int> stdvec = vec->conv2std<int>();
	detail::mem_matrix_store::const_ptr mem_m;
	if (mat->is_in_mem()) {
		dense_matrix::ptr tmp = dense_matrix::create(mat->get_raw_store());
		tmp->materialize_self();
		mem_m = detail::mem_matrix_store::cast(tmp->get_raw_store());
	}
	else {
		dense_matrix::ptr tmp = dense_matrix::create(mat->get_raw_store());
		tmp = tmp->conv_store(true, -1);
		mem_m = detail::mem_matrix_store::cast(tmp->get_raw_store());
		assert(mem_m);
	}
	for (size_t i = 0; i < stdvec.size(); i++) {
		int sum = 0;
		for (size_t j = 0; j < mat->get_num_cols(); j++)
			sum += mem_m->get<int>(i, j);
		assert(sum == stdvec[i]);
	}

	printf("test col sum on %s %s matrix\n", mat->is_wide() ? "wide" : "tall",
			mat->store_layout() == matrix_layout_t::L_COL ? "column" : "row");
	sum = mat->col_sum();
	sum->materialize_self();
	vec = col_vec::create(sum);
	assert(vec->get_length() == mat->get_num_cols());
	stdvec = vec->conv2std<int>();
	for (size_t i = 0; i < stdvec.size(); i++) {
		int sum = 0;
		for (size_t j = 0; j < mat->get_num_rows(); j++)
			sum += mem_m->get<int>(j, i);
		assert(sum == stdvec[i]);
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
			matrix_layout_t::L_NONE);
	dense_matrix::ptr vmat3 = orig_mat1->multiply(*smat5,
			matrix_layout_t::L_NONE);
	dense_matrix::ptr vmat4 = vmat3->multiply(*smat2,
			matrix_layout_t::L_NONE);
	dense_matrix::ptr vmat5 = orig_mat2->multiply(*smat3,
			matrix_layout_t::L_NONE);
	dense_matrix::ptr vmat6 = vmat4->multiply(*smat4,
			matrix_layout_t::L_NONE)->add(*vmat5);
	dense_matrix::ptr res = vmat6->transpose()->multiply(*vmat1,
			matrix_layout_t::L_NONE);
	detail::matrix_stats.print_diff(orig_stats);

	printf("materialize every matrix operations individually\n");
	detail::matrix_stats_t orig_stats1 = detail::matrix_stats;
	dense_matrix::ptr mat1 = orig_mat1->multiply(*smat1,
			matrix_layout_t::L_NONE);
	mat1->materialize_self();
	dense_matrix::ptr mat3 = orig_mat1->multiply(*smat5,
			matrix_layout_t::L_NONE);
	mat3->materialize_self();
	dense_matrix::ptr mat4 = vmat3->multiply(*smat2,
			matrix_layout_t::L_NONE);
	mat4->materialize_self();
	dense_matrix::ptr mat5 = orig_mat2->multiply(*smat3,
			matrix_layout_t::L_NONE);
	mat5->materialize_self();
	dense_matrix::ptr mat6 = vmat4->multiply(*smat4,
			matrix_layout_t::L_NONE)->add(*vmat5);
	mat6->materialize_self();
	dense_matrix::ptr res1 = mat6->transpose()->multiply(*mat1,
			matrix_layout_t::L_NONE);
	detail::matrix_stats.print_diff(orig_stats1);

	if (type == get_scalar_type<int>())
		verify_result(*res, *res1, equal_func<int>());
	else {
		assert(type == get_scalar_type<double>());
		verify_result(*res, *res1, approx_equal_func<double>());
	}
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

	std::vector<off_t> cols1(out[0]->get_num_cols());
	std::vector<off_t> cols2(out[1]->get_num_cols());
	for (size_t i = 0; i < cols1.size(); i++)
		cols1[i] = i;
	for (size_t i = 0; i < cols2.size(); i++)
		cols2[i] = cols1.back() + 1 + i;

	dense_matrix::ptr out_mat0 = dense_matrix::create(out[0]);
	dense_matrix::ptr out_mat1 = dense_matrix::create(out[1]);
	dense_matrix::ptr in_mat = dense_matrix::create(in[0]);
	dense_matrix::ptr in_mat0 = in_mat->get_cols(cols1);
	dense_matrix::ptr in_mat1 = in_mat->get_cols(cols2);
	dense_matrix::ptr diff0 = out_mat0->minus(*in_mat0);
	dense_matrix::ptr diff1 = out_mat1->minus(*in_mat1);
	scalar_variable::ptr agg0 = diff0->sum();
	scalar_variable::ptr agg1 = diff1->sum();
	assert(*(double *) agg0->get_raw() == 0);
	assert(*(double *) agg1->get_raw() == 0);
}

void test_min()
{
	printf("test min\n");
	bulk_operate::const_ptr op = bulk_operate::conv2ptr(
			*get_scalar_type<int>().get_basic_ops().get_op(
			basic_ops::op_idx::MIN));

	dense_matrix::ptr mat = dense_matrix::create_seq<int>(1, 1, 2000, 1,
			matrix_layout_t::L_COL, true);
	scalar_variable::ptr res = mat->aggregate(op);
	assert(*(int *) res->get_raw() == 1);

	mat = dense_matrix::create_seq<int>(1, 1, long_dim, 1,
			matrix_layout_t::L_COL, true);
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
			num_rows, num_cols, get_scalar_type<size_t>()) {
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &_ins,
			detail::local_matrix_store &out) const {
		std::vector<detail::local_matrix_store::const_ptr> ins(_ins.size());
		for (size_t i = 0; i < ins.size(); i++) {
			detail::local_buf_col_matrix_store::ptr in(
					new detail::local_buf_col_matrix_store(
						_ins[i]->get_global_start_row(),
						_ins[i]->get_global_start_col(),
						_ins[i]->get_num_rows(), _ins[i]->get_num_cols(),
						_ins[i]->get_type(), _ins[i]->get_node_id()));
			in->copy_from(*_ins[i]);
			ins[i] = in;
		}

		assert(out.store_layout() == matrix_layout_t::L_COL);
		if (ins.size() == 1) {
			assert(ins[0]->store_layout() == matrix_layout_t::L_COL);
			detail::local_matrix_store::const_ptr in = ins[0]->get_portion(0, 0,
					out.get_num_rows(), out.get_num_cols());
			assert(in);
			out.copy_from(*in);
		}
		else {
			assert(ins.size() >= 2);
			assert(ins[0]->store_layout() == matrix_layout_t::L_COL);
			assert(ins[1]->store_layout() == matrix_layout_t::L_COL);
			detail::local_matrix_store::const_ptr in0 = ins[0]->get_portion(0, 0,
					out.get_num_rows(), out.get_num_cols());
			detail::local_matrix_store::const_ptr in1 = ins[1]->get_portion(0, 0,
					out.get_num_rows(), out.get_num_cols());
			assert(in0 && in1);
			detail::part_dim_t dim = get_out_num_rows() > get_out_num_cols()
				? detail::part_dim_t::PART_DIM1 : detail::part_dim_t::PART_DIM2;
			detail::mapply2(*in0, *in1, out.get_type().get_basic_ops().get_add(),
					dim, out);
			for (size_t i = 2; i < ins.size(); i++) {
				assert(ins[i]->store_layout() == matrix_layout_t::L_COL);
				detail::local_matrix_store::const_ptr in = ins[i]->get_portion(0, 0,
						out.get_num_rows(), out.get_num_cols());
				assert(in);
				detail::mapply2(out, *in,
						out.get_type().get_basic_ops().get_add(), dim, out);
			}
		}
	}

	virtual portion_mapply_op::const_ptr transpose() const {
		return portion_mapply_op::const_ptr();
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return std::string();
	}

	virtual bool is_agg() const {
		return true;
	}
};

class t_add_portion_op: public detail::portion_mapply_op
{
public:
	t_add_portion_op(size_t num_rows, size_t num_cols): detail::portion_mapply_op(
			num_rows, num_cols, get_scalar_type<size_t>()) {
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &_ins,
			detail::local_matrix_store &out) const {
		std::vector<detail::local_matrix_store::const_ptr> ins(_ins.size());
		for (size_t i = 0; i < ins.size(); i++) {
			detail::local_buf_row_matrix_store::ptr in(
					new detail::local_buf_row_matrix_store(
						_ins[i]->get_global_start_row(),
						_ins[i]->get_global_start_col(),
						_ins[i]->get_num_rows(), _ins[i]->get_num_cols(),
						_ins[i]->get_type(), _ins[i]->get_node_id()));
			in->copy_from(*_ins[i]);
			ins[i] = in;
		}

		assert(out.store_layout() == matrix_layout_t::L_ROW);
		if (ins.size() == 1) {
			assert(ins[0]->store_layout() == matrix_layout_t::L_ROW);
			detail::local_matrix_store::const_ptr in = ins[0]->get_portion(0, 0,
					out.get_num_rows(), out.get_num_cols());
			assert(in);
			out.copy_from(*in);
		}
		else {
			assert(ins.size() >= 2);
			assert(ins[0]->store_layout() == matrix_layout_t::L_ROW);
			assert(ins[1]->store_layout() == matrix_layout_t::L_ROW);
			detail::local_matrix_store::const_ptr in0 = ins[0]->get_portion(0, 0,
					out.get_num_rows(), out.get_num_cols());
			detail::local_matrix_store::const_ptr in1 = ins[1]->get_portion(0, 0,
					out.get_num_rows(), out.get_num_cols());
			assert(in0 && in1);
			detail::part_dim_t dim = get_out_num_rows() > get_out_num_cols()
				? detail::part_dim_t::PART_DIM1 : detail::part_dim_t::PART_DIM2;
			detail::mapply2(*in0, *in1, out.get_type().get_basic_ops().get_add(),
					dim, out);
			for (size_t i = 2; i < ins.size(); i++) {
				assert(ins[i]->store_layout() == matrix_layout_t::L_ROW);
				detail::local_matrix_store::const_ptr in = ins[i]->get_portion(0, 0,
						out.get_num_rows(), out.get_num_cols());
				assert(in);
				detail::mapply2(out, *in,
						out.get_type().get_basic_ops().get_add(), dim, out);
			}
		}
	}

	virtual portion_mapply_op::const_ptr transpose() const {
		return portion_mapply_op::const_ptr();
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return std::string();
	}

	virtual bool is_agg() const {
		return true;
	}
};

void test_mapply_mixed(int num_nodes)
{
	if (!safs::is_safs_init())
		return;

	printf("test serial and parallel mapply\n");
	detail::portion_mapply_op::const_ptr op(new add_portion_op(
				long_dim, 10));
	detail::portion_mapply_op::const_ptr t_op(new t_add_portion_op(
				10, long_dim));

	// Test tall and col-major matrices.
	std::vector<dense_matrix::ptr> mats;
	mats.push_back(dense_matrix::create_randu<size_t>(0, 1000,
				long_dim, 10, matrix_layout_t::L_COL, -1, true));
	mats.push_back(dense_matrix::create_randu<size_t>(0, 1000,
				long_dim, 11, matrix_layout_t::L_COL, -1, false));
	mats.push_back(dense_matrix::create_randu<size_t>(0, 1000,
				long_dim, 12, matrix_layout_t::L_COL, -1, false));
	mats.push_back(dense_matrix::create_randu<size_t>(0, 1000,
				long_dim, 13, matrix_layout_t::L_COL, num_nodes, true));
	mats.push_back(dense_matrix::create_randu<size_t>(0, 1000,
				long_dim, 14, matrix_layout_t::L_COL, -1, false));

	std::vector<detail::matrix_store::const_ptr> stores(mats.size());
	for (size_t i = 0; i < stores.size(); i++)
		stores[i] = mats[i]->get_raw_store();
	dense_matrix::ptr par_res = dense_matrix::create(detail::__mapply_portion(
				stores, op, matrix_layout_t::L_COL, true));
	dense_matrix::ptr serial_res = dense_matrix::create(detail::__mapply_portion(
				stores, op, matrix_layout_t::L_COL, false));
	scalar_variable::ptr max_diff = par_res->minus(*serial_res)->abs()->max();
	printf("max diff: %ld\n", *(size_t *) max_diff->get_raw());
	assert(*(size_t *) max_diff->get_raw() == 0);

	// Test wide and row-major matrices.
	for (size_t i = 0; i < mats.size(); i++)
		mats[i] = mats[i]->transpose();
	stores.resize(mats.size());
	for (size_t i = 0; i < stores.size(); i++)
		stores[i] = mats[i]->get_raw_store();
	par_res = dense_matrix::create(detail::__mapply_portion(
				stores, t_op, matrix_layout_t::L_ROW, true));
	serial_res = dense_matrix::create(detail::__mapply_portion(
				stores, t_op, matrix_layout_t::L_ROW, false));
	max_diff = par_res->minus(*serial_res)->abs()->max();
	printf("max diff: %ld\n", *(size_t *) max_diff->get_raw());
	assert(*(size_t *) max_diff->get_raw() == 0);

	// Test tall and row-major matrices.
	mats.clear();
	mats.push_back(dense_matrix::create_randu<size_t>(0, 1000,
				long_dim, 10, matrix_layout_t::L_ROW, -1, true));
	mats.push_back(dense_matrix::create_randu<size_t>(0, 1000,
				long_dim, 11, matrix_layout_t::L_ROW, -1, false));
	mats.push_back(dense_matrix::create_randu<size_t>(0, 1000,
				long_dim, 12, matrix_layout_t::L_ROW, -1, false));
	mats.push_back(dense_matrix::create_randu<size_t>(0, 1000,
				long_dim, 13, matrix_layout_t::L_ROW, num_nodes, true));
	mats.push_back(dense_matrix::create_randu<size_t>(0, 1000,
				long_dim, 14, matrix_layout_t::L_ROW, -1, false));

	stores.resize(mats.size());
	for (size_t i = 0; i < stores.size(); i++)
		stores[i] = mats[i]->get_raw_store();
	par_res = dense_matrix::create(detail::__mapply_portion(
				stores, op, matrix_layout_t::L_COL, true));
	serial_res = dense_matrix::create(detail::__mapply_portion(
				stores, op, matrix_layout_t::L_COL, false));
	max_diff = par_res->minus(*serial_res)->abs()->max();
	printf("max diff: %ld\n", *(size_t *) max_diff->get_raw());
	assert(*(size_t *) max_diff->get_raw() == 0);

	// Test wide and col-major matrices.
	for (size_t i = 0; i < mats.size(); i++)
		mats[i] = mats[i]->transpose();
	stores.resize(mats.size());
	for (size_t i = 0; i < stores.size(); i++)
		stores[i] = mats[i]->get_raw_store();
	par_res = dense_matrix::create(detail::__mapply_portion(
				stores, t_op, matrix_layout_t::L_ROW, true));
	serial_res = dense_matrix::create(detail::__mapply_portion(
				stores, t_op, matrix_layout_t::L_ROW, false));
	max_diff = par_res->minus(*serial_res)->abs()->max();
	printf("max diff: %ld\n", *(size_t *) max_diff->get_raw());
	assert(*(size_t *) max_diff->get_raw() == 0);
}

namespace fm
{
namespace eigen
{
extern size_t num_cached_mats;
}
}

void test_sub_matrix()
{
	printf("test sub tall col-matrix\n");
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
			long_dim, num_cols, 2, get_scalar_type<double>(), in_mem, false);
	mv1->set_block(0, sub_mat1);
	mv1->set_block(1, sub_mat2);
	mv1->set_block(2, mat3);

	// Create the block MV from the in-mem matrices.
	eigen::block_multi_vector::ptr mv2 = eigen::block_multi_vector::create(
			long_dim, num_cols, 2, get_scalar_type<double>(), in_mem, false);
	mv2->set_block(0, copy1);
	mv2->set_block(1, copy2);
	mv2->set_block(2, mat3);

	// Create the right matrix for GEMM.
	dense_matrix::ptr right_mat = dense_matrix::create_randu<double>(0, 1,
			num_cols, 2, matrix_layout_t::L_COL, -1, true);
	detail::mem_col_matrix_store::const_ptr B
		= detail::mem_col_matrix_store::cast(right_mat->get_raw_store());
	scalar_variable_impl<double> alpha(2);
	scalar_variable_impl<double> beta(3);

	// Perform GEMM on the first block MV.
	eigen::block_multi_vector::ptr res1;
	dense_matrix::ptr res_mat1;
	fm::eigen::num_cached_mats = 0;
	{
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		res1 = eigen::block_multi_vector::create(
				long_dim, mv1->get_block_size(), mv1->get_block_size(),
				get_scalar_type<double>(), false, false);
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
				get_scalar_type<double>(), true, false);
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

void test_setdata(int num_nodes)
{
	dense_matrix::ptr mat = create_seq_matrix(long_dim, 10,
		matrix_layout_t::L_ROW, num_nodes, get_scalar_type<int>(), in_mem);
	size_t num_portions = mat->get_data().get_num_portions();
	for (size_t k = 0; k < num_portions; k++) {
		detail::local_matrix_store::const_ptr part
			= mat->get_data().get_portion(k);
		for (size_t i = 0; i < part->get_num_rows(); i++)
			for (size_t j = 0; j < part->get_num_cols(); j++) {
				size_t row_idx = part->get_global_start_row() + i;
				size_t col_idx = j;
				int expected = row_idx * mat->get_num_cols() + col_idx;
				assert(part->get<int>(i, j) == expected);
			}
	}
}

void _test_groupby_eles();
void _test_groupby(dense_matrix::ptr mat);

void test_groupby()
{
	dense_matrix::ptr mat;
	_test_groupby_eles();

	mat = create_matrix(10, long_dim, matrix_layout_t::L_ROW,
			-1, get_scalar_type<int>());
	mat->materialize_self();
	assert(!mat->is_virtual());
	_test_groupby(mat);

	mat = mat->abs();
	assert(mat->is_virtual());
	_test_groupby(mat);

	mat = create_matrix(long_dim, 10, matrix_layout_t::L_ROW,
			-1, get_scalar_type<int>());
	mat->materialize_self();
	assert(!mat->is_virtual());
	_test_groupby(mat);

	mat = mat->abs();
	assert(mat->is_virtual());
	_test_groupby(mat);
}

void verify_cum(dense_matrix::ptr res, dense_matrix::ptr mat, bool byrow)
{
	assert(res->get_num_rows() == mat->get_num_rows());
	assert(res->get_num_cols() == mat->get_num_cols());
	assert(res->get_type() == mat->get_type());
	res = dense_matrix::create(res->get_raw_store());
	mat = dense_matrix::create(mat->get_raw_store());
	res = res->conv_store(true, -1);
	mat = mat->conv_store(true, -1);
	assert(mat->is_in_mem());
	assert(res->is_in_mem());

	detail::mem_matrix_store::const_ptr mem_mat
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				mat->get_raw_store());
	detail::mem_matrix_store::const_ptr mem_res
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				res->get_raw_store());
	assert(mem_mat);
	assert(mem_res);
	if (byrow) {
		for (size_t i = 0; i < res->get_num_rows(); i++) {
			assert(mem_res->get<int>(i, 0) == mem_mat->get<int>(i, 0));
			for (size_t j = 1; j < res->get_num_cols(); j++)
				assert(mem_res->get<int>(i, j) - mem_res->get<int>(i, j - 1)
						== mem_mat->get<int>(i, j));
		}
	}
	else {
		for (size_t j = 0; j < res->get_num_cols(); j++) {
			assert(mem_res->get<int>(0, j) == mem_mat->get<int>(0, j));
			for (size_t i = 1; i < res->get_num_rows(); i++)
				assert(mem_res->get<int>(i, j) - mem_res->get<int>(i - 1, j)
						== mem_mat->get<int>(i, j));
		}
	}
}

void test_cum(int num_nodes)
{
	dense_matrix::ptr mat1 = create_matrix(long_dim, 10,
			matrix_layout_t::L_COL, num_nodes, get_scalar_type<int>());
	auto op = mat1->get_type().get_agg_ops().get_op(agg_ops::op_idx::SUM);
	printf("test cum on rows\n");
	dense_matrix::ptr res = mat1->cum(matrix_margin::MAR_ROW, op);
	verify_cum(res, mat1, true);
#if 0
	res = mat1->cum(matrix_margin::MAR_COL, op);
	verify_cum(res, mat1, false);
#endif

	printf("test cum on cols\n");
	mat1 = mat1->transpose();
	res = mat1->cum(matrix_margin::MAR_ROW, op);
	verify_cum(res, mat1, true);
#if 0
	res = mat1->cum(matrix_margin::MAR_COL, op);
	verify_cum(res, mat1, false);
#endif
}

void _test_EM_matrix()
{
	if (!safs::is_safs_init())
		return;

	printf("test EM matrix\n");
	in_mem = false;

	test_cum(-1);
	test_groupby();
	test_agg(-1, matrix_layout_t::L_ROW);
	test_agg(-1, matrix_layout_t::L_COL);
	test_setdata(-1);
	test_sub_matrix();
	test_mapply_chain(-1, get_scalar_type<double>());
	test_mapply_chain(-1, get_scalar_type<int>());
	test_multiply<double>(-1);
	test_multiply<float>(-1);
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
	test_multiply_matrix(-1);
	test_agg_sub_col(-1);
	test_agg_sub_row(-1);
	test_sum_row_col(-1);
#if 0
	test_rand_init();
	test_conv_row_col();
	test_flatten();
#endif
}

void test_EM_matrix()
{
	block_size = 3;
	matrix_val = matrix_val_t::SEQ;
	_test_EM_matrix();
	block_size = 0;
	_test_EM_matrix();
}

void _test_mem_matrix(int num_nodes)
{
	printf("test mem matrix\n");
	in_mem = true;

	test_cum(-1);
	test_cum(num_nodes);
	test_groupby();
	test_agg(-1, matrix_layout_t::L_COL);
	test_agg(num_nodes, matrix_layout_t::L_COL);
	test_agg(-1, matrix_layout_t::L_ROW);
	test_agg(num_nodes, matrix_layout_t::L_ROW);
	test_setdata(-1);
	test_setdata(num_nodes);
	test_copy(-1, true);
	test_apply_scalar();
	test_min();
	test_mul_output(-1);
	test_mul_output(num_nodes);
	test_mapply_chain(-1, get_scalar_type<double>());
	test_mapply_chain(-1, get_scalar_type<int>());
	test_mapply_chain(num_nodes, get_scalar_type<int>());
	test_multiply<double>(-1);
	test_multiply<float>(-1);
	test_cast();
	test_write2file();
	test_apply();
	test_conv_vec2mat();

	test_sum_row_col(-1);
	test_sum_row_col(num_nodes);
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
	test_multiply_matrix(-1);
	test_multiply_matrix(num_nodes);
	test_agg_sub_col(-1);
	test_agg_sub_row(-1);
#if 0
	test_rand_init();
	test_conv_row_col();
	test_flatten();
#endif
}

void test_mem_matrix(int num_nodes)
{
	auto orig = matrix_val;
	matrix_val = matrix_val_t::SPARSE;
	_test_mem_matrix(num_nodes);
	block_size = 3;
	matrix_val = matrix_val_t::SEQ_MATER;
	_test_mem_matrix(num_nodes);
	matrix_val = matrix_val_t::SEQ;
	_test_mem_matrix(num_nodes);
	block_size = 0;
	_test_mem_matrix(num_nodes);
	matrix_val = matrix_val_t::DEFAULT;
	_test_mem_matrix(num_nodes);
	matrix_val = orig;
}

void test_conv_store()
{
	if (!safs::is_safs_init())
		return;

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

class rand_set: public type_set_operate<factor_value_t>
{
	factor f;
public:
	rand_set(const factor &_f): f(_f) {
	}

	void set(factor_value_t *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = random() % f.get_num_levels();
			assert(f.is_valid_level(arr[i]));
		}
	}
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

void _test_groupby_eles()
{
	dense_matrix::ptr mat = dense_matrix::create_randu<int>(0, 1000, long_dim,
			10, matrix_layout_t::L_COL);
	printf("test groupby elements\n");
	agg_operate::const_ptr count_agg = mat->get_type().get_agg_ops().get_count();
	data_frame::ptr res = mat->groupby(count_agg, true);
	printf("#entries in res: %ld\n", res->get_num_entries());

	std::map<int, size_t> ele_counts;
	const detail::mem_matrix_store &vstore
		= dynamic_cast<const detail::mem_matrix_store &>(mat->get_data());
	for (size_t i = 0; i < mat->get_num_rows(); i++)
		for (size_t j = 0; j < mat->get_num_cols(); j++) {
			int val = vstore.get<int>(i, j);
			auto it = ele_counts.find(val);
			if (it == ele_counts.end())
				ele_counts.insert(std::pair<int, size_t>(val, 1));
			else
				it->second++;
		}

	detail::smp_vec_store::ptr vals = detail::smp_vec_store::cast(
			res->get_vec("val"));
	detail::smp_vec_store::ptr aggs = detail::smp_vec_store::cast(
			res->get_vec("agg"));
	assert(vals->get_length() == aggs->get_length());
	for (size_t i = 0; i < vals->get_length(); i++) {
		int val = vals->get<int>(i);
		size_t count = aggs->get<size_t>(i);
		auto it = ele_counts.find(val);
		assert(it != ele_counts.end());
		assert(it->second == count);
	}
}

void _test_groupby(dense_matrix::ptr mat)
{
	printf("group by rows on matrix: %ld, %ld\n", mat->get_num_rows(),
			mat->get_num_cols());
	factor f(10);
	factor_col_vector::ptr rand_factors = factor_col_vector::create(f,
			mat->get_num_rows(), -1, mat->is_in_mem(), rand_set(f));
	bulk_operate::const_ptr add = bulk_operate::conv2ptr(
			mat->get_type().get_basic_ops().get_add());
	dense_matrix::ptr group_sum = mat->groupby_row(rand_factors, add);
	group_sum->materialize_self();

	dense_matrix::ptr tmp = dense_matrix::create(mat->get_raw_store());
	tmp = tmp->conv_store(true, -1);
	detail::mem_matrix_store::const_ptr mem_mat
		= detail::mem_matrix_store::cast(tmp->get_raw_store());
	dense_matrix::ptr tmp_factors = rand_factors->conv_store(true, -1);
	detail::mem_matrix_store::const_ptr mem_factors
		= detail::mem_matrix_store::cast(tmp_factors->get_raw_store());
	std::vector<std::vector<int> > agg_res(f.get_num_levels());
	for (size_t i = 0; i < f.get_num_levels(); i++)
		agg_res[i].resize(mat->get_num_cols());
	printf("compute groupby manually\n");
	for (size_t i = 0; i < mat->get_num_rows(); i++) {
		factor_value_t label = mem_factors->get<factor_value_t>(i, 0);
		assert(agg_res.size() > label);
		assert(agg_res[label].size() == mat->get_num_cols());
		for (size_t j = 0; j < mat->get_num_cols(); j++)
			agg_res[label][j] += mem_mat->get<factor_value_t>(i, j);
	}

	printf("check groupby\n");
	group_sum = group_sum->conv_store(true, -1);
	detail::mem_matrix_store::const_ptr mem_res
		= detail::mem_matrix_store::cast(group_sum->get_raw_store());
	assert(group_sum->get_num_rows() == agg_res.size());
	for (size_t i = 0; i < mem_res->get_num_rows(); i++) {
		assert(agg_res[i].size() == group_sum->get_num_cols());
		for (size_t j = 0; j < mem_res->get_num_cols(); j++) {
			assert(agg_res[i][j] == mem_res->get<int>(i, j));
		}
	}

	dense_matrix::ptr group_sum1 = mat->groupby_row(rand_factors, add);
	group_sum1 = group_sum1->transpose();
	group_sum1->materialize_self();

	dense_matrix::ptr group_sum_t = group_sum->transpose();
	dense_matrix::ptr diff = group_sum_t->minus(*group_sum1);
	scalar_variable::ptr sum = diff->sum();
	assert(scalar_variable::get_val<int>(*sum) == 0);

	// Test mapply on groupby matrix.
	dense_matrix::ptr sink = mat->groupby_row(rand_factors, add);
	if (sink->get_raw_store()->is_sink()) {
		test_mapply_sink(sink);
		sink = mat->groupby_row(rand_factors, add);
		test_mapply_sink_t(sink);
	}
}

dense_matrix::ptr _test_get_rows(dense_matrix::ptr mat, size_t get_nrow)
{
	std::vector<off_t> idxs;
	idxs.resize(get_nrow);
	for (size_t i = 0; i < idxs.size(); i++)
		idxs[i] = random() % mat->get_num_rows();
	dense_matrix::ptr res = mat->get_rows(idxs);
	assert(res != NULL);
	mat = dense_matrix::create(mat->get_raw_store());
	dense_matrix::ptr mem_mat = mat->conv_store(true, -1);
	res = dense_matrix::create(res->get_raw_store());
	dense_matrix::ptr mem_res = res->conv_store(true, -1);
	detail::mem_matrix_store::const_ptr orig_store
		= detail::mem_matrix_store::cast(mem_mat->get_raw_store());
	detail::mem_matrix_store::const_ptr res_store
		= detail::mem_matrix_store::cast(mem_res->get_raw_store());
	for (size_t i = 0; i < res->get_num_rows(); i++)
		for (size_t j = 0; j < res->get_num_cols(); j++)
			assert(res_store->get<int>(i, j) == orig_store->get<int>(idxs[i], j));
	return res;
}

dense_matrix::ptr _test_get_rows_bool(dense_matrix::ptr mat)
{
	printf("get rows in bool idx\n");
	dense_matrix::ptr res;
	col_vec::ptr bool_idx;
	while (res == NULL) {
		bool_idx = col_vec::create_randu<bool>(mat->get_num_rows());
		res = mat->get_rows(bool_idx);
	}
	assert(res->get_num_cols() == mat->get_num_cols());

	std::vector<bool> bools = bool_idx->conv2std<bool>();
	size_t num_trues = 0;
	for (size_t i = 0; i < bools.size(); i++)
		if (bools[i])
			num_trues++;
	printf("want %ld rows and get %ld rows\n", num_trues, res->get_num_rows());
	assert(res->get_num_rows() == num_trues);

	// Make sure it's row-oriented.
	mat = mat->conv2(matrix_layout_t::L_ROW);
	mat = mat->conv_store(true, -1);
	mat->materialize_self();
	res = res->conv2(matrix_layout_t::L_ROW);
	res = res->conv_store(true, -1);

	detail::mem_matrix_store::const_ptr store
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				mat->get_raw_store());
	detail::mem_matrix_store::const_ptr res_store
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				res->get_raw_store());
	size_t res_i = 0;
	for (size_t i = 0; i < store->get_num_rows(); i++)
		if (bools[i]) {
			assert(memcmp(res_store->get_row(res_i), store->get_row(i),
					store->get_num_cols() * store->get_entry_size()) == 0);
			res_i++;
		}
	return res;
}

dense_matrix::ptr test_get_rows(dense_matrix::ptr mat)
{
	_test_get_rows_bool(mat);
	return _test_get_rows(mat, std::max(mat->get_num_rows() / 5, 2UL));
}

dense_matrix::ptr _test_get_cols(dense_matrix::ptr mat, size_t get_ncol)
{
	std::vector<off_t> idxs;
	idxs.resize(get_ncol);
	for (size_t i = 0; i < idxs.size(); i++)
		idxs[i] = (random() + mat->get_num_cols() / 2) % mat->get_num_cols();
	dense_matrix::ptr res = mat->get_cols(idxs);
	assert(res != NULL);
	res = dense_matrix::create(res->get_raw_store());
	dense_matrix::ptr mem_res = res->conv_store(true, -1);
	mat = dense_matrix::create(mat->get_raw_store());
	dense_matrix::ptr mem_mat = mat->conv_store(true, -1);
	detail::mem_matrix_store::const_ptr orig_store
		= detail::mem_matrix_store::cast(mem_mat->get_raw_store());
	assert(orig_store);
	detail::mem_matrix_store::const_ptr res_store
		= detail::mem_matrix_store::cast(mem_res->get_raw_store());
	assert(res_store);
	for (size_t i = 0; i < res->get_num_rows(); i++)
		for (size_t j = 0; j < res->get_num_cols(); j++)
			assert(res_store->get<int>(i, j) == orig_store->get<int>(i, idxs[j]));
	return res;
}

dense_matrix::ptr _test_get_cols_bool(dense_matrix::ptr mat)
{
	printf("get cols in bool idx, mat: %ld,%ld\n", mat->get_num_rows(),
			mat->get_num_cols());
	dense_matrix::ptr res;
	col_vec::ptr bool_idx;
	while (res == NULL) {
		bool_idx = col_vec::create_randu<bool>(mat->get_num_cols());
		res = mat->get_cols(bool_idx);
	}
	assert(res->get_num_rows() == mat->get_num_rows());

	std::vector<bool> bools = bool_idx->conv2std<bool>();
	size_t num_trues = 0;
	for (size_t i = 0; i < bools.size(); i++)
		if (bools[i])
			num_trues++;
	assert(res->get_num_cols() == num_trues);

	mat = mat->conv2(matrix_layout_t::L_COL);
	mat = mat->conv_store(true, -1);
	mat->materialize_self();
	res = res->conv2(matrix_layout_t::L_COL);
	res = res->conv_store(true, -1);

	detail::mem_matrix_store::const_ptr store
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				mat->get_raw_store());
	detail::mem_matrix_store::const_ptr res_store
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				res->get_raw_store());
	size_t res_i = 0;
	for (size_t i = 0; i < store->get_num_cols(); i++)
		if (bools[i]) {
			assert(memcmp(res_store->get_col(res_i), store->get_col(i),
					store->get_num_rows() * store->get_entry_size()) == 0);
			res_i++;
		}
	return res;
}

dense_matrix::ptr test_get_cols(dense_matrix::ptr mat)
{
	_test_get_cols_bool(mat);
	return _test_get_cols(mat, std::max(mat->get_num_cols() / 5, 2UL));
}

void _test_get_rowcols(dense_matrix::ptr mat)
{
	if (mat->is_wide())
		_test_get_cols(mat, 5);
	else
		_test_get_rows(mat, 5);
	dense_matrix::ptr tmp = test_get_rows(mat);
	test_get_cols(tmp);
	test_get_rows(tmp);
	tmp = test_get_cols(mat);
	test_get_cols(tmp);
	test_get_rows(tmp);
}

void test_get_rowcols(int num_nodes)
{
	block_size = 3;
	dense_matrix::ptr mat, tmp;

	if (safs::is_safs_init()) {
		bool orig_in_mem = in_mem;
		in_mem = false;

		printf("test on col-major tall dense matrix in disks\n");
		mat = create_matrix(long_dim, 10, matrix_layout_t::L_COL, -1,
				get_scalar_type<int>());
		tmp = test_get_rows(mat);
		assert(tmp->is_in_mem());
		tmp = test_get_cols(mat);
		assert(!tmp->is_in_mem());
		test_get_cols(tmp);
		test_get_rows(tmp);
		tmp = tmp->add(*tmp);
		test_get_cols(tmp);
		test_get_rows(tmp);

		printf("test on row-major tall dense matrix in disks\n");
		mat = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, -1,
				get_scalar_type<int>());
		tmp = test_get_rows(mat);
		assert(tmp->is_in_mem());
		test_get_cols(mat);

		printf("test on row-major wide dense matrix in disks\n");
		mat = create_matrix(10, long_dim, matrix_layout_t::L_ROW, -1,
				get_scalar_type<int>());
		tmp = test_get_cols(mat);
		assert(tmp->is_in_mem());
		tmp = test_get_rows(mat);
		test_get_cols(tmp);
		test_get_rows(tmp);
		tmp = tmp->add(*tmp);
		test_get_cols(tmp);
		test_get_rows(tmp);

		in_mem = orig_in_mem;
	}

	printf("test on a row-major tall matrix in SMP\n");
	mat = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, -1,
			get_scalar_type<int>());
	_test_get_rowcols(mat);

	printf("test on a col-major tall matrix in SMP\n");
	mat = create_matrix(long_dim, 10, matrix_layout_t::L_COL, -1,
			get_scalar_type<int>());
	_test_get_rowcols(mat);

	printf("test on a row-major tall matrix in NUMA\n");
	mat = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<int>());
	_test_get_rowcols(mat);

	printf("test on a col-major tall matrix in NUMA\n");
	mat = create_matrix(long_dim, 10, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<int>());
	_test_get_rowcols(mat);

	printf("test on a row-major wide matrix in NUMA\n");
	mat = create_matrix(10, long_dim, matrix_layout_t::L_ROW, num_nodes,
			get_scalar_type<int>());
	_test_get_rowcols(mat);

	printf("test on a col-major wide matrix in NUMA\n");
	mat = create_matrix(10, long_dim, matrix_layout_t::L_COL, num_nodes,
			get_scalar_type<int>());
	_test_get_rowcols(mat);

	printf("test on a row-major tall virtual matrix in SMP\n");
	mat = create_matrix(long_dim, 10, matrix_layout_t::L_ROW, -1,
			get_scalar_type<int>());
	mat = mat->add(*mat);
	_test_get_rowcols(mat);

	printf("test on a col-major tall virtual matrix in SMP\n");
	mat = create_matrix(long_dim, 10, matrix_layout_t::L_COL, -1,
			get_scalar_type<int>());
	mat = mat->add(*mat);
	_test_get_rowcols(mat);

	block_size = 0;
}

void test_set_cols_tall()
{
	printf("test set cols in a tall matrix\n");
	detail::mem_col_matrix_store::ptr data = detail::mem_col_matrix_store::create(
			long_dim, 32, get_scalar_type<double>());
	detail::mem_col_matrix_store::ptr res = detail::mem_col_matrix_store::create(
			long_dim, 32, get_scalar_type<double>());
	detail::mem_col_matrix_store::ptr new_data = detail::mem_col_matrix_store::create(
			long_dim, 10, get_scalar_type<double>());
	data->init_randu<double>(0, 1);
	new_data->init_randu<double>(0, 1);

	// Set cols in a contiguous range in a tall matrix.
	std::vector<off_t> idxs(new_data->get_num_cols());
	for (size_t i = 0; i < idxs.size(); i++)
		idxs[i] = 5 + i;
	memcpy(res->get_raw_arr(), data->get_raw_arr(),
			data->get_num_rows() * data->get_num_cols() * data->get_entry_size());
	for (size_t i = 0; i < idxs.size(); i++)
		memcpy(res->get_col(idxs[i]), new_data->get_col(i),
				data->get_num_rows() * data->get_entry_size());

	dense_matrix::ptr res_mat = dense_matrix::create(res);
	dense_matrix::ptr m1 = dense_matrix::create(data);
	dense_matrix::ptr m2 = dense_matrix::create(new_data);
	dense_matrix::ptr res1 = m1->set_cols(idxs, m2);
	verify_result(*res1, *res_mat, equal_func<double>());

	// Set cols in a range with gaps in a tall matrix.
	memcpy(res->get_raw_arr(), data->get_raw_arr(),
			data->get_num_rows() * data->get_num_cols() * data->get_entry_size());
	for (size_t i = 2; i < 22; i += 2)
		memcpy(res->get_col(i), new_data->get_col((i - 2)/2),
				data->get_num_rows() * data->get_entry_size());

	res_mat = dense_matrix::create(res);
	m1 = dense_matrix::create(data);
	m2 = dense_matrix::create(new_data);
	res1 = m1->set_cols(2, 22, 2, m2);
	verify_result(*res1, *res_mat, equal_func<double>());

	// Set individual cols in a sorted order in a tall matrix.
	idxs[0] = 3;
	idxs[1] = 4;
	idxs[2] = 5;
	idxs[3] = 10;
	idxs[4] = 11;
	idxs[5] = 12;
	idxs[6] = 13;
	idxs[7] = 14;
	idxs[8] = 20;
	idxs[9] = 21;
	memcpy(res->get_raw_arr(), data->get_raw_arr(),
			data->get_num_rows() * data->get_num_cols() * data->get_entry_size());
	for (size_t i = 0; i < idxs.size(); i++)
		memcpy(res->get_col(idxs[i]), new_data->get_col(i),
				data->get_num_rows() * data->get_entry_size());

	res_mat = dense_matrix::create(res);
	m1 = dense_matrix::create(data);
	m2 = dense_matrix::create(new_data);
	res1 = m1->set_cols(idxs, m2);
	verify_result(*res1, *res_mat, equal_func<double>());

	// Set cols in a random order in a tall matrix.
	for (size_t i = 0; i < idxs.size(); i++)
		idxs[i] = random() % data->get_num_cols();
	memcpy(res->get_raw_arr(), data->get_raw_arr(),
			data->get_num_rows() * data->get_num_cols() * data->get_entry_size());
	for (size_t i = 0; i < idxs.size(); i++)
		memcpy(res->get_col(idxs[i]), new_data->get_col(i),
				data->get_num_rows() * data->get_entry_size());

	res_mat = dense_matrix::create(res);
	m1 = dense_matrix::create(data);
	m2 = dense_matrix::create(new_data);
	res1 = m1->set_cols(idxs, m2);
	verify_result(*res1, *res_mat, equal_func<double>());
}

void test_set_rows_tall()
{
	printf("test set rows in a tall matrix\n");
	detail::mem_row_matrix_store::ptr data = detail::mem_row_matrix_store::create(
			long_dim, 32, get_scalar_type<double>());
	detail::mem_row_matrix_store::ptr res = detail::mem_row_matrix_store::create(
			long_dim, 32, get_scalar_type<double>());
	detail::mem_row_matrix_store::ptr new_data = detail::mem_row_matrix_store::create(
			long_dim / 100, 32, get_scalar_type<double>());
	data->init_randu<double>(0, 1);
	new_data->init_randu<double>(0, 1);

	// Set rows in a sorted order in a tall matrix.
	std::vector<off_t> idxs(new_data->get_num_rows());
	for (size_t i = 0; i < idxs.size(); i++)
		idxs[i] = random() % data->get_num_rows();
	std::sort(idxs.begin(), idxs.end());
	memcpy(res->get_raw_arr(), data->get_raw_arr(),
			data->get_num_rows() * data->get_num_cols() * data->get_entry_size());
	for (size_t i = 0; i < idxs.size(); i++)
		memcpy(res->get_row(idxs[i]), new_data->get_row(i),
				data->get_num_cols() * data->get_entry_size());

	dense_matrix::ptr res_mat = dense_matrix::create(res);
	dense_matrix::ptr m1 = dense_matrix::create(data);
	dense_matrix::ptr m2 = dense_matrix::create(new_data);
	dense_matrix::ptr res1 = m1->set_rows(idxs, m2);
	verify_result(*res1, *res_mat, equal_func<double>());

	// Set rows in a random order in a tall matrix.
	for (size_t i = 0; i < idxs.size(); i++)
		idxs[i] = random() % data->get_num_rows();
	memcpy(res->get_raw_arr(), data->get_raw_arr(),
			data->get_num_rows() * data->get_num_cols() * data->get_entry_size());
	for (size_t i = 0; i < idxs.size(); i++)
		memcpy(res->get_row(idxs[i]), new_data->get_row(i),
				data->get_num_cols() * data->get_entry_size());

	res_mat = dense_matrix::create(res);
	m1 = dense_matrix::create(data);
	m2 = dense_matrix::create(new_data);
	res1 = m1->set_rows(idxs, m2);
	verify_result(*res1, *res_mat, equal_func<double>());

	// Set rows in a tall col matrix
	res_mat = dense_matrix::create(res);
	m1 = dense_matrix::create(data);
	m1 = m1->conv2(matrix_layout_t::L_COL);
	m2 = dense_matrix::create(new_data);
	m2 = m2->conv2(matrix_layout_t::L_COL);
	res1 = m1->set_rows(idxs, m2);
	verify_result(*res1, *res_mat, equal_func<double>());

	// Set cols in a random order in a wide matrix.
	res_mat = dense_matrix::create(res->transpose());
	m1 = dense_matrix::create(data->transpose());
	m2 = dense_matrix::create(new_data->transpose());
	res1 = m1->set_cols(idxs, m2);
	verify_result(*res1, *res_mat, equal_func<double>());

	// Set rows in a tall col matrix
	res_mat = dense_matrix::create(res->transpose());
	m1 = dense_matrix::create(data->transpose());
	m1 = m1->conv2(matrix_layout_t::L_ROW);
	m2 = dense_matrix::create(new_data->transpose());
	m2 = m2->conv2(matrix_layout_t::L_ROW);
	res1 = m1->set_cols(idxs, m2);
	verify_result(*res1, *res_mat, equal_func<double>());
}

void test_set_eles()
{
	detail::mem_matrix_store::ptr data = detail::mem_row_matrix_store::create(
			1000, 32, get_scalar_type<double>());
	detail::mem_matrix_store::ptr col_idx = detail::mem_col_matrix_store::create(
			1000, 1, get_scalar_type<off_t>());
	detail::mem_matrix_store::ptr new_data = detail::mem_row_matrix_store::create(
			1000, 1, get_scalar_type<double>());
	detail::mem_matrix_store::ptr res = detail::mem_row_matrix_store::create(
			1000, 32, get_scalar_type<double>());
	data->init_randu<double>(0, 1);
	col_idx->init_randu<off_t>(0, data->get_num_cols() - 1);
	new_data->init_randu<double>(0, 1);

	memcpy(res->get_raw_arr(), data->get_raw_arr(),
			data->get_num_rows() * data->get_num_cols() * data->get_entry_size());
	for (size_t i = 0; i < data->get_num_rows(); i++) {
		off_t idx = *(off_t *) col_idx->get(i, 0);
		memcpy(res->get(i, idx), new_data->get(i, 0), res->get_type().get_size());
	}

	printf("test set elements on a tall matrix\n");
	std::vector<dense_matrix::ptr> idx_vec(2);
	idx_vec[0] = dense_matrix::create_seq<off_t>(0, 1,
			data->get_num_rows(), 1, matrix_layout_t::L_COL, false);
	idx_vec[1] = dense_matrix::create(col_idx);
	dense_matrix::ptr idxs = dense_matrix::cbind(idx_vec);
	dense_matrix::ptr res_mat = dense_matrix::create(res);
	dense_matrix::ptr m1 = dense_matrix::create(data);
	col_vec::ptr m2 = col_vec::create(new_data);
	dense_matrix::ptr res1 = m1->set_eles(idxs, m2);
	verify_result(*res1, *res_mat, equal_func<double>());

	printf("test set elements on a wide matrix\n");
	res1 = m1->set_eles(idxs, m2);
	res1 = res1->transpose();
	verify_result(*res1, *res_mat->transpose(), equal_func<double>());
}

void test_set_rowcols()
{
	test_set_eles();
	test_set_cols_tall();
	test_set_rows_tall();
}

void _test_repeat_rowcols(dense_matrix::ptr mat, size_t long_dim)
{
	col_vec::ptr idxs = col_vec::create(dense_matrix::create_randu<size_t>(0,
				mat->get_num_rows() - 1, long_dim, 1, matrix_layout_t::L_COL));
	dense_matrix::ptr tmp = mat->get_rows(idxs);
	assert(tmp->get_type() == mat->get_type());
	assert(tmp->get_num_rows() == long_dim);
	assert(tmp->get_num_cols() == mat->get_num_cols());
	dense_matrix::ptr rsum = tmp->row_sum();
	assert(rsum->get_type() == get_scalar_type<int>());
	dense_matrix::ptr scale_idxs
		= idxs->multiply_scalar<size_t>(tmp->get_num_cols());
	scalar_variable::ptr diff_sum = rsum->minus(*scale_idxs)->abs()->sum();
	assert(diff_sum->get_type() == get_scalar_type<size_t>());
	assert(scalar_variable::get_val<size_t>(*diff_sum) == 0);

	tmp = tmp->transpose();
	dense_matrix::ptr csum = tmp->col_sum();
	csum = csum->transpose();
	scale_idxs = scale_idxs->transpose();
	diff_sum = csum->minus(*scale_idxs)->abs()->sum();
	assert(diff_sum->get_type() == get_scalar_type<size_t>());
	assert(scalar_variable::get_val<size_t>(*diff_sum) == 0);
}

void test_repeat_rowcols()
{
	col_vec::ptr seq;
	dense_matrix::ptr mat;

	seq = col_vec::create(dense_matrix::create_seq<int>(0, 1,
				10, 1, matrix_layout_t::L_COL, false));
	mat = dense_matrix::create_repeat(seq, seq->get_length(), 10,
			matrix_layout_t::L_COL, false);
	assert(mat->get_num_rows() == 10);
	assert(mat->get_num_cols() == 10);
	_test_repeat_rowcols(mat, 1000);
	_test_repeat_rowcols(mat, long_dim);

	mat = dense_matrix::create_repeat(seq, seq->get_length(), 100,
			matrix_layout_t::L_COL, false);
	assert(mat->get_num_rows() == 10);
	assert(mat->get_num_cols() == 100);
	_test_repeat_rowcols(mat, 1000);
	_test_repeat_rowcols(mat, long_dim);

	seq = col_vec::create(dense_matrix::create_seq<int>(0, 1,
				100000, 1, matrix_layout_t::L_COL, false));
	mat = dense_matrix::create_repeat(seq, seq->get_length(), 100,
			matrix_layout_t::L_COL, false);
	assert(mat->get_num_rows() == 100000);
	assert(mat->get_num_cols() == 100);
	_test_repeat_rowcols(mat, 1000);
	_test_repeat_rowcols(mat, long_dim);
}

void test_materialize(int num_nodes)
{
	// Test in-memory tall matrix
	printf("Test full materialization on in-mem tall matrix\n");
	matrix_val_t matrix_val = matrix_val_t::SEQ;
	dense_matrix::ptr m1 = create_matrix(long_dim, 1, matrix_layout_t::L_COL,
			num_nodes, get_scalar_type<int>());
	dense_matrix::ptr tmp = m1->add(*m1);
	tmp = tmp->cast_ele_type(get_scalar_type<size_t>());
	tmp->set_materialize_level(materialize_level::MATER_FULL);
	scalar_variable::ptr res = tmp->sum();
	assert(*(size_t *) res->get_raw() == ((long_dim - 1) * long_dim));

	detail::mapply_matrix_store::const_ptr vstore
		= std::dynamic_pointer_cast<const detail::mapply_matrix_store>(
				tmp->get_raw_store());
	assert(vstore);
	assert(vstore->has_materialized());
	dense_matrix::ptr materialize_res = dense_matrix::create(vstore->materialize(
			vstore->is_in_mem(), vstore->get_num_nodes()));
	res = materialize_res->sum();
	assert(*(size_t *) res->get_raw() == ((long_dim - 1) * long_dim));

	// Test in-memory wide matrix
	printf("Test full materialization on in-mem wide matrix\n");
	m1 = create_matrix(1, long_dim, matrix_layout_t::L_COL,
			num_nodes, get_scalar_type<int>());
	tmp = m1->add(*m1);
	tmp = tmp->cast_ele_type(get_scalar_type<size_t>());
	tmp->set_materialize_level(materialize_level::MATER_FULL);
	res = tmp->sum();
	assert(*(size_t *) res->get_raw() == ((long_dim - 1) * long_dim));

	vstore = std::dynamic_pointer_cast<const detail::mapply_matrix_store>(
				tmp->get_raw_store());
	assert(vstore);
	assert(vstore->has_materialized());
	materialize_res = dense_matrix::create(vstore->materialize(
			vstore->is_in_mem(), vstore->get_num_nodes()));
	res = materialize_res->sum();
	assert(*(size_t *) res->get_raw() == ((long_dim - 1) * long_dim));

	// Test EM tall matrix
	if (safs::is_safs_init()) {
		printf("Test full materialization on EM tall matrix\n");
		bool orig_in_mem = in_mem;
		in_mem = false;
		m1 = create_matrix(long_dim, 1, matrix_layout_t::L_COL,
				num_nodes, get_scalar_type<int>());
		dense_matrix::ptr m2 = create_matrix(long_dim, 1, matrix_layout_t::L_COL,
				num_nodes, get_scalar_type<int>());
		tmp = m1->add(*m2);
		tmp = tmp->cast_ele_type(get_scalar_type<size_t>());
		assert(tmp->get_raw_store()->get_portion_size().first
				== detail::EM_matrix_store::CHUNK_SIZE);
		tmp->set_materialize_level(materialize_level::MATER_FULL);
		res = tmp->sum();
		assert(*(size_t *) res->get_raw() == ((long_dim - 1) * long_dim));

		vstore = std::dynamic_pointer_cast<const detail::mapply_matrix_store>(
				tmp->get_raw_store());
		assert(vstore->get_portion_size().first
				== detail::EM_matrix_store::CHUNK_SIZE);
		assert(vstore);
		assert(!vstore->is_in_mem());
		assert(vstore->has_materialized());
		materialize_res = dense_matrix::create(vstore->materialize(
					vstore->is_in_mem(), vstore->get_num_nodes()));
		res = materialize_res->sum();
		assert(*(size_t *) res->get_raw() == ((long_dim - 1) * long_dim));

		// Test EM wide matrix
		printf("Test full materialization on EM wide matrix\n");
		m1 = create_matrix(1, long_dim, matrix_layout_t::L_COL,
				num_nodes, get_scalar_type<int>());
		m2 = create_matrix(1, long_dim, matrix_layout_t::L_COL,
				num_nodes, get_scalar_type<int>());
		tmp = m1->add(*m2);
		tmp = tmp->cast_ele_type(get_scalar_type<size_t>());
		assert(tmp->get_raw_store()->get_portion_size().second
				== detail::EM_matrix_store::CHUNK_SIZE);
		tmp->set_materialize_level(materialize_level::MATER_FULL);
		res = tmp->sum();
		assert(*(size_t *) res->get_raw() == ((long_dim - 1) * long_dim));

		vstore = std::dynamic_pointer_cast<const detail::mapply_matrix_store>(
				tmp->get_raw_store());
		assert(vstore->get_portion_size().second
				== detail::EM_matrix_store::CHUNK_SIZE);
		assert(vstore);
		assert(!vstore->is_in_mem());
		assert(vstore->has_materialized());
		materialize_res = dense_matrix::create(vstore->materialize(
					vstore->is_in_mem(), vstore->get_num_nodes()));
		res = materialize_res->sum();
		assert(*(size_t *) res->get_raw() == ((long_dim - 1) * long_dim));

		in_mem = orig_in_mem;
	}
}

template<class T>
void print_mat(detail::matrix_store::const_ptr mat)
{
	detail::mem_matrix_store::const_ptr mem_mat
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(mat);
	assert(mem_mat);
	for (size_t i = 0; i < mat->get_num_rows(); i++) {
		for (size_t j = 0; j < mat->get_num_cols(); j++)
			std::cout << mem_mat->get<T>(i, j) << " ";
		std::cout << std::endl;
	}
}

void _test_materialize_all(int num_nodes)
{
	std::vector<dense_matrix::ptr> mats;
	scalar_variable::ptr res;
	dense_matrix::ptr tmp1, tmp2, tmp3;
	dense_matrix::ptr agg1, agg2;
	bulk_operate::const_ptr add = bulk_operate::conv2ptr(
			get_scalar_type<size_t>().get_basic_ops().get_add());
	agg_operate::const_ptr sum = agg_operate::create(add);

	matrix_val_t matrix_val = matrix_val_t::SEQ;
	dense_matrix::ptr m1 = create_matrix(long_dim, 10, matrix_layout_t::L_COL,
			num_nodes, get_scalar_type<int>());
	dense_matrix::ptr m2 = create_matrix(long_dim, 10, matrix_layout_t::L_COL,
			num_nodes, get_scalar_type<int>());
	size_t num_eles = m1->get_num_rows() * m1->get_num_cols();

	printf("compute sum\n");
	tmp1 = m1->add(*m1);
	tmp1 = tmp1->cast_ele_type(get_scalar_type<size_t>());
	res = tmp1->sum();
	printf("compute %ld, expect %ld\n", *(size_t *) res->get_raw(),
			((num_eles - 1) * num_eles));
	assert(*(size_t *) res->get_raw() == ((num_eles - 1) * num_eles));

	printf("materialize one tall matrix\n");
	tmp1 = m1->add(*m1);
	tmp1 = tmp1->cast_ele_type(get_scalar_type<size_t>());
	mats.resize(1);
	mats[0] = tmp1;
	materialize(mats);
	assert(!tmp1->is_virtual());
	res = tmp1->sum();
	printf("compute %ld, expect %ld\n", *(size_t *) res->get_raw(),
			((num_eles - 1) * num_eles));
	assert(*(size_t *) res->get_raw() == ((num_eles - 1) * num_eles));

	printf("materialize two tall matrices\n");
	tmp2 = m2->add(*m2);
	tmp2 = tmp2->cast_ele_type(get_scalar_type<size_t>());
	tmp1 = m1->add(*m1);
	tmp1 = tmp1->cast_ele_type(get_scalar_type<size_t>());

	mats.resize(2);
	mats[0] = tmp1;
	mats[1] = tmp2;
	materialize(mats);
	assert(!tmp1->is_virtual());
	assert(!tmp2->is_virtual());
	res = tmp1->sum();
	assert(*(size_t *) res->get_raw() == ((num_eles - 1) * num_eles));
	res = tmp2->sum();
	assert(*(size_t *) res->get_raw() == ((num_eles - 1) * num_eles));

	printf("materialize one aggregation\n");
	tmp1 = m1->add(*m1);
	tmp1 = tmp1->cast_ele_type(get_scalar_type<size_t>());
	agg1 = tmp1->aggregate(matrix_margin::BOTH, sum);
	assert(agg1->get_num_rows() == 1 && agg1->get_num_cols() == 1);
	assert(agg1->is_virtual());
	mats.resize(1);
	mats[0] = agg1;
	materialize(mats);
	assert(!agg1->is_virtual());
	{
		const detail::mem_matrix_store &mem_agg1
			= dynamic_cast<const detail::mem_matrix_store &>(agg1->get_data());
		num_eles = tmp1->get_num_rows() * tmp1->get_num_cols();
		assert(*(const size_t *) mem_agg1.get_raw_arr() == ((num_eles - 1) * num_eles));
	}

	printf("materialize two aggregations\n");
	tmp1 = m1->add(*m1);
	tmp1 = tmp1->cast_ele_type(get_scalar_type<size_t>());
	tmp2 = m2->add(*m2);
	tmp2 = tmp2->cast_ele_type(get_scalar_type<size_t>());
	agg1 = tmp1->aggregate(matrix_margin::BOTH, sum);
	agg2 = tmp2->aggregate(matrix_margin::BOTH, sum);
	assert(agg1->get_num_rows() == 1 && agg1->get_num_cols() == 1);
	assert(agg2->get_num_rows() == 1 && agg2->get_num_cols() == 1);
	assert(agg1->is_virtual());
	assert(agg2->is_virtual());

	mats.resize(2);
	mats[0] = agg1;
	mats[1] = agg2;
	materialize(mats);
	assert(!agg1->is_virtual());
	assert(!agg2->is_virtual());
	{
		const detail::mem_matrix_store &mem_agg1
			= dynamic_cast<const detail::mem_matrix_store &>(agg1->get_data());
		const detail::mem_matrix_store &mem_agg2
			= dynamic_cast<const detail::mem_matrix_store &>(agg2->get_data());
		num_eles = tmp1->get_num_rows() * tmp1->get_num_cols();
		assert(*(const size_t *) mem_agg1.get_raw_arr() == ((num_eles - 1) * num_eles));
		assert(*(const size_t *) mem_agg2.get_raw_arr() == ((num_eles - 1) * num_eles));
	}

	printf("materialize a tall virtual matrix and two aggregation together\n");
	tmp1 = m1->add(*m1);
	tmp1 = tmp1->cast_ele_type(get_scalar_type<size_t>());
	tmp2 = m2->add(*m2);
	tmp2 = tmp2->cast_ele_type(get_scalar_type<size_t>());
	agg1 = tmp1->aggregate(matrix_margin::BOTH, sum);
	agg2 = tmp2->aggregate(matrix_margin::BOTH, sum);
	assert(agg1->get_num_rows() == 1 && agg1->get_num_cols() == 1);
	assert(agg2->get_num_rows() == 1 && agg2->get_num_cols() == 1);
	assert(agg1->is_virtual());
	assert(agg2->is_virtual());
	tmp3 = m1->add(*m2);
	tmp3 = tmp3->cast_ele_type(get_scalar_type<size_t>());

	mats.resize(3);
	mats[0] = agg1;
	mats[1] = agg2;
	mats[2] = tmp3;
	materialize(mats);
	assert(!agg1->is_virtual());
	assert(!agg2->is_virtual());
	assert(!tmp3->is_virtual());
	{
		const detail::mem_matrix_store &mem_agg1
			= dynamic_cast<const detail::mem_matrix_store &>(agg1->get_data());
		const detail::mem_matrix_store &mem_agg2
			= dynamic_cast<const detail::mem_matrix_store &>(agg2->get_data());
		num_eles = tmp1->get_num_rows() * tmp1->get_num_cols();
		assert(*(const size_t *) mem_agg1.get_raw_arr() == ((num_eles - 1) * num_eles));
		assert(*(const size_t *) mem_agg2.get_raw_arr() == ((num_eles - 1) * num_eles));

		size_t num_eles = tmp3->get_num_rows() * tmp3->get_num_cols();
		res = tmp3->sum();
		printf("sum: %ld, expect: %ld\n", *(size_t *) res->get_raw(),
				((num_eles - 1) * num_eles));
		assert(*(size_t *) res->get_raw() == ((num_eles - 1) * num_eles));
	}
}

void test_materialize_all(int num_nodes)
{
	// Test the in-mem case.
	printf("test materialization on multiple virtual matrices in memory\n");
	_test_materialize_all(num_nodes);
	printf("test materialization on multiple virtual block matrices in memory\n");
	block_size = 3;
	_test_materialize_all(num_nodes);
	// Test the EM case.
	if (safs::is_safs_init()) {
		in_mem = false;
		block_size = 0;
		printf("test materialization on multiple virtual matrices in EM\n");
		_test_materialize_all(num_nodes);
		printf("test materialization on multiple virtual block matrices in EM\n");
		block_size = 3;
		_test_materialize_all(num_nodes);
		block_size = 0;
		in_mem = true;
	}
}

void test_bmv_multiply_tall()
{
	if (!safs::is_safs_init())
		return;

	bool in_mem = false;
	printf("gemm tall on block multi-vector\n");
	eigen::block_multi_vector::ptr mv = eigen::block_multi_vector::create(
			long_dim, 14, 2, get_scalar_type<double>(), in_mem, false);
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
				get_scalar_type<double>(), true, false);
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
				get_scalar_type<double>(), true, false);
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
				get_scalar_type<double>(), true, false);
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
				get_scalar_type<double>(), true, false);
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
	if (!safs::is_safs_init())
		return;

	bool in_mem = false;
	printf("gemm wide on block multi-vector\n");
	eigen::block_multi_vector::ptr mv1 = eigen::block_multi_vector::create(
			long_dim, 14, 2, get_scalar_type<double>(), in_mem, false);
	for (size_t i = 0; i < mv1->get_num_blocks(); i++) {
		dense_matrix::ptr mat = dense_matrix::create_randu<int>(0, 10, long_dim,
				mv1->get_block_size(), matrix_layout_t::L_COL, -1, in_mem);
		mv1->set_block(i, mat->cast_ele_type(get_scalar_type<double>()));
	}

	eigen::block_multi_vector::ptr mv2 = eigen::block_multi_vector::create(
			long_dim, 2, 2, get_scalar_type<double>(), in_mem, false);
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

void test_bind(int num_nodes)
{
	printf("test matrix rbind\n");
	std::vector<dense_matrix::ptr> mats(3);
	std::vector<scalar_variable::ptr> sums(mats.size());
	int tot_sum = 0;
	std::vector<const void *> rows;
	for (size_t i = 0; i < mats.size(); i++) {
		mats[i] = dense_matrix::create_seq<int>(0, 1, random() % 1000000, 10,
				matrix_layout_t::L_ROW, true, num_nodes);
		sums[i] = mats[i]->sum();
		tot_sum += scalar_variable::get_val<int>(*sums[i]);
		mats[i]->materialize_self();
		detail::mem_matrix_store::const_ptr store
			= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
					mats[i]->get_raw_store());
		assert(store);
		for (size_t j = 0; j < store->get_num_rows(); j++)
			rows.push_back(store->get_row(j));
	}
	dense_matrix::ptr combined = dense_matrix::rbind(mats);
	assert(combined->store_layout() == matrix_layout_t::L_ROW);
	detail::mem_matrix_store::const_ptr combined_store
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				combined->get_raw_store());
	assert(combined_store);
	scalar_variable::ptr sum = combined->sum();
	assert(tot_sum == scalar_variable::get_val<int>(*sum));
	assert(rows.size() == combined->get_num_rows());
	for (size_t i = 0; i < rows.size(); i++)
		assert(memcmp(rows[i], combined_store->get_row(i),
					combined->get_num_cols() * combined->get_entry_size()) == 0);

	printf("test matrix cbined\n");
	for (size_t i = 0; i < mats.size(); i++)
		mats[i] = mats[i]->transpose();
	combined = dense_matrix::cbind(mats);
	assert(combined->store_layout() == matrix_layout_t::L_COL);
	combined_store = std::dynamic_pointer_cast<const detail::mem_matrix_store>(
			combined->get_raw_store());
	assert(combined_store);
	sum = combined->sum();
	assert(tot_sum == scalar_variable::get_val<int>(*sum));
	assert(rows.size() == combined->get_num_cols());
	for (size_t i = 0; i < rows.size(); i++)
		assert(memcmp(rows[i], combined_store->get_col(i),
					combined->get_num_rows() * combined->get_entry_size()) == 0);

	printf("test block matrix rbind\n");
	set_operate::const_ptr rand_init = create_urand_init<int>(0, 1000);
	tot_sum = 0;
	for (size_t i = 0; i < mats.size(); i++) {
		mats[i] = block_matrix::create(random() % 1000000 + 100, 9, 4,
				get_scalar_type<int>(), *rand_init, num_nodes);
		sums[i] = mats[i]->sum();
		tot_sum += scalar_variable::get_val<int>(*sums[i]);
	}
	combined = dense_matrix::rbind(mats);
	sum = combined->sum();
	assert(tot_sum == scalar_variable::get_val<int>(*sum));

	printf("test block matrix cbind\n");
	tot_sum = 0;
	for (size_t i = 0; i < mats.size(); i++) {
		mats[i] = mats[i]->transpose();
		sums[i] = mats[i]->sum();
		tot_sum += scalar_variable::get_val<int>(*sums[i]);
	}
	combined = dense_matrix::cbind(mats);
	sum = combined->sum();
	assert(tot_sum == scalar_variable::get_val<int>(*sum));
}

void _test_seq_byrow(dense_matrix::ptr mat)
{
	// Test rows
	size_t num_tests = std::min(mat->get_num_rows(), 100UL);
	for (size_t i = 0; i < num_tests; i++) {
		off_t row_idx = random() % mat->get_num_rows();
		std::vector<off_t> rows(1, row_idx);
		dense_matrix::ptr row = mat->get_rows(rows);
		auto sum = row->sum();
		size_t start = row_idx * mat->get_num_cols();
		size_t end = start + (mat->get_num_cols() - 1);
		assert(scalar_variable::get_val<size_t>(*sum)
				== (start + end) * mat->get_num_cols() / 2);
	}
	// Test cols
	num_tests = std::min(mat->get_num_cols(), 100UL);
	for (size_t i = 0; i < num_tests; i++) {
		off_t col_idx = random() % mat->get_num_cols();
		std::vector<off_t> cols(1, col_idx);
		dense_matrix::ptr col = mat->get_cols(cols);
		auto sum = col->sum();
		size_t start = col_idx;
		size_t end = col_idx + (mat->get_num_rows() - 1) * mat->get_num_cols();
		assert(scalar_variable::get_val<size_t>(*sum)
				== (start + end) * mat->get_num_rows() / 2);
	}
}

void _test_seq_bycol(dense_matrix::ptr mat)
{
	// Test rows
	size_t num_tests = std::min(mat->get_num_rows(), 100UL);
	for (size_t i = 0; i < num_tests; i++) {
		off_t row_idx = random() % mat->get_num_rows();
		std::vector<off_t> rows(1, row_idx);
		dense_matrix::ptr row = mat->get_rows(rows);
		auto sum = row->sum();
		size_t start = row_idx;
		size_t end = start + (mat->get_num_cols() - 1) * mat->get_num_rows();
		assert(scalar_variable::get_val<size_t>(*sum)
				== (start + end) * mat->get_num_cols() / 2);
	}
	// Test cols
	num_tests = std::min(mat->get_num_cols(), 100UL);
	for (size_t i = 0; i < num_tests; i++) {
		off_t col_idx = random() % mat->get_num_cols();
		std::vector<off_t> cols(1, col_idx);
		dense_matrix::ptr col = mat->get_cols(cols);
		auto sum = col->sum();
		size_t start = col_idx * mat->get_num_rows();
		size_t end = start + (mat->get_num_rows() - 1);
		assert(scalar_variable::get_val<size_t>(*sum)
				== (start + end) * mat->get_num_rows() / 2);
	}
}

void test_seq_matrix()
{
	dense_matrix::ptr mat, tmp;

	printf("test matrices with sequence numbers\n");
	mat = dense_matrix::create_seq<size_t>(0, 1, long_dim, 10,
			matrix_layout_t::L_COL, true);
	_test_seq_byrow(mat);
	mat = mat->transpose();
	_test_seq_bycol(mat);

	mat = dense_matrix::create_seq<size_t>(0, 1, long_dim, 10,
			matrix_layout_t::L_ROW, true);
	_test_seq_byrow(mat);
	mat = mat->transpose();
	_test_seq_bycol(mat);

	mat = dense_matrix::create_seq<size_t>(0, 1, long_dim, 10,
			matrix_layout_t::L_COL, false);
	_test_seq_bycol(mat);
	mat = mat->transpose();
	_test_seq_byrow(mat);

	mat = dense_matrix::create_seq<size_t>(0, 1, long_dim, 10,
			matrix_layout_t::L_ROW, false);
	_test_seq_bycol(mat);
	mat = mat->transpose();
	_test_seq_byrow(mat);

	printf("test block matrices with sequence numbers\n");
	printf("create tall block matrix1\n");
	mat = dense_matrix::create_seq<size_t>(0, 1, long_dim / 10, 100,
			matrix_layout_t::L_COL, true);
	assert(mat->get_data().is_virtual());
	_test_seq_byrow(mat);
	mat = mat->transpose();
	_test_seq_bycol(mat);

	printf("create tall block matrix2\n");
	mat = dense_matrix::create_seq<size_t>(0, 1, long_dim / 10, 100,
			matrix_layout_t::L_COL, false);
	assert(mat->get_data().is_virtual());
	_test_seq_bycol(mat);
	mat = mat->transpose();
	_test_seq_byrow(mat);

	printf("create wide block matrix1\n");
	mat = dense_matrix::create_seq<size_t>(0, 1, 100, long_dim / 10,
			matrix_layout_t::L_COL, true);
	assert(mat->get_data().is_virtual());
	_test_seq_byrow(mat);
	mat = mat->transpose();
	_test_seq_bycol(mat);

	printf("create wide block matrix2\n");
	mat = dense_matrix::create_seq<size_t>(0, 1, 100, long_dim / 10,
			matrix_layout_t::L_COL, false);
	assert(mat->get_data().is_virtual());
	_test_seq_bycol(mat);
	mat = mat->transpose();
	_test_seq_byrow(mat);
}

void _test_repeat_byrow(dense_matrix::ptr mat)
{
	// Test rows
	size_t num_tests = std::min(mat->get_num_rows(), 100UL);
	for (size_t i = 0; i < num_tests; i++) {
		off_t row_idx = random() % mat->get_num_rows();
		std::vector<off_t> rows(1, row_idx);
		dense_matrix::ptr row = mat->get_rows(rows);
		auto sum = row->sum();
		assert(scalar_variable::get_val<size_t>(*sum)
				== (mat->get_num_cols() - 1) * mat->get_num_cols() / 2);
	}
	// Test cols
	num_tests = std::min(mat->get_num_cols(), 100UL);
	for (size_t i = 0; i < num_tests; i++) {
		off_t col_idx = random() % mat->get_num_cols();
		std::vector<off_t> cols(1, col_idx);
		dense_matrix::ptr col = mat->get_cols(cols);
		auto sum = col->sum();
		assert(scalar_variable::get_val<size_t>(*sum)
				== col_idx * mat->get_num_rows());
	}
}

void _test_repeat_bycol(dense_matrix::ptr mat)
{
	// Test rows
	size_t num_tests = std::min(mat->get_num_rows(), 100UL);
	for (size_t i = 0; i < num_tests; i++) {
		off_t row_idx = random() % mat->get_num_rows();
		std::vector<off_t> rows(1, row_idx);
		dense_matrix::ptr row = mat->get_rows(rows);
		auto sum = row->sum();
		assert(scalar_variable::get_val<size_t>(*sum)
				== row_idx * mat->get_num_cols());
	}
	// Test cols
	num_tests = std::min(mat->get_num_cols(), 100UL);
	for (size_t i = 0; i < num_tests; i++) {
		off_t col_idx = random() % mat->get_num_cols();
		std::vector<off_t> cols(1, col_idx);
		dense_matrix::ptr col = mat->get_cols(cols);
		auto sum = col->sum();
		assert(scalar_variable::get_val<size_t>(*sum)
				== (mat->get_num_rows() - 1) * mat->get_num_rows() / 2);
	}
}

void test_repeat()
{
	col_vec::ptr vec;
	dense_matrix::ptr mat;

	printf("test matrices with repeated vectors\n");
	// skinny matrix.
	vec = col_vec::create(dense_matrix::create_seq<size_t>(0, 1, 10, 1,
				matrix_layout_t::L_COL, false));
	mat = dense_matrix::create_repeat(vec, long_dim, vec->get_length(),
			matrix_layout_t::L_ROW, true);
	_test_repeat_byrow(mat);
	mat = mat->transpose();
	_test_repeat_bycol(mat);

	mat = dense_matrix::create_repeat(vec, long_dim, vec->get_length(),
			matrix_layout_t::L_COL, true);
	_test_repeat_byrow(mat);
	mat = mat->transpose();
	_test_repeat_bycol(mat);

	vec = col_vec::create(dense_matrix::create_seq<size_t>(0, 1, long_dim, 1,
				matrix_layout_t::L_COL, false));
	mat = dense_matrix::create_repeat(vec, vec->get_length(), 10,
			matrix_layout_t::L_ROW, false);
	_test_repeat_bycol(mat);
	mat = mat->transpose();
	_test_repeat_byrow(mat);

	mat = dense_matrix::create_repeat(vec, vec->get_length(), 10,
			matrix_layout_t::L_COL, false);
	_test_repeat_bycol(mat);
	mat = mat->transpose();
	_test_repeat_byrow(mat);

	printf("test block matrices with repeated vectors\n");
	// block matrix
	vec = col_vec::create(dense_matrix::create_seq<size_t>(0, 1, 100, 1,
				matrix_layout_t::L_COL, false));
	mat = dense_matrix::create_repeat(vec, long_dim, vec->get_length(),
			matrix_layout_t::L_ROW, true);
	_test_repeat_byrow(mat);
	mat = mat->transpose();
	_test_repeat_bycol(mat);

	mat = dense_matrix::create_repeat(vec, long_dim, vec->get_length(),
			matrix_layout_t::L_COL, true);
	_test_repeat_byrow(mat);
	mat = mat->transpose();
	_test_repeat_bycol(mat);

	vec = col_vec::create(dense_matrix::create_seq<size_t>(0, 1, long_dim, 1,
				matrix_layout_t::L_COL, false));
	mat = dense_matrix::create_repeat(vec, vec->get_length(), 100,
			matrix_layout_t::L_ROW, false);
	_test_repeat_bycol(mat);
	mat = mat->transpose();
	_test_repeat_byrow(mat);

	mat = dense_matrix::create_repeat(vec, vec->get_length(), 100,
			matrix_layout_t::L_COL, false);
	_test_repeat_bycol(mat);
	mat = mat->transpose();
	_test_repeat_byrow(mat);
}

void _test_share_data(dense_matrix::ptr mat1, dense_matrix::ptr mat2)
{
	std::vector<off_t> idxs(3);
	for (size_t i = 0; i < idxs.size(); i++)
		idxs[i] = random() % 10;

	assert(mat1->get_raw_store()->share_data(*mat1->get_raw_store()));
	assert(!mat1->get_raw_store()->share_data(*mat2->get_raw_store()));
	mat2 = mat1->transpose();
	assert(mat1->get_raw_store()->share_data(*mat2->get_raw_store()));
	mat2 = mat1->get_cols(idxs);
	assert(!mat1->get_raw_store()->share_data(*mat2->get_raw_store()));
	assert(!mat2->get_raw_store()->share_data(*mat1->get_raw_store()));
	dense_matrix::ptr mat3 = mat2->transpose();
	assert(mat3->get_raw_store()->share_data(*mat2->get_raw_store()));
	assert(mat2->get_raw_store()->share_data(*mat3->get_raw_store()));
}

void test_share_data()
{
	dense_matrix::ptr mat1, mat2;

	// Test EM matrix.
	mat1 = dense_matrix::create_randu<size_t>(0, 1000, long_dim, 10,
			matrix_layout_t::L_ROW, -1, false);
	mat2 = dense_matrix::create_randu<size_t>(0, 1000, long_dim, 10,
			matrix_layout_t::L_ROW, -1, false);
	printf("test share data in EM matrix\n");
	_test_share_data(mat1, mat2);

	// Test SMP matrix.
	mat1 = dense_matrix::create_randu<size_t>(0, 1000, long_dim, 10,
			matrix_layout_t::L_ROW, -1, true);
	mat2 = dense_matrix::create_randu<size_t>(0, 1000, long_dim, 10,
			matrix_layout_t::L_ROW, -1, true);
	printf("test share data in SMP matrix\n");
	_test_share_data(mat1, mat2);

	// Test NUMA matrix.
	mat1 = dense_matrix::create_randu<size_t>(0, 1000, long_dim, 10,
			matrix_layout_t::L_ROW, matrix_conf.get_num_nodes(), true);
	mat2 = dense_matrix::create_randu<size_t>(0, 1000, long_dim, 10,
			matrix_layout_t::L_ROW, matrix_conf.get_num_nodes(), true);
	printf("test share data in NUMA matrix\n");
	_test_share_data(mat1, mat2);

	// Test mapply matrix.
	mat1 = mat1->conv2(matrix_layout_t::L_COL);
	mat2 = mat2->conv2(matrix_layout_t::L_COL);
	assert(mat1->get_raw_store()->is_virtual());
	assert(mat2->get_raw_store()->is_virtual());
	printf("test share data in mapply matrix\n");
	_test_share_data(mat1, mat2);

	// Test set_data matrix.
	mat1 = dense_matrix::create_seq<size_t>(0, 1, long_dim, 10,
			matrix_layout_t::L_ROW, true, -1, true);
	mat2 = dense_matrix::create_seq<size_t>(0, 1, long_dim, 10,
			matrix_layout_t::L_ROW, true, -1, true);
	printf("test share data in set_data matrix\n");
	_test_share_data(mat1, mat2);

	// Test one-val matrix.
	mat1 = dense_matrix::create_const<size_t>(0, long_dim, 10,
			matrix_layout_t::L_ROW, -1, true);
	mat2 = dense_matrix::create_const<size_t>(1, long_dim, 10,
			matrix_layout_t::L_ROW, -1, true);
	printf("test share data in one-val matrix\n");
	_test_share_data(mat1, mat2);
	mat2 = dense_matrix::create_const<size_t>(0, long_dim, 10,
			matrix_layout_t::L_ROW, -1, true);
	assert(mat1->get_raw_store()->share_data(*mat2->get_raw_store()));

	// Test block matrix
	mat1 = dense_matrix::create_randu<size_t>(0, 1000, long_dim, 100,
			matrix_layout_t::L_ROW, matrix_conf.get_num_nodes(), true);
	mat2 = dense_matrix::create_randu<size_t>(0, 1000, long_dim, 100,
			matrix_layout_t::L_ROW, matrix_conf.get_num_nodes(), true);
	printf("test share data in block matrix\n");
	_test_share_data(mat1, mat2);
}

void test_factor()
{
	dense_matrix::ptr mat = dense_matrix::create_randu<int>(1, 1000,
			long_dim, 1, matrix_layout_t::L_COL);
	factor_col_vector::ptr fvec = factor_col_vector::create(mat);
	assert(fvec);
	assert(fvec->get_num_levels() <= 1000);
}

void test_cross_prod()
{
	dense_matrix::ptr mat1 = dense_matrix::create_randu<double>(1, 1000,
			long_dim, 10, matrix_layout_t::L_COL);
	dense_matrix::ptr mat2 = dense_matrix::create_randu<double>(1, 1000,
			9, long_dim, matrix_layout_t::L_COL);
	dense_matrix::ptr res = mat2->multiply(*mat1);
	dense_matrix::ptr tres = res->transpose();
	assert(res->get_raw_store()->get_data_id()
			== tres->get_raw_store()->get_data_id());
	printf("materialize tres\n");
	tres->materialize_self();
	printf("materialize res\n");
	res->materialize_self();
	printf("complete crossprod\n");
}

void test_ref_cnts(int num_nodes)
{
	std::vector<dense_matrix::ptr> mats;

	dense_matrix::ptr mat1 = dense_matrix::create_randu<double>(1, 1000,
			long_dim, 10, matrix_layout_t::L_COL, num_nodes);
	dense_matrix::ptr mat2 = dense_matrix::create_randu<double>(1, 1000,
			long_dim, 10, matrix_layout_t::L_COL, num_nodes);
	dense_matrix::ptr mat3 = mat1->add(*mat2);
	printf("mat3: %s\n", mat3->get_data().get_name().c_str());
	mats.push_back(mat1);
	mats.push_back(mat2);
	mats.push_back(mat3);

	std::vector<dense_matrix::ptr> tmp;
	tmp.resize(3);
	tmp[0] = mat1;
	tmp[1] = mat2;
	tmp[2] = mat3;
	dense_matrix::ptr mat4 = dense_matrix::cbind(tmp);
	printf("mat4: %s\n", mat4->get_data().get_name().c_str());
	mats.push_back(mat4);

	dense_matrix::ptr mat5 = mat4->row_sum();
	printf("mat5: %s\n", mat5->get_data().get_name().c_str());
	mats.push_back(mat5);
	const_cast<detail::matrix_store &>(mat5->get_data()).inc_dag_ref(
			detail::INVALID_MAT_ID);

	dense_matrix::ptr mat6 = mat4->col_sum();
	const_cast<detail::matrix_store &>(mat6->get_data()).inc_dag_ref(
			detail::INVALID_MAT_ID);
	dense_matrix::ptr mat7 = mat1->col_sum();
	const_cast<detail::matrix_store &>(mat7->get_data()).inc_dag_ref(
			detail::INVALID_MAT_ID);
	for (size_t i = 0; i < mats.size(); i++)
		printf("mat %ld: %ld (%s)\n", i, mats[i]->get_data().get_dag_ref(),
				mats[i]->get_data().get_name().c_str());
	assert(mat1->get_data().get_dag_ref() == 3);
	assert(mat2->get_data().get_dag_ref() == 2);
	assert(mat3->get_data().get_dag_ref() == 1);
	assert(mat4->get_data().get_dag_ref() == 2);
	assert(mat5->get_data().get_dag_ref() == 0);

	const_cast<detail::matrix_store &>(mat6->get_data()).reset_dag_ref();
	for (size_t i = 0; i < mats.size(); i++)
		printf("mat %ld: %ld\n", i, mats[i]->get_data().get_dag_ref());
	assert(mat1->get_data().get_dag_ref() == 0);
	assert(mat2->get_data().get_dag_ref() == 0);
	assert(mat3->get_data().get_dag_ref() == 0);
	assert(mat4->get_data().get_dag_ref() == 0);
	assert(mat5->get_data().get_dag_ref() == 0);
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

	test_ref_cnts(num_nodes);
	test_set_rowcols();
	test_cross_prod();
	test_factor();
	test_share_data();
	test_repeat_rowcols();
	test_repeat();
	test_seq_matrix();
	test_bind(num_nodes);
	test_materialize_all(num_nodes);
	test_materialize(num_nodes);
	test_get_rowcols(num_nodes);
	test_block_mv();
	test_conv_store();
	test_mapply_mixed(num_nodes);
	test_mem_matrix(num_nodes);
	test_EM_matrix();
	long_dim = 9999;
	test_mem_matrix(num_nodes);

	destroy_flash_matrix();
}
