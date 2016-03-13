/*
 * Copyright 2016 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "sparse_matrix.h"
#include "block_matrix.h"
#include "combined_matrix_store.h"

using namespace fm;

template<class T>
void check_mat_equal(dense_matrix::ptr m1, dense_matrix::ptr m2)
{
	assert(m1->get_num_rows() == m2->get_num_rows());
	assert(m1->get_num_cols() == m2->get_num_cols());
	block_matrix::ptr block_m = std::dynamic_pointer_cast<block_matrix>(m1);
	if (block_m) {
		m1 = dense_matrix::create(block_m->get_raw_store());
		m1->materialize_self();
	}
	scalar_variable::const_ptr res = m1->minus(*m2)->abs()->max();
	assert(scalar_variable::get_val<T>(*res) == 0);
}

template<class T>
void check_mat_approx(dense_matrix::ptr m1, dense_matrix::ptr m2)
{
	assert(m1->get_num_rows() == m2->get_num_rows());
	assert(m1->get_num_cols() == m2->get_num_cols());
	block_matrix::ptr block_m = std::dynamic_pointer_cast<block_matrix>(m1);
	if (block_m) {
		m1 = dense_matrix::create(block_m->get_raw_store());
		m1->materialize_self();
	}
	dense_matrix::ptr diff = m1->minus(*m2)->abs();
	dense_matrix::ptr max = m1->abs()->pmax(*m2->abs());
	scalar_variable::const_ptr res = diff->div(*max)->max();
	std::cout << "max diff: " << scalar_variable::get_val<T>(*res) << std::endl;
	assert(scalar_variable::get_val<T>(*res) < 10e-14);
}

template<class T>
void check_vec_equal(vector::ptr v1, vector::ptr v2)
{
	assert(v1->get_length() == v2->get_length());
	assert(v1->equals(*v2));
}

bool is_block_matrix(dense_matrix::ptr mat)
{
	block_matrix::ptr block_m = std::dynamic_pointer_cast<block_matrix>(mat);
	return block_m != NULL;
}

void test_mapply2()
{
	printf("test mapply2\n");
	dense_matrix::ptr mat1 = block_matrix::create(10000, 10, 3,
			get_scalar_type<int>(), const_set_operate<int>(1));
	dense_matrix::ptr mat2 = block_matrix::create(10000, 10, 3,
			get_scalar_type<int>(), const_set_operate<int>(1));
	dense_matrix::ptr res = mat1->add(*mat2);
	assert(is_block_matrix(res));
	res->materialize_self();

	dense_matrix::ptr tmp = dense_matrix::create(10000, 10,
			matrix_layout_t::L_COL, get_scalar_type<int>(),
			const_set_operate<int>(2));
	check_mat_equal<int>(res, tmp);
}

void test_get_rc()
{
	printf("test getting rows and columns\n");
	set_operate::const_ptr init_op = create_urand_init<double>(0, 1);

	// check getting columns.
	dense_matrix::ptr mat1 = block_matrix::create(10000, 10, 3,
			init_op->get_type(), *init_op);
	dense_matrix::ptr mat2 = dense_matrix::create(mat1->get_raw_store());
	mat2->materialize_self();

	assert(is_block_matrix(mat1));
	assert(!is_block_matrix(mat2));

	vector::ptr vec1 = mat1->get_col(4);
	vector::ptr vec2 = mat2->get_col(4);
	check_vec_equal<double>(vec1, vec2);

	std::vector<off_t> offs(3);
	offs[0] = 1;
	offs[1] = 3;
	offs[2] = 9;
	dense_matrix::ptr sub_m1 = mat1->get_cols(offs);
	dense_matrix::ptr sub_m2 = mat2->get_cols(offs);
	check_mat_equal<double>(sub_m1, sub_m2);

	// check getting rows.
	mat1 = block_matrix::create(10, 10000, 3, init_op->get_type(), *init_op);
	mat2 = dense_matrix::create(mat1->get_raw_store());
	mat2->materialize_self();

	assert(is_block_matrix(mat1));
	assert(!is_block_matrix(mat2));

	vec1 = mat1->get_row(4);
	vec2 = mat2->get_row(4);
	check_vec_equal<double>(vec1, vec2);

	sub_m1 = mat1->get_rows(offs);
	sub_m2 = mat2->get_rows(offs);
	check_mat_equal<double>(sub_m1, sub_m2);
}

void test_transpose()
{
	printf("test transpose\n");
	set_operate::const_ptr init_op = create_urand_init<double>(0, 1);

	dense_matrix::ptr mat1 = block_matrix::create(10000, 10, 3,
			init_op->get_type(), *init_op);
	dense_matrix::ptr mat2 = dense_matrix::create(mat1->get_raw_store());
	mat2->materialize_self();

	mat1 = mat1->transpose();
	mat2 = mat2->transpose();

	assert(is_block_matrix(mat1));
	assert(!is_block_matrix(mat2));
	assert(mat1->get_num_rows() == mat2->get_num_rows());
	assert(mat1->get_num_cols() == mat2->get_num_cols());
	check_mat_equal<double>(mat1, mat2);
}

template<class T>
class vec_init: public type_set_vec_operate<T>
{
public:
	virtual void set(T *arr, size_t num_eles, off_t start_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = start_idx + i;
	}
};

void test_mapply_rows()
{
	printf("test mapply rows\n");
	set_operate::const_ptr init_op = create_urand_init<double>(0, 1);
	const bulk_operate &op = init_op->get_type().get_basic_ops().get_add();

	// Tall matrix
	dense_matrix::ptr mat1 = block_matrix::create(10000, 10, 3,
			init_op->get_type(), *init_op);
	dense_matrix::ptr mat2 = dense_matrix::create(mat1->get_raw_store());
	mat2->materialize_self();

	vector::ptr vec = vector::create(mat1->get_num_cols(),
			get_scalar_type<double>(), -1, true, vec_init<double>());
	mat1 = mat1->mapply_rows(vec, bulk_operate::conv2ptr(op));
	mat2 = mat2->mapply_rows(vec, bulk_operate::conv2ptr(op));
	assert(is_block_matrix(mat1));
	assert(!is_block_matrix(mat2));
	check_mat_equal<double>(mat1, mat2);

	// Wide matrix
	mat1 = block_matrix::create(10, 10000, 3, init_op->get_type(), *init_op);
	mat2 = dense_matrix::create(mat1->get_raw_store());
	mat2->materialize_self();

	vec = vector::create(mat1->get_num_cols(), get_scalar_type<double>(),
			-1, true, vec_init<double>());
	mat1 = mat1->mapply_rows(vec, bulk_operate::conv2ptr(op));
	mat2 = mat2->mapply_rows(vec, bulk_operate::conv2ptr(op));
	assert(is_block_matrix(mat1));
	assert(!is_block_matrix(mat2));
	check_mat_equal<double>(mat1, mat2);
}

void test_mapply_cols()
{
	printf("test mapply columns\n");
	set_operate::const_ptr init_op = create_urand_init<double>(0, 1);
	const bulk_operate &op = init_op->get_type().get_basic_ops().get_add();

	// Tall matrix
	dense_matrix::ptr mat1 = block_matrix::create(10000, 10, 3,
			init_op->get_type(), *init_op);
	dense_matrix::ptr mat2 = dense_matrix::create(mat1->get_raw_store());
	mat2->materialize_self();

	vector::ptr vec = vector::create(mat1->get_num_rows(),
			get_scalar_type<double>(), -1, true, vec_init<double>());
	mat1 = mat1->mapply_cols(vec, bulk_operate::conv2ptr(op));
	mat2 = mat2->mapply_cols(vec, bulk_operate::conv2ptr(op));
	assert(is_block_matrix(mat1));
	assert(!is_block_matrix(mat2));
	check_mat_equal<double>(mat1, mat2);

	// Wide matrix
	mat1 = block_matrix::create(10, 10000, 3, init_op->get_type(), *init_op);
	mat2 = dense_matrix::create(mat1->get_raw_store());
	mat2->materialize_self();

	vec = vector::create(mat1->get_num_rows(), get_scalar_type<double>(),
			-1, true, vec_init<double>());
	mat1 = mat1->mapply_cols(vec, bulk_operate::conv2ptr(op));
	mat2 = mat2->mapply_cols(vec, bulk_operate::conv2ptr(op));
	assert(is_block_matrix(mat1));
	assert(!is_block_matrix(mat2));
	check_mat_equal<double>(mat1, mat2);
}

void test_sapply()
{
	printf("test sapply\n");
	set_operate::const_ptr init_op = create_urand_init<double>(0, 1);
	const bulk_uoperate &op
		= *init_op->get_type().get_basic_uops().get_op(basic_uops::op_idx::NEG);

	dense_matrix::ptr mat1 = block_matrix::create(10000, 10, 3,
			init_op->get_type(), *init_op);
	dense_matrix::ptr mat2 = dense_matrix::create(mat1->get_raw_store());
	mat2->materialize_self();

	mat1 = mat1->sapply(bulk_uoperate::conv2ptr(op));
	mat2 = mat2->sapply(bulk_uoperate::conv2ptr(op));
	assert(is_block_matrix(mat1));
	assert(!is_block_matrix(mat2));
	check_mat_equal<double>(mat1, mat2);
}

template<class T>
void print_mat(detail::matrix_store::const_ptr mat)
{
	detail::mem_matrix_store::const_ptr mem_mat
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(mat);
	assert(mem_mat);
	for (size_t i = 0; i < mem_mat->get_num_rows(); i++) {
		for (size_t j = 0; j < mem_mat->get_num_cols(); j++)
			std::cout << mem_mat->get<T>(i, j) << ", ";
		std::cout << std::endl;
	}
}

void test_inner_prod()
{
	set_operate::const_ptr init_op = create_urand_init<size_t>(0, 10000);
	const bulk_operate &mul_op = init_op->get_type().get_basic_ops().get_multiply();
	const bulk_operate &add_op = init_op->get_type().get_basic_ops().get_add();

	printf("test inner product on tall matrix\n");
	dense_matrix::ptr mat1 = block_matrix::create(10000, 10, 3,
			init_op->get_type(), *init_op);
	dense_matrix::ptr mat2 = dense_matrix::create_randu<size_t>(0, 10000,
			10, 11, matrix_layout_t::L_ROW);
	dense_matrix::ptr res1 = mat1->inner_prod(*mat2,
			bulk_operate::conv2ptr(mul_op), bulk_operate::conv2ptr(add_op));

	dense_matrix::ptr mat3 = dense_matrix::create(mat1->get_raw_store());
	mat3->materialize_self();
	dense_matrix::ptr res2 = mat3->inner_prod(*mat2,
			bulk_operate::conv2ptr(mul_op), bulk_operate::conv2ptr(add_op));

	assert(is_block_matrix(res1));
	assert(!is_block_matrix(res2));
	check_mat_equal<size_t>(res1, res2);

	printf("test inner product on wide matrix\n");
	mat1 = block_matrix::create(10, 10000, 3, init_op->get_type(), *init_op);
	mat2 = block_matrix::create(10000, 5, 3, init_op->get_type(), *init_op);
	res1 = mat1->inner_prod(*mat2,
			bulk_operate::conv2ptr(mul_op), bulk_operate::conv2ptr(add_op));

	mat3 = dense_matrix::create(mat1->get_raw_store());
	mat3->materialize_self();
	dense_matrix::ptr mat4 = dense_matrix::create(mat2->get_raw_store());
	mat4->materialize_self();
	res2 = mat3->inner_prod(*mat4,
			bulk_operate::conv2ptr(mul_op), bulk_operate::conv2ptr(add_op));
	res2->materialize_self();
	check_mat_equal<size_t>(res1, res2);
}

void test_multiply()
{
	set_operate::const_ptr init_op = create_urand_init<double>(0, 1);

	printf("test multiply on tall matrix\n");
	dense_matrix::ptr mat1 = block_matrix::create(10000, 10, 3,
			init_op->get_type(), *init_op);
	dense_matrix::ptr mat2 = dense_matrix::create_randu<double>(0, 1,
			10, 11, matrix_layout_t::L_ROW);
	dense_matrix::ptr res1 = mat1->multiply(*mat2);

	dense_matrix::ptr mat3 = dense_matrix::create(mat1->get_raw_store());
	mat3->materialize_self();
	dense_matrix::ptr res2 = mat3->multiply(*mat2);

	assert(is_block_matrix(res1));
	assert(!is_block_matrix(res2));
	check_mat_approx<double>(res1, res2);

	printf("test multiply on wide matrix\n");
	mat1 = block_matrix::create(10, 10000, 3, init_op->get_type(), *init_op);
	mat2 = block_matrix::create(10000, 5, 3, init_op->get_type(), *init_op);
	res1 = mat1->multiply(*mat2);

	mat3 = dense_matrix::create(mat1->get_raw_store());
	mat3->materialize_self();
	dense_matrix::ptr mat4 = dense_matrix::create(mat2->get_raw_store());
	mat4->materialize_self();
	res2 = mat3->multiply(*mat4);
	res2->materialize_self();
	check_mat_approx<double>(res1, res2);
}

int main()
{
	init_flash_matrix(NULL);
	test_mapply2();
	test_get_rc();
	test_transpose();
	test_mapply_rows();
	test_mapply_cols();
	test_sapply();
	test_inner_prod();
	test_multiply();
	destroy_flash_matrix();
}
