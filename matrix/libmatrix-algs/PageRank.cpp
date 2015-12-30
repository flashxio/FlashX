/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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

#include "matrix_algs.h"

namespace fm
{

namespace alg
{

namespace
{

class mat_init_one: public type_set_operate<float>
{
public:
	virtual void set(float *arr, size_t num_eles, off_t row_idx,
			            off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = 1;
	}
};

class mat_init_operate: public type_set_operate<float>
{
	size_t num_vertices;
public:
	mat_init_operate(size_t num_vertices) {
		this->num_vertices = num_vertices;
	}

	virtual void set(float *arr, size_t num_eles, off_t row_idx,
			            off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = 1.0 / num_vertices;
	}
};

detail::matrix_store::ptr get_out_degree(sparse_matrix::ptr mat, bool in_mem)
{
	int num_nodes = matrix_conf.get_num_nodes();
	detail::matrix_store::ptr out_deg = detail::matrix_store::create(
			mat->get_num_rows(), 1, matrix_layout_t::L_ROW,
			get_scalar_type<float>(), num_nodes, in_mem);
	detail::matrix_store::ptr one = detail::matrix_store::create(
			mat->get_num_cols(), 1, matrix_layout_t::L_ROW,
			get_scalar_type<float>(), num_nodes, true);
	one->set_data(mat_init_one());
	mat->multiply<float, float>(one, out_deg);
	return out_deg;
}

}

dense_matrix::ptr PageRank(sparse_matrix::ptr mat, size_t max_niters,
		float damping_factor, size_t num_in_mem)
{
	int num_nodes = matrix_conf.get_num_nodes();
	assert(num_in_mem >= 1);
	// If we allow 3 or more vectors in memory, we will keep
	// the out-degree vector in memory.
	dense_matrix::ptr out_deg = dense_matrix::create(get_out_degree(mat,
				num_in_mem >= 3));
	dense_matrix::ptr one = dense_matrix::create_const<float>(1,
			out_deg->get_num_rows(), out_deg->get_num_cols(),
			matrix_layout_t::L_ROW);
	out_deg = out_deg->pmax(*one);

	detail::matrix_store::ptr tmp = detail::matrix_store::create(
			mat->get_num_rows(), 1, matrix_layout_t::L_ROW,
			get_scalar_type<float>(), num_nodes, num_in_mem >= 2);
	tmp->set_data(mat_init_operate(tmp->get_num_rows()));
	dense_matrix::ptr pr1 = dense_matrix::create(tmp);

	detail::matrix_store::ptr out = detail::matrix_store::create(
			mat->get_num_rows(), 1, matrix_layout_t::L_ROW,
			// If we allow 2 or more vectors in memory, we will keep
			// the out-degree vector in memory.
			get_scalar_type<float>(), num_nodes, num_in_mem >= 2);

	size_t converge = 0;
	size_t num_iters = 0;
	const size_t num_vertices = mat->get_num_rows();
	mat = mat->transpose();
	printf("after transpose\n");
	const float eps = 0.001 / num_vertices;
	struct timeval start, end;
	while (converge < num_vertices && num_iters < max_niters) {
		struct timeval start1, end1;
		gettimeofday(&start, NULL);
		start1 = start;
		dense_matrix::ptr in = pr1->div(*out_deg);
		in = in->cast_ele_type(get_scalar_type<float>());
		// This is guaranteed to be in memory.
		in->move_store(true, num_nodes);
		pr1 = in->multiply_ele(*out_deg);
		gettimeofday(&end1, NULL);
		printf("generate input takes %.3f seconds\n", time_diff(start1, end1));
		start1 = end1;
		mat->multiply<float, float>(in->get_raw_store(), out);
		gettimeofday(&end1, NULL);
		printf("SpMM takes %.3f seconds\n", time_diff(start1, end1));
		start1 = end1;
		// pr2 is a virtual matrix, all computation below has to be computed
		// twice. Once to test the converged vertices, once to compute pr1.
		// TODO we may instruct to materialize pr2 directly.
		dense_matrix::ptr pr2 = dense_matrix::create(out);
		pr2 = pr2->multiply_scalar<float>(damping_factor);
		pr2 = pr2->add_scalar<float>((1.0 - damping_factor) / num_vertices);
		dense_matrix::ptr diff = pr1->minus(*pr2)->abs();
		gettimeofday(&end1, NULL);
		printf("virtual computation takes %.3f seconds\n", time_diff(start1, end1));
		scalar_variable::ptr convg = diff->lt_scalar<float>(eps)->sum();
		assert(convg->get_type() == get_scalar_type<size_t>());
		converge = *(const size_t *) convg->get_raw();
		pr1 = pr2;
		num_iters++;
		gettimeofday(&end, NULL);
		printf("computing converge takes %.3f seconds\n", time_diff(start1, end));
		printf("iter %ld takes %.3f seconds, not converged: %ld\n",
				num_iters, time_diff(start, end), num_vertices - converge);
	}

	return pr1;
}

}

}
